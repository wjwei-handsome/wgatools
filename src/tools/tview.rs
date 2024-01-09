use crate::tools::index::MafIndex;
use crate::{errors::WGAError, parser::maf::MAFReader};
use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use itertools::enumerate;
use ratatui::{prelude::*, widgets::*};
use regex::Regex;
use rust_lapper::{Interval, Lapper};
use std::{
    fs::File,
    io::{self, BufReader, Read, Seek},
    rc::Rc,
    time::{Duration, Instant},
};

const WINDOW_SIZE: usize = 20;

// ref to <https://en.wikipedia.org/wiki/ANSI_escape_code#8-bit>
const OPTION_8BIT_COLOR: [u8; 10] = [2, 14, 3, 4, 5, 1, 8, 27, 99, 36];

#[derive(Default)]
struct Scroll<'a> {
    scroll: usize,
    scroll_state: ScrollbarState,
    para_lines: Vec<Line<'a>>,
    ref_name: String,
    ref_start: u64,
    destpos: u64,
    seek: u64,
    block_size: usize,
}

impl Scroll<'_> {
    fn scroll_left(&mut self, step: usize) {
        self.scroll = self.scroll.saturating_sub(step);
        self.scroll_state = self.scroll_state.position(self.scroll);
    }

    fn scroll_right(&mut self, step: usize) {
        if self.scroll + step > self.block_size {
            self.scroll = self.block_size;
        } else {
            self.scroll = self.scroll.saturating_add(step);
        }
        // self.scroll = self.scroll.saturating_add(step);
        self.scroll_state = self.scroll_state.position(self.scroll);
    }

    fn scroll_init(&mut self) {
        self.scroll = 0;
        self.scroll_state = self.scroll_state.position(self.scroll);
    }
}

type Iv = Interval<u64, u64>;

struct Navigation {
    show: bool,
    input: String,
    cursor_position: usize,
    input_valid: bool,
    cddt_name: Vec<String>,
    select_name_idx: usize,
    cddt_region: Vec<String>,
    select_region_idx: usize,
    all_regions: Vec<Vec<Iv>>,
    select_region: bool,
}

impl Navigation {
    fn select_up(&mut self) {
        if self.show {
            if !self.select_region {
                if self.select_name_idx == 0 {
                    self.select_name_idx = self.cddt_name.len();
                }
                self.select_name_idx = self.select_name_idx.saturating_sub(1);
                self.cddt_region = ivvec2strvec(&self.all_regions[self.select_name_idx]);
                self.select_region_idx = 0;
                self.update_input();
            } else {
                if self.select_region_idx == 0 {
                    self.select_region_idx = self.cddt_region.len();
                }
                self.select_region_idx = self.select_region_idx.saturating_sub(1);
                self.update_input();
            }
        }
    }

    fn select_down(&mut self) {
        if self.show {
            if !self.select_region {
                if self.select_name_idx + 1 == self.cddt_name.len() {
                    self.select_name_idx = 0;
                } else {
                    self.select_name_idx = self.select_name_idx.saturating_add(1);
                }
                // update regions
                self.cddt_region = ivvec2strvec(&self.all_regions[self.select_name_idx]);
                self.select_region_idx = 0;
                self.update_input();
            } else {
                if self.select_region_idx + 1 == self.cddt_region.len() {
                    self.select_region_idx = 0;
                } else {
                    self.select_region_idx = self.select_region_idx.saturating_add(1);
                }
                self.update_input();
            }
        }
    }

    fn update_input(&mut self) {
        let name = &self.cddt_name[self.select_name_idx];
        let region = &self.cddt_region[self.select_region_idx];
        let replace_text = format!("{}:{}", name, region);
        self.input.replace_range(6.., &replace_text);
        self.cursor_position = self.input.len();
    }

    fn move_cursor_left(&mut self) {
        let cursor_moved_left = self.cursor_position.saturating_sub(1);
        self.cursor_position = self.clamp_cursor(cursor_moved_left);
    }

    fn move_cursor_right(&mut self) {
        let cursor_moved_right = self.cursor_position.saturating_add(1);
        self.cursor_position = self.clamp_cursor(cursor_moved_right);
    }

    fn enter_char(&mut self, new_char: char) {
        self.input.insert(self.cursor_position, new_char);

        self.move_cursor_right();
    }

    fn delete_char(&mut self) {
        let is_not_cursor_leftmost = self.cursor_position != 6;
        if is_not_cursor_leftmost {
            let current_index = self.cursor_position;
            let from_left_to_current_index = current_index - 1;
            // Getting all characters before the selected character.
            let before_char_to_delete = self.input.chars().take(from_left_to_current_index);
            // Getting all characters after selected character.
            let after_char_to_delete = self.input.chars().skip(current_index);
            // Put all characters together except the selected one.
            // By leaving the selected one out, it is forgotten and therefore deleted.
            self.input = before_char_to_delete.chain(after_char_to_delete).collect();
            self.move_cursor_left();
        }
    }

    fn clamp_cursor(&self, new_cursor_pos: usize) -> usize {
        new_cursor_pos.clamp(6, self.input.len())
    }
}

struct MafViewApp<'a, R: Read + Send + Seek> {
    fixed: Vec<Line<'a>>,
    scroll: Scroll<'a>,
    navigation: Navigation,
    #[allow(dead_code)]
    wait: bool,
    filerdr: MAFReader<R>,
}

impl MafViewApp<'_, File> {
    fn gen_navigation(mafindex: MafIndex) -> Navigation {
        let mut all_regions = Vec::new();
        let mut cddt_names = Vec::new();
        for (name, mafindex_item) in mafindex {
            let mut region = Vec::new();
            for ivp in &mafindex_item.ivls {
                region.push(Iv {
                    start: ivp.start,
                    stop: ivp.end,
                    val: ivp.offset,
                });
            }
            cddt_names.push(name.to_string());
            all_regions.push(region);
        }
        let cddt_regions = &all_regions[0];
        Navigation {
            show: false,
            input: "Goto: ".to_string(),
            cursor_position: 6,
            input_valid: true,
            cddt_name: cddt_names,
            select_name_idx: 0,
            cddt_region: ivvec2strvec(cddt_regions),
            select_region_idx: 0,
            all_regions,
            select_region: false,
        }
    }

    fn new(input: &String) -> Result<Self, WGAError> {
        // creat reader
        let mut mafreader = MAFReader::from_path(input)?;
        // init scroll, fixed
        let mut scroll = Scroll::default();
        let mut fixed = vec![Line::from("pos:"), Line::from("|")];
        // read index
        let index_file = File::open(format!("{}.index", input))?;
        let mafindex: MafIndex = serde_json::from_reader(BufReader::new(index_file))?;
        // create navigation
        let mut navigation = Self::gen_navigation(mafindex);

        // init first record
        let init_maf_rec = mafreader.records().next().ok_or(WGAError::EmptyRecord)??;

        // init first line as ref line
        let init_sline = &init_maf_rec.slines[0];
        let init_ref_seq = &init_sline.name;
        scroll.ref_name = init_ref_seq.to_string();
        scroll.ref_start = init_maf_rec.slines[0].start;

        let ref_seq = &init_sline.seq;

        let (axis_text, indicator_text, len_count) =
            get_axis_idc_len(ref_seq, scroll.ref_start, WINDOW_SIZE);

        let mut para_lines = vec![
            Line::from(axis_text.red()),
            Line::from(indicator_text.yellow()),
        ];
        for (idx, sline) in enumerate(&init_maf_rec.slines) {
            let seq = &sline.seq;
            let name = &sline.name;
            let color = Color::Indexed(OPTION_8BIT_COLOR[idx % 10]);
            para_lines.push(Line::from(seq.to_string().fg(color)));
            fixed.push(Line::from(name.to_string().fg(color)));
        }

        scroll.block_size = len_count;
        scroll.scroll_state = scroll.scroll_state.content_length(len_count);
        scroll.para_lines = para_lines;
        navigation.update_input();

        let app = Self {
            fixed,
            scroll,
            navigation,
            wait: false,
            filerdr: mafreader,
        };

        Ok(app)
    }

    fn update(&mut self) -> Result<(), WGAError> {
        self.filerdr
            .inner
            .seek(std::io::SeekFrom::Start(self.scroll.seek))?;
        // new mafrec
        let mafrec = self
            .filerdr
            .records()
            .next()
            .ok_or(WGAError::EmptyRecord)??;
        // init scroll
        self.scroll.scroll_init();
        // change ref line
        let mut add_para_lines = Vec::new();
        let mut add_fixed_lines = Vec::new();
        // let mut first_3_para_lines = Vec::new();
        for (idx, sline) in enumerate(mafrec.slines) {
            let name = &sline.name;
            let seq = &sline.seq;
            let option_colors = OPTION_8BIT_COLOR;
            let first_color = Color::Indexed(option_colors[0]);
            let rest_option_color = option_colors.split_at(1).1;
            let color = Color::Indexed(rest_option_color[idx % 10]);
            if self.scroll.ref_name == *name {
                // change ...
                self.scroll.ref_start = sline.start;
                let ref_start = sline.start;
                let ref_seq = &sline.seq;

                let (axis_text, indicator_text, len_count) =
                    get_axis_idc_len(ref_seq, ref_start, WINDOW_SIZE);
                let first_3_para_lines = vec![
                    Line::from(axis_text.red()),
                    Line::from(indicator_text.yellow()),
                    Line::from(seq.to_string().fg(first_color)),
                ];
                let first_3_fixed_lines = vec![
                    Line::from("pos:"),
                    Line::from("|"),
                    Line::from(name.to_string().fg(first_color)),
                ];

                self.scroll.para_lines = first_3_para_lines;
                self.fixed = first_3_fixed_lines;
                self.scroll.block_size = len_count;
                self.scroll.scroll_state = self.scroll.scroll_state.content_length(len_count);
            } else {
                add_para_lines.push(Line::from(seq.to_string().fg(color)));
                add_fixed_lines.push(Line::from(name.to_string().fg(color)));
            }
        }
        self.scroll.para_lines.append(&mut add_para_lines);
        self.fixed.append(&mut add_fixed_lines);
        // scroll
        let current_pos = self.scroll.ref_start + self.scroll.scroll as u64;
        let scroll_size = self.scroll.destpos - current_pos;
        self.scroll.scroll_right(scroll_size as usize);
        self.navigation.show = false;
        Ok(())
    }
}

pub fn tview(input: &String, step: usize) -> Result<(), WGAError> {
    // creat app and fill init data
    let app = MafViewApp::new(input)?;

    // setup terminal
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    // run app
    let tick_rate = Duration::from_millis(250);
    run_app(&mut terminal, app, tick_rate, step)?;

    // restore terminal
    disable_raw_mode()?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    terminal.show_cursor()?;

    Ok(())
}

fn run_app<B: Backend>(
    terminal: &mut Terminal<B>,
    mut app: MafViewApp<'_, File>,
    tick_rate: Duration,
    step: usize,
) -> Result<(), WGAError> {
    let mut last_tick = Instant::now();
    loop {
        terminal.draw(|f| main_ui(f, &mut app))?;

        let timeout = tick_rate.saturating_sub(last_tick.elapsed());
        if crossterm::event::poll(timeout)? {
            if let Event::Key(key) = event::read()? {
                match key.code {
                    KeyCode::Left => {
                        if app.navigation.show {
                            app.navigation.move_cursor_left();
                        } else {
                            app.scroll.scroll_left(step);
                        }
                    }
                    KeyCode::Right => {
                        if app.navigation.show {
                            app.navigation.move_cursor_right();
                        } else {
                            app.scroll.scroll_right(step);
                        }
                    }
                    KeyCode::Up => {
                        app.navigation.select_up();
                    }
                    KeyCode::Down => {
                        app.navigation.select_down();
                    }
                    KeyCode::Esc => {
                        if app.navigation.show {
                            app.navigation.show = false;
                        }
                    }
                    KeyCode::Char(input_char) => {
                        if app.navigation.show {
                            app.navigation.enter_char(input_char);
                        } else if input_char == 'q' {
                            return Ok(());
                        } else if input_char == 'g' {
                            app.navigation.show = true;
                        }
                    }
                    KeyCode::Backspace => {
                        if app.navigation.show {
                            app.navigation.delete_char();
                        }
                    }
                    KeyCode::Tab => {
                        if app.navigation.show {
                            app.navigation.select_region = !app.navigation.select_region;
                        }
                    }
                    KeyCode::Enter => {
                        if app.navigation.show && app.navigation.input_valid {
                            app.update()?;
                        }
                    }
                    _ => {}
                }
            }
        }
        if last_tick.elapsed() >= tick_rate {
            last_tick = Instant::now();
        }
    }
}

fn main_ui(f: &mut Frame, app: &mut MafViewApp<'_, File>) {
    let size = f.size();

    let block = Block::default().black();
    f.render_widget(block, size);

    let main_layout = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(15), Constraint::Percentage(85)])
        .split(size);

    let create_block = |title| {
        Block::default()
            .borders(Borders::ALL)
            .cyan()
            .title(block::Title::from(title).alignment(Alignment::Center))
            .border_type(BorderType::Rounded)
    };

    let seqname_para = Paragraph::new(app.fixed.clone()).block(create_block("seq name"));
    f.render_widget(seqname_para, main_layout[0]);

    let paragraph = Paragraph::new(app.scroll.para_lines.clone())
        .block(create_block("Press ◄ ► to scroll"))
        .scroll((0, app.scroll.scroll as u16));
    f.render_widget(paragraph, main_layout[1]);
    f.render_stateful_widget(
        Scrollbar::default()
            .orientation(ScrollbarOrientation::HorizontalBottom)
            .thumb_symbol("░")
            .track_symbol(Some("─")),
        main_layout[1].inner(&Margin {
            vertical: 0,
            horizontal: 1,
        }),
        &mut app.scroll.scroll_state,
    );

    if app.navigation.show {
        // judge if input is valid
        app.navigation.input_valid = true;
        // update jump to seek/ Destination / ref name from **valid**-input
        let _ = input_valid_update(app);

        // gen 4 split areas
        let two_scroll_one_input_one_msg = gen_two_scroll_one_input_one_msg_area(f);
        // two scroll areas
        let two_scroll = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(30), Constraint::Percentage(70)])
            .split(two_scroll_one_input_one_msg[0]);
        let scroll_area_0 = two_scroll[0];
        let scroll_area_1 = two_scroll[1];
        // input area
        let input_area = two_scroll_one_input_one_msg[1];
        // msg area
        let msg_area = two_scroll_one_input_one_msg[2];

        // set selected background color
        let (name_bg_col, region_bg_col) = if app.navigation.select_region {
            (Color::LightYellow, Color::LightGreen)
        } else {
            (Color::LightGreen, Color::LightYellow)
        };

        // fill name candidate scroll
        gen_fill_scroll(
            f,
            scroll_area_0,
            &app.navigation.cddt_name,
            app.navigation.select_name_idx,
            name_bg_col,
            "Name",
        );
        // fill region candidate scroll
        gen_fill_scroll(
            f,
            scroll_area_1,
            &app.navigation.cddt_region,
            app.navigation.select_region_idx,
            region_bg_col,
            "Region",
        );

        let input = Paragraph::new(app.navigation.input.as_str()).magenta();
        f.render_widget(input, input_area);
        f.set_cursor(
            input_area.x + app.navigation.cursor_position as u16,
            input_area.y,
        );

        let message = if app.navigation.input_valid {
            "Press ▲ ▼ to select, <Tab> to switch between name and region, <Esc> to exit, <Enter> to jump"
        } else {
            "Invalid input, please re-select or enter"
        };
        let msg = Paragraph::new(message).light_blue();
        f.render_widget(msg, msg_area);
    }
}

fn get_axis_idc_len(seq: &str, start: u64, window_size: usize) -> (String, String, usize) {
    let start = start + 1; // MAF is 0-based
    let mut axis_text = String::new();
    let mut indicator_text = String::new();
    let mut idx = 0;
    let mut len_count = 0;
    for base in seq.chars() {
        len_count += 1;
        if base == '-' {
            axis_text.push(' ');
            indicator_text.push(' ');
        } else if idx % window_size == 0 {
            let pos = start + idx as u64;
            axis_text.push_str(&format!(
                "{:width$}",
                pos.to_string(),
                width = { window_size }
            ));
            indicator_text.push('|');
            idx += 1
        } else {
            indicator_text.push(' ');
            idx += 1
        }
    }
    axis_text.push('\n');
    indicator_text.push('\n');
    (axis_text, indicator_text, len_count)
}

fn ivvec2strvec(invec: &[Iv]) -> Vec<String> {
    invec
        .iter()
        .map(|i| format!("{}-{}", i.start, i.stop))
        .collect::<Vec<String>>()
}

fn input_valid_update(app: &mut MafViewApp<'_, File>) -> Result<(), WGAError> {
    let re = Regex::new(r"^[a-zA-Z0-9.@-_#]+:[0-9]+-[0-9]+?$")?; // NO ERROR
    match re.is_match(&app.navigation.input[6..]) {
        true => {
            let name = &app.navigation.input[6..].split(':').collect::<Vec<&str>>()[0];
            match app.navigation.cddt_name.iter().position(|i| i == name) {
                Some(name_idx) => {
                    let region = &app.navigation.input[6..].split(':').collect::<Vec<&str>>()[1];
                    let start = match region.split('-').collect::<Vec<&str>>()[0].parse::<usize>() {
                        Ok(i) => i,
                        Err(_) => {
                            app.navigation.input_valid = false;
                            0
                        }
                    };
                    let end = match region.contains('-') {
                        true => {
                            match region.split('-').collect::<Vec<&str>>()[1].parse::<usize>() {
                                Ok(i) => i,
                                Err(_) => {
                                    app.navigation.input_valid = false;
                                    0
                                }
                            }
                        }
                        false => start,
                    };
                    if start > end {
                        app.navigation.input_valid = false;
                    } else {
                        let cddt_regions: &Vec<Iv> = &app.navigation.all_regions[name_idx];
                        let lapper = Lapper::new(cddt_regions.clone());
                        let find = lapper
                            .find(start as u64, start as u64 + 1)
                            .collect::<Vec<&Iv>>();
                        if find.is_empty() {
                            app.navigation.input_valid = false;
                        } else {
                            let dest_block = find[0];
                            app.scroll.seek = dest_block.val;
                            app.scroll.destpos = start as u64;
                            app.scroll.ref_name = name.to_string();
                        }
                    }
                }
                None => {
                    app.navigation.input_valid = false;
                }
            }
        }
        false => {
            app.navigation.input_valid = false;
        }
    }
    Ok(())
}

fn gen_two_scroll_one_input_one_msg_area(f: &mut Frame) -> Rc<[Rect]> {
    let percent_y = 20;
    let percent_x = 60;
    let ver_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(f.size());

    let hor_layout = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(ver_layout[1]);
    Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage(80),
            Constraint::Percentage(10),
            Constraint::Percentage(10),
        ])
        .split(hor_layout[1])
}

fn gen_fill_scroll(
    f: &mut Frame,
    area: Rect,
    cddt: &Vec<String>,
    idx: usize,
    select_bg_col: Color,
    title: &str,
) {
    let scroll_item = cddt
        .iter()
        .map(|i| {
            ListItem::new(
                Text::from(i.clone()),
                // lines.into_iter(),
            )
            .style(Style::default().fg(Color::Black).bg(Color::White))
        })
        .collect::<Vec<ListItem>>();
    let scroll_item = List::new(scroll_item)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(block::Title::from(title).alignment(Alignment::Center))
                .border_type(BorderType::Rounded),
        )
        .highlight_style(
            Style::default()
                .bg(select_bg_col)
                .add_modifier(Modifier::BOLD)
                .add_modifier(Modifier::ITALIC),
        )
        .highlight_symbol(">> ");
    let mut scroll_state = ListState::default().with_selected(Some(idx));
    f.render_stateful_widget(scroll_item, area, &mut scroll_state);
    let mut scrollbar_state = ScrollbarState::default()
        .content_length(cddt.len())
        .position(idx);
    f.render_stateful_widget(
        Scrollbar::default()
            .thumb_symbol("░")
            .track_symbol(Some("─")),
        area,
        &mut scrollbar_state,
    );
}
