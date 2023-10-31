use std::{
    error::Error,
    io,
    rc::Rc,
    time::{Duration, Instant},
};

use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{prelude::*, widgets::*};

use crate::parser::{common::AlignRecord, maf::MAFReader};

struct Scroll<'a> {
    axis_scroll_state: ScrollbarState,
    axis_scroll: usize,
    indicator_scroll_state: ScrollbarState,
    indicator_scroll: usize,
    sline_scroll_state: ScrollbarState,
    sline_scroll: usize,
    para_lines: Vec<Line<'a>>,
}

impl Scroll<'_> {
    fn scroll_left(&mut self, step: usize) {
        self.axis_scroll = self.axis_scroll.saturating_sub(step);
        self.axis_scroll_state = self.axis_scroll_state.position(self.axis_scroll);
        self.indicator_scroll = self.indicator_scroll.saturating_sub(step);
        self.indicator_scroll_state = self.indicator_scroll_state.position(self.indicator_scroll);
        self.sline_scroll = self.sline_scroll.saturating_sub(step);
        self.sline_scroll_state = self.sline_scroll_state.position(self.sline_scroll);
    }

    fn scroll_right(&mut self, step: usize) {
        self.axis_scroll = self.axis_scroll.saturating_add(step);
        self.axis_scroll_state = self.axis_scroll_state.position(self.axis_scroll);
        self.indicator_scroll = self.indicator_scroll.saturating_add(step);
        self.indicator_scroll_state = self.indicator_scroll_state.position(self.indicator_scroll);
        self.sline_scroll = self.sline_scroll.saturating_add(step);
        self.sline_scroll_state = self.sline_scroll_state.position(self.sline_scroll);
    }
}

struct Navigation {
    show: bool,
    all_seqs: Vec<String>,
    candidate_seqs: Vec<String>,
    select_index: usize,
    input: String,
    cursor_position: usize,
}

impl Navigation {
    fn select_up(&mut self) {
        if self.show {
            if self.select_index == 0 {
                self.select_index = self.candidate_seqs.len();
            }
            self.select_index = self.select_index.saturating_sub(1);
            self.input
                .replace_range(6.., &self.candidate_seqs[self.select_index]);
            self.input.push(':');
            self.cursor_position = self.input.len();
        }
    }

    fn select_down(&mut self) {
        if self.show {
            if self.select_index + 1 == self.candidate_seqs.len() {
                self.select_index = 0;
                return;
            }
            self.select_index = self.select_index.saturating_add(1);
            self.input
                .replace_range(6.., &self.candidate_seqs[self.select_index]);
            self.input.push(':');
            self.cursor_position = self.input.len();
        }
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

struct MafViewApp<'a> {
    fixed: Vec<Line<'a>>,
    indexed: bool,
    scroll: Scroll<'a>,
    navigation: Navigation,
}

impl Default for MafViewApp<'_> {
    fn default() -> Self {
        let item_list = vec!["Zm-B73v5.chr8".to_string(), "Zm-CML333.chr8".to_string()];
        let scroll = Scroll {
            axis_scroll_state: ScrollbarState::default(),
            axis_scroll: 0,
            indicator_scroll_state: ScrollbarState::default(),
            indicator_scroll: 0,
            sline_scroll_state: ScrollbarState::default(),
            sline_scroll: 0,
            para_lines: Vec::new(),
        };
        let navigation = Navigation {
            show: false,
            all_seqs: item_list.clone(),
            candidate_seqs: item_list,
            select_index: 0,
            input: "Goto: ".to_string(),
            cursor_position: 6,
        };
        let fixed = vec![
            Line::from("pos: "),
            Line::from("|"),
            Line::from("target"),
            Line::from("query"),
        ];
        Self {
            fixed,
            indexed: false,
            scroll,
            navigation,
        }
    }
}

pub fn test_tview() -> Result<(), Box<dyn Error>> {
    // setup terminal
    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableMouseCapture)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;

    // create app and run it
    let tick_rate = Duration::from_millis(250);
    let mut app = MafViewApp::default();
    app.navigation.candidate_seqs = app.navigation.all_seqs.clone();

    // impute fake data
    let mut mafreader = MAFReader::from_path("sub.maf").unwrap();
    let mafrec = mafreader.records().next().unwrap().unwrap();

    let mut target_seq = mafrec.target_seq().to_string();
    target_seq.push('\n');
    let mut query_seq = mafrec.query_seq().to_string();
    query_seq.push('\n');

    let t_start = mafrec.target_start();
    let t_end = mafrec.target_end();

    let t_name = mafrec.target_name();
    let q_name = mafrec.query_name();

    let mut axis_text = String::new();
    let mut indicator_text = String::new();
    let window_size = 20;

    for pos in (t_start..t_end).step_by(window_size) {
        axis_text.push_str(&format!("{}", pos));
        indicator_text.push('|');
        indicator_text.push_str(&" ".repeat(window_size - 1));
        // count pos digits
        let pos_digits = pos.to_string().len();
        let fill_size = window_size - pos_digits;
        axis_text.push_str(&" ".repeat(fill_size));
    }

    axis_text.push('\n');
    indicator_text.push('\n');

    app.scroll.axis_scroll_state = app.scroll.axis_scroll_state.content_length(axis_text.len());
    app.scroll.indicator_scroll_state = app
        .scroll
        .indicator_scroll_state
        .content_length(indicator_text.len());
    app.scroll.sline_scroll_state = app
        .scroll
        .sline_scroll_state
        .content_length(target_seq.len());
    let text = vec![
        Line::from(axis_text.red()),
        Line::from(indicator_text.yellow()),
        Line::from(target_seq.green()),
        Line::from(query_seq.blue()),
    ];
    app.scroll.para_lines = text;

    app.fixed[2] = Line::from(t_name.green());
    app.fixed[3] = Line::from(q_name.blue());

    let res = run_app(&mut terminal, app, tick_rate);

    // restore terminal
    disable_raw_mode()?;
    execute!(
        terminal.backend_mut(),
        LeaveAlternateScreen,
        DisableMouseCapture
    )?;
    terminal.show_cursor()?;

    if let Err(err) = res {
        println!("{err:?}");
    }

    Ok(())
}

fn run_app<B: Backend>(
    terminal: &mut Terminal<B>,
    mut app: MafViewApp,
    tick_rate: Duration,
) -> io::Result<()> {
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
                            app.scroll.scroll_left(10);
                        }
                    }
                    KeyCode::Right => {
                        if app.navigation.show {
                            app.navigation.move_cursor_right();
                        } else {
                            app.scroll.scroll_right(10);
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
                    _ => {}
                }
            }
        }
        if last_tick.elapsed() >= tick_rate {
            last_tick = Instant::now();
        }
    }
}

fn main_ui(f: &mut Frame, app: &mut MafViewApp) {
    let size = f.size();

    let block = Block::default().black();
    f.render_widget(block, size);

    let main_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Min(1), Constraint::Percentage(99)])
        .split(size);

    let inner_layout = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(10), Constraint::Percentage(90)])
        .split(main_layout[1]);

    let create_block = |title| {
        Block::default()
            .borders(Borders::ALL)
            .gray()
            .title(Span::styled(
                title,
                Style::default().add_modifier(Modifier::BOLD),
            ))
    };

    let title = Block::default()
        .title("Use ◄ ► to scroll")
        .title_alignment(Alignment::Center);
    f.render_widget(title, main_layout[0]);

    let seqname_para = Paragraph::new(app.fixed.clone())
        .gray()
        .block(create_block("seq name"));
    f.render_widget(seqname_para, inner_layout[0]);

    let paragraph = Paragraph::new(app.scroll.para_lines.clone())
        .gray()
        .block(create_block(
            "Horizontal scrollbar with only begin arrow & custom thumb symbol",
        ))
        .scroll((0, app.scroll.axis_scroll as u16));
    f.render_widget(paragraph, inner_layout[1]);
    f.render_stateful_widget(
        Scrollbar::default()
            .orientation(ScrollbarOrientation::HorizontalBottom)
            .thumb_symbol("*")
            .end_symbol(None),
        inner_layout[1].inner(&Margin {
            vertical: 0,
            horizontal: 1,
        }),
        &mut app.scroll.axis_scroll_state,
    );

    if app.navigation.show {
        let mut matched = true;

        let area = centered_rect(60, 20, size);
        f.render_widget(Clear, area[0]);
        f.render_widget(Clear, area[1]);

        let match_text = &app.navigation.input[6..];
        if !match_text.ends_with(':') {
            let mut new_candidate_seqs = Vec::new();
            for item in app.navigation.all_seqs.iter() {
                if item.contains(match_text) {
                    new_candidate_seqs.push(item.clone());
                }
            }
            if !new_candidate_seqs.is_empty() {
                app.navigation.candidate_seqs = new_candidate_seqs;
                app.navigation.select_index = 0;
            } else {
                matched = false;
            }
        }
        if matched {
            let items: Vec<ListItem> = app
                .navigation
                .candidate_seqs
                .iter()
                .map(|i| {
                    let lines = vec![i.clone().italic().into()];
                    ListItem::new(lines).style(Style::default().fg(Color::Black).bg(Color::White))
                })
                .collect();
            // Create a List from all list items and highlight the currently selected one
            let items = List::new(items)
                .block(Block::default().borders(Borders::ALL).title("List"))
                .highlight_style(
                    Style::default()
                        .bg(Color::LightGreen)
                        .add_modifier(Modifier::BOLD),
                )
                .highlight_symbol(">> ");

            let mut state = ListState::default().with_selected(Some(app.navigation.select_index));

            f.render_stateful_widget(items, area[0], &mut state);

            let mut scrollbar_state = ScrollbarState::default()
                .content_length(app.navigation.candidate_seqs.len())
                .position(app.navigation.select_index);

            f.render_stateful_widget(
                Scrollbar::default()
                    // .orientation(ScrollbarOrientation::VerticalRight)
                    .thumb_symbol("░")
                    .track_symbol(Some("─")),
                area[0],
                &mut scrollbar_state,
            );
        } else {
            let text = vec![Line::from("No match!")];
            let paragraph = Paragraph::new(text).gray();
            f.render_widget(paragraph, area[0]);
        }

        let input = Paragraph::new(app.navigation.input.as_str()).gray();
        f.render_widget(input, area[1]);
        f.set_cursor(area[1].x + app.navigation.cursor_position as u16, area[1].y);
    }
}

/// helper function to create a centered rect using up certain percentage of the available rect `r`
fn centered_rect(percent_x: u16, percent_y: u16, r: Rect) -> Rc<[Rect]> {
    let ver_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(r);

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
        .constraints([Constraint::Percentage(80), Constraint::Percentage(20)])
        .split(hor_layout[1])
}
