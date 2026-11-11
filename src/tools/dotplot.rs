use crate::{
    errors::WGAError,
    parser::{
        cigar::{parse_cigar_to_base_plotdata, parse_maf_to_base_plotdata},
        common::{AlignRecord, DotplotMode, DotplotoutFormat, FileFormat, Strand},
        maf::MAFReader,
        paf::PAFReader,
    },
};
use minijinja::{context, Environment};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::HashMap,
    io::{BufRead, Read, Write},
};

fn parse_color_config(color_str: &str) -> Result<HashMap<String, String>, WGAError> {
    let mut color_map = HashMap::new();
    for pair in color_str.split(',') {
        let parts: Vec<&str> = pair.split(':').collect();
        if parts.len() != 2 {
            return Err(WGAError::Other(anyhow::anyhow!(
                "Invalid color format: {}. Expected format: M:#FF0000",
                pair
            )));
        }
        // Validate hex color code
        if !parts[1].starts_with('#') || parts[1].len() != 7 {
            return Err(WGAError::Other(anyhow::anyhow!(
                "Invalid hex color code: {}. Expected format: #RRGGBB",
                parts[1]
            )));
        }
        color_map.insert(parts[0].to_string(), parts[1].to_string());
    }
    Ok(color_map)
}

/// value of the chromosome selectors meaning "show the whole genome"
const ALL_CHROMS: &str = "(all)";
/// default plot size (px); the HTML page resizes it to fit the window
const DEFAULT_SIZE: u64 = 800;
/// a chromosome is labelled only if it takes at least this fraction of the genome
const LABEL_MIN_FRACTION: f64 = 0.005;

const VEGA_TEMP: &str = r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>wgatools dotplot</title>
<script src="https://cdn.jsdelivr.net/npm/vega@5.30.0"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-lite@5.21.0"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-embed@6.26.0"></script>
<style>
  :root { --fg: #1f2933; --muted: #616e7c; --line: #e4e7eb; --accent: #1d4e89; }
  * { box-sizing: border-box; }
  body { margin: 0; padding: 16px 24px; font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; color: var(--fg); background: #fff; }
  header { display: flex; flex-wrap: wrap; align-items: baseline; gap: 4px 20px; padding-bottom: 10px; border-bottom: 1px solid var(--line); margin-bottom: 12px; }
  h1 { font-size: 18px; margin: 0; }
  .stats { display: flex; flex-wrap: wrap; gap: 4px 16px; color: var(--muted); }
  .stats b { color: var(--fg); font-weight: 600; }
  .help { color: var(--muted); font-size: 13px; margin: 0 0 8px; }
  .help kbd { font: 12px ui-monospace, monospace; background: #f5f7fa; border: 1px solid var(--line); border-radius: 3px; padding: 0 4px; }
  #view .vega-bindings { display: flex; flex-wrap: wrap; gap: 8px 24px; margin-top: 10px; padding: 10px 12px; background: #f5f7fa; border-radius: 6px; }
  #view .vega-bind { display: flex; align-items: center; }
  #view .vega-bind-name { font-weight: 600; margin-right: 6px; }
  #view .vega-bind select { max-width: 220px; padding: 2px 4px; }
  #view .vega-bind label { margin-right: 8px; }
  #view .vega-bind input[type=range] { width: 180px; accent-color: var(--accent); }
  #error { color: #b42318; white-space: pre-wrap; }
</style>
</head>
<body>
<header>
  <h1>wgatools dotplot</h1>
  <div class="stats">
    <span>Mode <b>{{ mode }}</b></span>
    <span><b>{{ n_records }}</b> {{ record_kind }}</span>
    <span><b>{{ n_ref }}</b> reference chromosomes</span>
    <span><b>{{ n_query }}</b> query chromosomes</span>
  </div>
</header>
<p class="help">
  <kbd>scroll</kbd> zoom · <kbd>drag</kbd> pan · <kbd>double-click</kbd> reset view · hover a line for details ·
  use the controls below the plot to focus on one chromosome or filter alignments{{ extra_help }} ·
  <kbd>…</kbd> menu (top right) saves PNG/SVG
</p>
<div id="view"></div>
<div id="error"></div>
<script>
  const spec = {{ vl_json | safe }};
  // fit the (square) plot to the window
  const size = Math.max(320, Math.min(window.innerWidth - 220, window.innerHeight - 240, 1400));
  spec.width = size;
  spec.height = size;
  vegaEmbed('#view', spec, { renderer: 'canvas' }).then(res => { window.vegaView = res.view; }).catch(err => {
    console.error(err);
    document.getElementById('error').textContent = 'Failed to render the plot: ' + err;
  });
</script>
</body>
</html>
"#;

#[derive(Debug, Deserialize, Serialize)]
struct AllPlotdata {
    ref_start: u64,
    ref_end: u64,
    query_start: u64,
    query_end: u64,
    identity: f64,
    ref_chro: String,
    query_chro: String,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct BasePlotdata {
    pub ref_start: u64,
    pub ref_end: u64,
    pub query_start: u64,
    pub query_end: u64,
    pub cigar: char,
    pub ref_chro: String,
    pub query_chro: String,
}

/// chromosome lengths seen in the alignments
#[derive(Debug, Default)]
struct ChromLens {
    refs: HashMap<String, u64>,
    queries: HashMap<String, u64>,
}

impl ChromLens {
    fn add<T: AlignRecord>(&mut self, rec: &T) {
        if !self.refs.contains_key(rec.target_name()) {
            self.refs
                .insert(rec.target_name().to_string(), rec.target_length());
        }
        if !self.queries.contains_key(rec.query_name()) {
            self.queries
                .insert(rec.query_name().to_string(), rec.query_length());
        }
    }

    fn merge(mut self, other: ChromLens) -> ChromLens {
        for (k, v) in other.refs {
            self.refs.entry(k).or_insert(v);
        }
        for (k, v) in other.queries {
            self.queries.entry(k).or_insert(v);
        }
        self
    }
}

/// chromosomes laid end to end in natural order
struct ChromLayout {
    table: Vec<Value>,
    names: Vec<String>,
    total: u64,
}

fn chrom_layout(lens: &HashMap<String, u64>) -> ChromLayout {
    let mut chroms: Vec<(&String, &u64)> = lens.iter().collect();
    chroms.sort_by(|a, b| natord::compare(a.0, b.0));
    let total: u64 = chroms.iter().map(|(_, l)| **l).sum();
    let mut offset = 0;
    let mut table = Vec::with_capacity(chroms.len());
    for (name, len) in &chroms {
        table.push(json!({
            "name": name,
            "len": len,
            "offset": offset,
            "end": offset + **len,
            "mid": offset + **len / 2,
            "label": total > 0 && (**len as f64 / total as f64) >= LABEL_MIN_FRACTION,
        }));
        offset += **len;
    }
    ChromLayout {
        table,
        names: chroms.iter().map(|(n, _)| n.to_string()).collect(),
        total,
    }
}

/// build the Vega-Lite spec: one genome-wide dotplot with chromosome
/// boundaries, chromosome selectors and filters
fn build_spec(
    lens: &ChromLens,
    max_len: u64,
    aligned_len: u64,
    mode: DotplotMode,
    min_identity: Option<f64>,
    color_map: Option<HashMap<String, String>>,
) -> Value {
    let ref_layout = chrom_layout(&lens.refs);
    let query_layout = chrom_layout(&lens.queries);
    let all = |names: &[String]| {
        let mut opts = vec![ALL_CHROMS.to_string()];
        opts.extend(names.iter().cloned());
        opts
    };
    let len_step = (max_len / 200).max(1);
    // alignments covering a tiny part of the genome would be invisible at
    // full-chromosome scale: fit the view to the alignments by default
    let fit_default =
        if ref_layout.total > 0 && (aligned_len as f64) < 0.01 * ref_layout.total as f64 {
            "alignments"
        } else {
            "chromosomes"
        };

    let mut params = vec![
        json!({
            "name": "ref_sel", "value": ALL_CHROMS,
            "bind": {"input": "select", "options": all(&ref_layout.names), "name": "Reference"}
        }),
        json!({
            "name": "query_sel", "value": ALL_CHROMS,
            "bind": {"input": "select", "options": all(&query_layout.names), "name": "Query"}
        }),
        json!({
            "name": "fit", "value": fit_default,
            "bind": {"input": "radio", "options": ["chromosomes", "alignments"], "name": "Fit view to"}
        }),
        json!({
            "name": "min_len", "value": 0,
            "bind": {"input": "range", "min": 0, "max": max_len, "step": len_step, "name": "Min length (bp)"}
        }),
    ];
    let is_overview = matches!(mode, DotplotMode::Overview);
    if is_overview && min_identity.is_some() {
        params.push(json!({
            "name": "min_identity", "value": 0,
            "bind": {"input": "range", "min": 0, "max": 1, "step": 0.01, "name": "Min identity"}
        }));
    }

    // every layer uses identical scale/axis settings so they share one view
    let pos = |field: &str, title: &str| {
        json!({
            "field": field, "type": "quantitative",
            "scale": {"zero": false, "nice": false},
            "axis": {"title": title, "format": "~s", "grid": false, "tickCount": 8}
        })
    };
    let xpos = |field: &str| pos(field, "Reference");
    let ypos = |field: &str| pos(field, "Query");

    let ref_all = format!("ref_sel == '{ALL_CHROMS}'");
    let query_all = format!("query_sel == '{ALL_CHROMS}'");

    // ---- alignment layer ----
    let mut transform = vec![
        json!({"lookup": "ref_chro", "from": {"data": {"name": "ref_chroms"}, "key": "name", "fields": ["offset"]}, "as": ["ref_off"]}),
        json!({"lookup": "query_chro", "from": {"data": {"name": "query_chroms"}, "key": "name", "fields": ["offset"]}, "as": ["query_off"]}),
        json!({"filter": format!("({ref_all} || datum.ref_chro == ref_sel) && ({query_all} || datum.query_chro == query_sel)")}),
        json!({"calculate": "abs(datum.ref_end - datum.ref_start)", "as": "ref_len"}),
        json!({"calculate": "abs(datum.query_end - datum.query_start)", "as": "query_len"}),
        json!({"filter": "max(datum.ref_len, datum.query_len) >= min_len"}),
        json!({"calculate": format!("({ref_all} ? datum.ref_off : 0) + datum.ref_start"), "as": "x"}),
        json!({"calculate": format!("({ref_all} ? datum.ref_off : 0) + datum.ref_end"), "as": "x2"}),
        json!({"calculate": format!("({query_all} ? datum.query_off : 0) + datum.query_start"), "as": "y"}),
        json!({"calculate": format!("({query_all} ? datum.query_off : 0) + datum.query_end"), "as": "y2"}),
        json!({"calculate": "datum.query_start > datum.query_end ? '-' : '+'", "as": "strand"}),
        json!({"calculate": "datum.ref_chro + ':' + format(datum.ref_start, ',') + '-' + format(datum.ref_end, ',')", "as": "ref"}),
        json!({"calculate": "datum.query_chro + ':' + format(min(datum.query_start, datum.query_end), ',') + '-' + format(max(datum.query_start, datum.query_end), ',')", "as": "query"}),
    ];
    let mut tooltip = vec![
        json!({"field": "ref", "title": "Reference"}),
        json!({"field": "query", "title": "Query"}),
        json!({"field": "strand", "title": "Strand"}),
        json!({"field": "ref_len", "type": "quantitative", "format": ",", "title": "Ref length"}),
        json!({"field": "query_len", "type": "quantitative", "format": ",", "title": "Query length"}),
    ];
    let mut layer_params = vec![json!({"name": "zoom", "select": "interval", "bind": "scales"})];
    let mut encoding = json!({
        "x": xpos("x"),
        "x2": {"field": "x2"},
        "y": ypos("y"),
        "y2": {"field": "y2"},
    });

    if is_overview {
        if let Some(min_identity) = min_identity {
            // fixed domain so colors stay stable while filtering, and a
            // constant identity is not drawn with the "low identity" color
            let domain_min = ((min_identity * 100.0).floor() / 100.0).min(0.99);
            transform.push(json!({"filter": "datum.identity >= min_identity"}));
            tooltip.push(json!({"field": "identity", "type": "quantitative", "format": ".2%", "title": "Identity"}));
            encoding["color"] = json!({
                "field": "identity", "type": "quantitative", "title": "Identity",
                "scale": {"domain": [domain_min, 1.0], "range": ["#f28e2b", "#1d4e89"], "interpolate": "hcl", "clamp": true},
                "legend": {"format": ".1%", "gradientLength": 200}
            });
        } else {
            encoding["color"] = json!({"value": "#1d4e89"});
        }
        encoding["strokeWidth"] = json!({"value": 1.5});
    } else {
        // base-level: color by CIGAR type, click legend to highlight one type
        transform.push(json!({
            "calculate": "datum.cigar == 'M' && datum.query_start > datum.query_end ? 'M_R' : datum.cigar",
            "as": "cigar"
        }));
        tooltip.push(json!({"field": "cigar", "title": "CIGAR"}));
        layer_params.push(json!({
            "name": "cigartype", "select": {"type": "point", "fields": ["cigar"]}, "bind": "legend"
        }));
        let mut scale = json!({"scheme": "category10"});
        if let Some(color_map) = color_map {
            let (domain, range): (Vec<String>, Vec<String>) = color_map.into_iter().unzip();
            scale = json!({"domain": domain, "range": range});
        }
        encoding["color"] = json!({
            "field": "cigar", "type": "nominal", "title": "CIGAR (click)", "scale": scale,
            "legend": {"symbolType": "stroke", "symbolStrokeWidth": 4}
        });
        encoding["opacity"] =
            json!({"condition": {"param": "cigartype", "value": 1}, "value": 0.15});
        encoding["strokeWidth"] =
            json!({"condition": {"param": "cigartype", "value": 2.5}, "value": 1.5});
    }
    encoding["tooltip"] = json!(tooltip);

    let alignment_layer = json!({
        "params": layer_params,
        "mark": {"type": "rule", "strokeCap": "round", "clip": true},
        "transform": transform,
        "encoding": encoding,
    });

    // ---- chromosome boundaries and labels (whole-genome view only) ----
    // in the whole-genome view: chromosome boundaries; when one chromosome is
    // selected: its start and end, so the view spans the full chromosome
    let boundary = |data: &str, sel: &str, channel: &str| {
        let filter =
            format!("fit == 'chromosomes' && ({sel} == '{ALL_CHROMS}' || datum.name == {sel})");
        let start = format!("{sel} == '{ALL_CHROMS}' ? datum.offset : 0");
        let end = format!("{sel} == '{ALL_CHROMS}' ? datum.end : datum.len");
        let def = if channel == "x" {
            xpos("pos")
        } else {
            ypos("pos")
        };
        json!({
            "data": {"name": data},
            "transform": [
                {"filter": filter},
                {"calculate": start, "as": "start"},
                {"calculate": end, "as": "stop"},
                {"fold": ["start", "stop"], "as": ["kind", "pos"]}
            ],
            "mark": {"type": "rule", "color": "#9aa5b1", "strokeWidth": 0.6, "strokeDash": [3, 3], "clip": true},
            "encoding": {channel: def}
        })
    };
    // many / long reference chromosome names: vertical labels so they don't overlap
    let ref_label_chars: usize = ref_layout
        .table
        .iter()
        .filter(|c| c["label"] == true)
        .filter_map(|c| c["name"].as_str())
        .map(|n| n.len() + 2)
        .sum();
    let ref_label_mark = if ref_label_chars > 60 {
        json!({"type": "text", "angle": 90, "align": "left", "baseline": "middle", "dy": 0, "dx": 4, "fontSize": 11, "color": "#3e4c59", "clip": true})
    } else {
        json!({"type": "text", "baseline": "top", "dy": 4, "fontSize": 11, "color": "#3e4c59", "clip": true})
    };
    let ref_labels = json!({
        "data": {"name": "ref_chroms"},
        "transform": [{"filter": format!("fit == 'chromosomes' && {ref_all} && datum.label")}],
        "mark": ref_label_mark,
        "encoding": {
            "x": xpos("mid"),
            "y": {"value": 0},
            "text": {"field": "name"},
            "tooltip": [{"field": "name", "title": "Reference"}, {"field": "len", "type": "quantitative", "format": ",", "title": "Length"}]
        }
    });
    let query_labels = json!({
        "data": {"name": "query_chroms"},
        "transform": [{"filter": format!("fit == 'chromosomes' && {query_all} && datum.label")}],
        "mark": {"type": "text", "align": "left", "baseline": "middle", "dx": 4, "fontSize": 11, "color": "#3e4c59", "clip": true},
        "encoding": {
            "y": ypos("mid"),
            "x": {"value": 0},
            "text": {"field": "name"},
            "tooltip": [{"field": "name", "title": "Query"}, {"field": "len", "type": "quantitative", "format": ",", "title": "Length"}]
        }
    });

    json!({
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "width": DEFAULT_SIZE,
        "height": DEFAULT_SIZE,
        "data": {"values": []},
        "datasets": {
            "ref_chroms": ref_layout.table,
            "query_chroms": query_layout.table,
        },
        "params": params,
        "layer": [
            boundary("ref_chroms", "ref_sel", "x"),
            boundary("query_chroms", "query_sel", "y"),
            alignment_layer,
            ref_labels,
            query_labels,
        ],
        "resolve": {"scale": {"color": "independent"}},
        "config": {
            "view": {"stroke": "#cbd2d9"},
            "axis": {"labelFontSize": 11, "titleFontSize": 13, "titleFontWeight": 600},
            "legend": {"labelFontSize": 12, "titleFontSize": 12}
        },
        "usermeta": {"ref_total": ref_layout.total, "query_total": query_layout.total}
    })
}

/// sum over reference chromosomes of the span (first start to last end)
/// covered by alignments
fn aligned_span<'a, I: Iterator<Item = (&'a str, u64, u64)>>(iter: I) -> u64 {
    let mut spans: HashMap<&str, (u64, u64)> = HashMap::new();
    for (chro, start, end) in iter {
        let (lo, hi) = (start.min(end), start.max(end));
        let span = spans.entry(chro).or_insert((lo, hi));
        span.0 = span.0.min(lo);
        span.1 = span.1.max(hi);
    }
    spans.values().map(|(lo, hi)| hi - lo).sum()
}

/// information shown in the HTML page header
struct PageInfo {
    mode: &'static str,
    record_kind: &'static str,
    extra_help: &'static str,
    n_ref: usize,
    n_query: usize,
}

#[allow(clippy::too_many_arguments)]
pub fn dotplot(
    reader: Box<dyn BufRead + Send>,
    writer: &mut dyn Write,
    format: FileFormat,
    out_format: DotplotoutFormat,
    mode: DotplotMode,
    no_identity: bool,
    skip_cutoff: usize,
    query_name: Option<&str>,
    color_config: Option<&str>,
) -> Result<(), WGAError> {
    // match mode to generate data
    match mode {
        DotplotMode::Overview => {
            let (pair_stat_vec, lens) = match format {
                FileFormat::Maf => {
                    generate_maf_data(MAFReader::new(reader)?, no_identity, query_name)?
                }
                FileFormat::Paf => generate_paf_data(PAFReader::new(reader), no_identity)?,
                _ => {
                    return Err(WGAError::Other(anyhow::anyhow!(
                        "Only support MAF and PAF format"
                    )));
                }
            };
            let max_len = pair_stat_vec
                .iter()
                .map(|d| (d.ref_end - d.ref_start).max(d.query_end.abs_diff(d.query_start)))
                .max()
                .unwrap_or(0);
            let min_identity = (!no_identity).then(|| {
                pair_stat_vec
                    .iter()
                    .map(|d| d.identity)
                    .fold(1.0_f64, f64::min)
            });
            let aligned_len = aligned_span(
                pair_stat_vec
                    .iter()
                    .map(|d| (d.ref_chro.as_str(), d.ref_start, d.ref_end)),
            );
            let spec = build_spec(&lens, max_len, aligned_len, mode, min_identity, None);
            let page = PageInfo {
                mode: "overview",
                record_kind: "alignments",
                extra_help: "",
                n_ref: lens.refs.len(),
                n_query: lens.queries.len(),
            };
            render_output(pair_stat_vec, writer, out_format, spec, &page)?;
        }
        DotplotMode::BaseLevel => {
            let (pair_base_plot_vec, lens) = match format {
                FileFormat::Maf => {
                    generate_maf_basedata(MAFReader::new(reader)?, skip_cutoff, query_name)?
                }
                FileFormat::Paf => generate_paf_basedata(PAFReader::new(reader), skip_cutoff)?,
                _ => {
                    return Err(WGAError::Other(anyhow::anyhow!(
                        "Only support MAF and PAF format"
                    )));
                }
            };
            let final_base_plotdata = pair_base_plot_vec
                .into_par_iter()
                .flatten()
                .collect::<Vec<_>>();
            let max_len = final_base_plotdata
                .iter()
                .map(|d| {
                    d.ref_end
                        .abs_diff(d.ref_start)
                        .max(d.query_end.abs_diff(d.query_start))
                })
                .max()
                .unwrap_or(0);

            let color_map = color_config.map(parse_color_config).transpose()?;
            let aligned_len = aligned_span(
                final_base_plotdata
                    .iter()
                    .map(|d| (d.ref_chro.as_str(), d.ref_start, d.ref_end)),
            );
            let spec = build_spec(&lens, max_len, aligned_len, mode, None, color_map);
            let page = PageInfo {
                mode: "base-level",
                record_kind: "segments",
                extra_help:
                    " · click a CIGAR type in the legend to highlight it (shift-click for several)",
                n_ref: lens.refs.len(),
                n_query: lens.queries.len(),
            };
            render_output(final_base_plotdata, writer, out_format, spec, &page)?;
        }
    }
    Ok(())
}

/// render data output
fn render_output<S: Serialize>(
    data: Vec<S>,
    writer: &mut dyn Write,
    format: DotplotoutFormat,
    mut vega_spec: Value,
    page: &PageInfo,
) -> Result<(), WGAError> {
    match format {
        DotplotoutFormat::Json => {
            vega_spec["data"]["values"] = serde_json::to_value(&data)?;
            writeln!(writer, "{}", vega_spec)?;
        }
        DotplotoutFormat::Html => {
            let mut env = Environment::new();
            env.add_template("vega", VEGA_TEMP)?;
            let template = env.get_template("vega")?;
            let n_records = data.len();
            vega_spec["data"]["values"] = serde_json::to_value(&data)?;
            // escape `</` so chromosome names can never close the <script> tag
            let vl_json = serde_json::to_string(&vega_spec)?.replace("</", "<\\/");
            let rendered = template.render(context! {
                vl_json => vl_json,
                mode => page.mode,
                record_kind => page.record_kind,
                extra_help => page.extra_help,
                n_records => n_records,
                n_ref => page.n_ref,
                n_query => page.n_query,
            })?;
            writeln!(writer, "{}", rendered)?;
        }
        DotplotoutFormat::Csv => {
            let mut wtr = csv::Writer::from_writer(writer);
            for record in data {
                wtr.serialize(record)?;
            }
            wtr.flush()?;
        }
    }
    Ok(())
}

type Collected<T> = (Vec<T>, ChromLens);

fn merge_collected<T>(mut a: Collected<T>, mut b: Collected<T>) -> Result<Collected<T>, WGAError> {
    a.0.append(&mut b.0);
    Ok((a.0, a.1.merge(b.1)))
}

/// Generate Plotdatas from MAF records
fn generate_maf_data<R: Read + Send>(
    mut reader: MAFReader<R>,
    no_identity: bool,
    query_name: Option<&str>,
) -> Result<Collected<AllPlotdata>, WGAError> {
    reader
        .records()
        .par_bridge()
        .try_fold(Collected::default, |mut acc, rec| {
            let mut rec = rec?;
            if let Some(qname) = query_name {
                rec.set_query_idx_byname(qname)?;
            }
            acc.1.add(&rec);
            acc.0.push(rec_dot_data(&rec, no_identity)?);
            Ok::<Collected<AllPlotdata>, WGAError>(acc)
        })
        .try_reduce(Collected::default, merge_collected)
}

/// Generate Plotdatas from PAF records
fn generate_paf_data<R: Read + Send>(
    mut reader: PAFReader<R>,
    no_identity: bool,
) -> Result<Collected<AllPlotdata>, WGAError> {
    reader
        .records()
        .par_bridge()
        .try_fold(Collected::default, |mut acc, rec| {
            let rec = rec?;
            acc.1.add(&rec);
            acc.0.push(rec_dot_data(&rec, no_identity)?);
            Ok::<Collected<AllPlotdata>, WGAError>(acc)
        })
        .try_reduce(Collected::default, merge_collected)
}

/// Generate BasePlotdatas from PAF records
fn generate_paf_basedata<R: Read + Send>(
    mut reader: PAFReader<R>,
    cutoff: usize,
) -> Result<Collected<Vec<BasePlotdata>>, WGAError> {
    reader
        .records()
        .par_bridge()
        .try_fold(Collected::default, |mut acc, rec| {
            let rec = rec?;
            acc.1.add(&rec);
            acc.0.push(parse_cigar_to_base_plotdata(&rec, cutoff)?);
            Ok::<Collected<Vec<BasePlotdata>>, WGAError>(acc)
        })
        .try_reduce(Collected::default, merge_collected)
}

/// Generate BasePlotdatas from MAF records
fn generate_maf_basedata<R: Read + Send>(
    mut reader: MAFReader<R>,
    cutoff: usize,
    query_name: Option<&str>,
) -> Result<Collected<Vec<BasePlotdata>>, WGAError> {
    reader
        .records()
        .par_bridge()
        .try_fold(Collected::default, |mut acc, rec| {
            let mut rec = rec?;
            if let Some(qname) = query_name {
                rec.set_query_idx_byname(qname)?;
            }
            acc.1.add(&rec);
            acc.0.push(parse_maf_to_base_plotdata(&rec, cutoff)?);
            Ok::<Collected<Vec<BasePlotdata>>, WGAError>(acc)
        })
        .try_reduce(Collected::default, merge_collected)
}

// stat a record to generate a Plotdata
fn rec_dot_data<T: AlignRecord>(rec: &T, no_identity: bool) -> Result<AllPlotdata, WGAError> {
    // get pair
    let ref_start = rec.target_start();
    let mut query_start = rec.query_start();
    let ref_end = rec.target_end();
    let mut query_end = rec.query_end();
    let identity = if no_identity {
        1.0
    } else {
        calculate_identity(rec)?
    };
    let ref_chro = rec.target_name().to_string();
    let query_chro = rec.query_name().to_string();
    let strand = rec.query_strand();
    match strand {
        Strand::Positive => {}
        Strand::Negative => {
            // reverse the query start and end
            std::mem::swap(&mut query_start, &mut query_end);
        }
    }
    Ok(AllPlotdata {
        ref_start,
        ref_end,
        query_start,
        query_end,
        identity,
        ref_chro,
        query_chro,
    })
}

// calculate identity for a record
fn calculate_identity<T: AlignRecord>(rec: &T) -> Result<f64, WGAError> {
    let aligned_size = rec.target_align_size();
    let rec_stat = rec.get_stat()?;
    let matched = rec_stat.matched;
    let identity = matched as f64 / aligned_size as f64;
    Ok(identity)
}
