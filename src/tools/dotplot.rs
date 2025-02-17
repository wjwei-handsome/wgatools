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

const DOTPLOT_SPEC: &str = r#"
{
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "height": 800,
    "width": 800,
    "data": {
        "values": []
    },
    "params": [
        {
        "name": "zoom",
        "select": "interval",
        "bind": "scales"
        },
        {
        "name": "cigartype",
        "select": {"type": "point", "fields": ["cigar"]},
        "bind": "legend"
        }
    ],
    "mark": {
        "type": "rule",
        "tooltip": true,
        "strokeCap": "round"
    },
    "transform": [{
        "calculate": "datum.ref_chro+':'+toString(datum.ref_start)+'-'+toString(datum.ref_end)",
        "as": "ref"
    }, {
        "calculate": "datum.query_chro+':'+toString(datum.query_start)+'-'+toString(datum.query_end)",
        "as": "query"
    }, {
        "calculate": "abs(datum.ref_end-datum.ref_start)",
        "as": "ref_len"
    },{
        "calculate": "abs(datum.query_end-datum.query_start)",
        "as": "query_len"
    }, {
        "as": "cigar",
        "calculate": "datum.cigar == 'M' && datum.query_start > datum.query_end ? 'M_R' : datum.cigar"
    }],
    "encoding": {
        "x": {
            "field": "ref_start",
            "type": "quantitative",
            "title":null
        },
        "y": {
            "field": "query_start",
            "type": "quantitative",
            "title":null
        },
        "x2": {
            "field": "ref_end"
        },
        "y2": {
            "field": "query_end"
        },
        "color": {
            "field": "identity",
            "type": "quantitative",
            "scale": {
                "scheme": "blues"
              },
            "legend": {
                "labelFontSize": 20,
                "symbolSize": 10,
                "symbolStrokeWidth": 10,
                "symbolType": "square"
              }
        },
        "tooltip": [{
            "field": "ref",
            "type": "nominal"
        }, {
            "field": "query",
            "type": "nominal"
        }, {
            "field": "identity",
            "type": "nominal"
        }, {
            "field": "ref_len",
            "type": "quantitative"
        },{
            "field": "query_len",
            "type": "quantitative"
        }],
        "column": {
            "field": "ref_chro",
            "title": null
        },
        "row": {
            "field": "query_chro",
            "header": {
                "labelAngle": 0
            },
            "sort": "descending",
            "title": null
        },
        "opacity": {
            "condition": {"param": "cigartype", "value": 1},
            "value": 0.2
          },
        "strokeWidth": {
            "condition": {"param": "cigartype", "value": 5},
            "value": 2
        }
    },
    "resolve": {"scale": {"x": "independent", "y": "independent"}}
}"#;

const VEGA_TEMP: &str = r#"<head>
    <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
</head>

<body>
    <div id="view" style="display: flex; justify-content: space-evenly;"></div>
    <script>
        const spec = {{ vl_json | safe }};
        vegaEmbed(
            '#view',
            spec
        );
    </script>
</body>
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
    // init vega spec
    let mut vega_spec: Value = serde_json::from_str(DOTPLOT_SPEC)?;

    // match mode to generate data
    match mode {
        DotplotMode::Overview => {
            let pair_stat_vec = match format {
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
            render_output(pair_stat_vec, writer, out_format, vega_spec)?;
        }
        DotplotMode::BaseLevel => {
            let pair_base_plot_vec = match format {
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

            // change the vega spec
            vega_spec["encoding"]["x"]["scale"]["zero"] = false.into();
            vega_spec["encoding"]["y"]["scale"]["zero"] = false.into();
            vega_spec["encoding"]["color"]["scale"]["scheme"] = "category10".into();
            vega_spec["encoding"]["color"]["field"] = "cigar".into();
            vega_spec["encoding"]["color"]["type"] = "nominal".into();
            vega_spec["encoding"]["tooltip"][2]["field"] = "cigar".into();

            if let Some(color_str) = color_config {
                let color_map = parse_color_config(color_str)?;
                let domain: Vec<String> = color_map.keys().cloned().collect();
                let range: Vec<String> = color_map.values().cloned().collect();

                vega_spec["encoding"]["color"]["scale"]["domain"] = json!(domain);
                vega_spec["encoding"]["color"]["scale"]["range"] = json!(range);
            }

            render_output(final_base_plotdata, writer, out_format, vega_spec)?;
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
            vega_spec["data"]["values"] = serde_json::to_value(&data)?;
            let vl_json = serde_json::to_string(&vega_spec)?;
            let rendered = template.render(context! { vl_json => vl_json })?;
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

/// Generate Plotdatas from MAF records
fn generate_maf_data<R: Read + Send>(
    mut reader: MAFReader<R>,
    no_identity: bool,
    query_name: Option<&str>,
) -> Result<Vec<AllPlotdata>, WGAError> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .try_fold(Vec::new, |mut acc, rec| {
            let mut rec = rec?;
            if let Some(qname) = query_name {
                rec.set_query_idx_byname(qname)?;
            }
            acc.push(rec_dot_data(&rec, no_identity)?);
            Ok::<Vec<AllPlotdata>, WGAError>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut vec| {
            acc.append(&mut vec);
            Ok(acc)
        })?;
    Ok(pair_stat_vec)
}

/// Generate Plotdatas from PAF records
fn generate_paf_data<R: Read + Send>(
    mut reader: PAFReader<R>,
    no_identity: bool,
) -> Result<Vec<AllPlotdata>, WGAError> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .try_fold(Vec::new, |mut acc, rec| {
            acc.push(rec_dot_data(&rec?, no_identity)?);
            Ok::<Vec<AllPlotdata>, WGAError>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut vec| {
            acc.append(&mut vec);
            Ok(acc)
        })?;
    Ok(pair_stat_vec)
}

/// Generate BasePlotdatas from PAF records
fn generate_paf_basedata<R: Read + Send>(
    mut reader: PAFReader<R>,
    cutoff: usize,
) -> Result<Vec<Vec<BasePlotdata>>, WGAError> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .try_fold(Vec::new, |mut acc, rec| {
            acc.push(parse_cigar_to_base_plotdata(&rec?, cutoff)?);
            Ok::<Vec<Vec<BasePlotdata>>, WGAError>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut vec| {
            // join the nested vec
            acc.append(&mut vec);
            Ok(acc)
        })?;
    Ok(pair_stat_vec)
}

/// Generate BasePlotdatas from MAF records
fn generate_maf_basedata<R: Read + Send>(
    mut reader: MAFReader<R>,
    cutoff: usize,
    query_name: Option<&str>,
) -> Result<Vec<Vec<BasePlotdata>>, WGAError> {
    let pair_stat_vec = reader
        .records()
        .par_bridge()
        .try_fold(Vec::new, |mut acc, rec| {
            let mut rec = rec?;
            if let Some(qname) = query_name {
                rec.set_query_idx_byname(qname)?;
            }
            acc.push(parse_maf_to_base_plotdata(&rec, cutoff)?);
            Ok::<Vec<Vec<BasePlotdata>>, WGAError>(acc)
        })
        .try_reduce(Vec::new, |mut acc, mut vec| {
            // join the nested vec
            acc.append(&mut vec);
            Ok(acc)
        })?;
    Ok(pair_stat_vec)
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
