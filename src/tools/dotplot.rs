use base64::{engine::general_purpose, Engine as _};
use std::io::prelude::*;

pub fn chart() {
    let schema = r#"
    {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "data": {
          "values": [
            {"x": 0, "y": 0, "color": 10},
            {"x": 10, "y": 15, "color": 10},
            {"x": 10, "y": null, "color": 10},
            {"x": 12, "y": 27, "color": 20},
            {"x": 18, "y": 18, "color": 20},
            {"x": 18, "y": null, "color": 20},
            {"x": 20, "y": 30, "color": 30},
            {"x": 40, "y": 60, "color": 30},
            {"x": 40, "y": null, "color": 30}
          ]
        },
        "mark": "line",
        "encoding": {
          "x": {
            "field": "x",
            "type": "quantitative",
            "axis": {
              "values": [20, 40],
              "labelExpr": "{20: 'chrA', 40: 'chrb'}[datum.value]"

            }
          },
          "y": {"field": "y", "type": "quantitative", "axis": {"values": [30, 60]}},
          "color": {"field": "color", "type": "quantitative"}
        }
      }
    "#;
    let code = general_purpose::STANDARD.encode(schema);
    let get_url = format!("https://kroki.io/vegalite/svg/{}", code);
    let res = reqwest::blocking::get(get_url).unwrap().text().unwrap();
    let stdout = std::io::stdout();
    let mut handle = stdout.lock();
    handle.write_all(res.as_bytes()).unwrap();
}
