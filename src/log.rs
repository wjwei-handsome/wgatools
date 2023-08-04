use log::LevelFilter;
use log4rs::{
    append::console::{ConsoleAppender, Target},
    config::{Appender, Config, Root},
    encode::pattern::PatternEncoder,
    filter::threshold::ThresholdFilter,
};

pub fn init_logger(verbose: u8) {
    let log_level = match verbose {
        0 => LevelFilter::Off,
        1 => LevelFilter::Trace,
        2 => LevelFilter::Debug,
        _ => LevelFilter::Info,
    };
    // Build a stderr logger.
    let log_stderr = ConsoleAppender::builder()
        .target(Target::Stderr)
        .encoder(Box::new(PatternEncoder::new("{d} {h({l})} {m}{n}")))
        .build();
    let log_config = Config::builder()
        .appender(
            Appender::builder()
                .filter(Box::new(ThresholdFilter::new(log_level)))
                .build("stderr", Box::new(log_stderr)),
        )
        .build(Root::builder().appender("stderr").build(log_level))
        .unwrap();
    // init logger using config
    log4rs::init_config(log_config).unwrap();
}
