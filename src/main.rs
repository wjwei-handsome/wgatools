use log::error;
use wgalib::cli::{make_cli_parse, Commands};
use wgalib::log::init_logger;
use wgalib::tools::caller::maf_call;
use wgalib::tools::tview::tview;
use wgalib::utils::{
    chain2maf, chain2paf, maf2chain, maf2paf, maf2sam, paf2chain, paf2maf, wrap_build_index,
    wrap_maf_extract,
};

fn main() {
    let cli = make_cli_parse();
    let verbose = cli.verbose;

    init_logger(verbose);

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap(); // TODO: handle threadpool error

    let outfile = cli.outfile;
    let rewrite = cli.rewrite;

    // Info log
    // info!("Process: {:?}", &cli.command);

    match &cli.command {
        Commands::Maf2Paf { input } => {
            maf2paf(input, &outfile, rewrite);
        }
        Commands::Paf2Maf {
            input,
            target,
            query,
        } => {
            paf2maf(input, &outfile, target, query, rewrite);
        }
        Commands::Paf2Chain { input } => paf2chain(input, &outfile, rewrite),
        Commands::Chain2Paf { input } => chain2paf(input, &outfile, rewrite),
        Commands::Chain2Maf {
            input,
            target,
            query,
        } => chain2maf(input, &outfile, target, query, rewrite),
        Commands::Maf2Chain { input } => maf2chain(input, &outfile, rewrite),
        Commands::MafExtract {
            input,
            regions,
            file,
        } => {
            wrap_maf_extract(input, regions, file, &outfile, rewrite);
        }
        Commands::Call { input: _ } => {
            maf_call();
        }
        Commands::Maf2Sam { input } => maf2sam(input, &outfile, rewrite),
        Commands::MafIndex { input } => wrap_build_index(input, &outfile),
        Commands::Tview { input, step } => match tview(input, *step) {
            Ok(_) => {}
            Err(err) => {
                error!("Error: {}", err);
                std::process::exit(1);
            }
        },
    }
}
