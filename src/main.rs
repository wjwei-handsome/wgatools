use log::error;
use wgalib::cli::{make_cli_parse, Commands};
use wgalib::log::init_logger;
use wgalib::tools::tview::tview;
use wgalib::utils::{
    chain2maf, chain2paf, maf2chain, maf2paf, maf2sam, paf2chain, paf2maf, wrap_build_index,
    wrap_filter, wrap_maf_call, wrap_maf_extract, wrap_stat,
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
        Commands::Call {
            input,
            sample,
            snp,
            svlen,
        } => {
            wrap_maf_call(
                input,
                &outfile,
                rewrite,
                *snp,
                *svlen,
                false,
                sample.as_deref(),
            );
        }
        Commands::Maf2Sam { input } => maf2sam(input, &outfile, rewrite),
        Commands::MafIndex { input } => wrap_build_index(input, &outfile),
        Commands::Tview { input, step } => match tview(input, *step) {
            Ok(_) => {}
            Err(err) => {
                error!("{}", err);
                std::process::exit(1);
            }
        },
        Commands::Stat {
            input,
            format,
            each,
        } => {
            wrap_stat(*format, input, &outfile, rewrite, *each);
        }
        Commands::Dotplot {} => {
            wgalib::tools::dotplot::chart();
        }
        Commands::Filter {
            input,
            format,
            min_block_size,
            min_query_size,
        } => {
            wrap_filter(
                *format,
                input,
                &outfile,
                rewrite,
                *min_block_size,
                *min_query_size,
            );
        }
    }
}
