use wgalib::cli::{make_cli_parse, Commands};
use wgalib::log::init_logger;
use wgalib::utils::{chain2maf, chain2paf, maf2chain, maf2paf, paf2chain, paf2maf};

fn main() {
    init_logger();

    let cli = make_cli_parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()
        .unwrap(); // TODO: handle threadpool error

    let outfile = cli.outfile;
    let rewrite = cli.rewrite;

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
    }
}
