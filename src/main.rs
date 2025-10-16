use log::{error, info};
use wgalib::cli::{make_cli_parse, Commands};
use wgalib::errors::WGAError;
use wgalib::log::init_logger;
use wgalib::parser::common::FileFormat;
use wgalib::tools::tview::tview;
use wgalib::utils::{
    wrap_build_index, wrap_chain2maf, wrap_chain2paf, wrap_chunk, wrap_dotplot, wrap_filter,
    wrap_gencomp, wrap_maf2chain, wrap_maf2paf, wrap_maf2sam, wrap_maf_call, wrap_maf_extract,
    wrap_paf2chain, wrap_paf2maf, wrap_paf_call, wrap_paf_cov, wrap_paf_pesudo_maf,
    wrap_rename_maf, wrap_stat, wrap_validate,
};

fn main() {
    match main_entry() {
        Ok(_) => {}
        Err(e) => {
            error!("{}", e);
            std::process::exit(1);
        }
    }
}

fn main_entry() -> Result<(), WGAError> {
    let cli = make_cli_parse();
    let verbose = cli.verbose;

    init_logger(verbose);

    rayon::ThreadPoolBuilder::new()
        .num_threads(cli.threads)
        .build_global()?;

    let outfile = cli.outfile;
    let rewrite = cli.rewrite;

    // Info log
    info!("Command: {:?}", &cli.command);

    match &cli.command {
        Commands::Maf2Paf { input, query_name } => {
            wrap_maf2paf(input, &outfile, query_name.clone(), rewrite)?;
        }
        Commands::Paf2Maf {
            input,
            target,
            query,
        } => {
            wrap_paf2maf(input, &outfile, target, query, rewrite)?;
        }
        Commands::Paf2Chain { input } => {
            wrap_paf2chain(input, &outfile, rewrite)?;
        }
        Commands::Chain2Paf { input } => {
            wrap_chain2paf(input, &outfile, rewrite)?;
        }
        Commands::Chain2Maf {
            input,
            target,
            query,
        } => {
            wrap_chain2maf(input, &outfile, target, query, rewrite)?;
        }
        Commands::Maf2Chain { input, query_name } => {
            wrap_maf2chain(input, &outfile, rewrite, query_name.clone())?;
        }
        Commands::MafExtract {
            input,
            regions,
            file,
        } => {
            wrap_maf_extract(input, regions, file, &outfile, rewrite)?;
        }
        Commands::Call {
            input,
            sample,
            snp,
            inv,
            svlen,
            format,
            target,
            query,
            query_regex,
            chunk_size,
        } => match format {
            FileFormat::Maf => {
                wrap_maf_call(
                    input,
                    &outfile,
                    rewrite,
                    *snp,
                    *inv,
                    *svlen,
                    false,
                    sample.as_deref(),
                    query_regex,
                    *chunk_size,
                )?;
            }
            FileFormat::Paf => {
                let (target, query) = match (target, query) {
                    (Some(t), Some(q)) => (t, q),
                    _ => {
                        return Err(WGAError::Other(anyhow::anyhow!(
                            "target and query are necessary"
                        )));
                    }
                };
                wrap_paf_call(
                    input,
                    target,
                    query,
                    &outfile,
                    rewrite,
                    *snp,
                    *svlen,
                    true,
                    sample.as_deref(),
                )?;
            }
            _ => {
                return Err(WGAError::Other(anyhow::anyhow!("format is not supported")));
            }
        },
        Commands::Maf2Sam { input } => {
            wrap_maf2sam(input, &outfile, rewrite)?;
        }
        Commands::MafIndex { input } => {
            wrap_build_index(input, &outfile)?;
        }
        Commands::Tview { input, step } => {
            tview(input, *step)?;
        }
        Commands::Stat {
            input,
            format,
            each,
            query_name,
        } => wrap_stat(*format, input, &outfile, query_name.clone(), rewrite, *each)?,
        Commands::Dotplot {
            input,
            format,
            out_format,
            no_identity,
            length,
            mode,
            query_name,
            color,
        } => {
            wrap_dotplot(
                input,
                *format,
                *out_format,
                *mode,
                *no_identity,
                *length,
                &outfile,
                query_name.clone(),
                rewrite,
                color.clone(),
            )?;
        }
        Commands::Filter {
            input,
            format,
            min_block_size,
            min_query_size,
            min_align_size,
        } => {
            wrap_filter(
                *format,
                input,
                &outfile,
                rewrite,
                *min_block_size,
                *min_query_size,
                *min_align_size,
            )?;
        }
        Commands::Rename { input, prefixs } => {
            wrap_rename_maf(input, &outfile, rewrite, prefixs)?;
        }
        Commands::PafCov { input } => {
            wrap_paf_cov(input, &outfile, rewrite)?;
        }
        Commands::PafPseudo {
            input,
            fasta,
            target,
        } => {
            wrap_paf_pesudo_maf(input, &outfile, rewrite, fasta, target)?;
        } // Commands::TrimOvp { input } => {
        //     wrap_paf_trim_overlap(input, &outfile, rewrite)?;
        // }
        Commands::Chunk { input, length } => {
            wrap_chunk(input, &outfile, rewrite, *length)?;
        }
        Commands::GenCompletion { shell } => {
            wrap_gencomp(*shell, &outfile, rewrite)?;
        }
        Commands::Validate { input, fix } => {
            wrap_validate(input, fix, &outfile, rewrite)?;
        }
    }
    Ok(())
}
