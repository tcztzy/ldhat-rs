use clap::Parser;
use ldhat::commands::{Convert, Executable};
use ldhat::LDhatResult as Result;

#[derive(Parser)]
struct LDhatCLIOptions {
    #[command(subcommand)]
    action: LDhatAction,
    #[command(flatten)]
    verbosity: clap_verbosity_flag::Verbosity,
}

#[derive(Parser)]
enum LDhatAction {
    Convert(Convert),
}

impl Executable for LDhatAction {
    fn execute(&self) -> Result<()> {
        match self {
            Self::Convert(options) => options.execute(),
        }
    }
}

impl Executable for LDhatCLIOptions {
    fn execute(&self) -> Result<()> {
        env_logger::Builder::new()
            .filter_level(self.verbosity.log_level_filter())
            .init();
        self.action.execute()
    }
}

fn main() -> Result<()> {
    let args = LDhatCLIOptions::parse();
    args.execute()?;
    Ok(())
}
