use sage_core::database::Builder;
use sage_core::fasta::Fasta;
use sage_core::modification::VarModEntry;
use std::collections::HashMap;

const UBIQUITIN_SEQ: &str = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

fn build_fasta(decoy_tag: &str, generate_decoys: bool) -> Fasta {
    let fasta_str = format!(">sp|P0CG48|UBB_HUMAN Ubiquitin\n{}\n", UBIQUITIN_SEQ);
    Fasta::parse(fasta_str, decoy_tag, generate_decoys)
}

struct Config {
    label: &'static str,
    variable_mods: Option<HashMap<String, Vec<VarModEntry>>>,
    max_variable_mods: Option<usize>,
    max_combinations: Option<usize>,
}

fn run_config(cfg: &Config) -> (usize, usize) {
    let fasta = build_fasta("rev_", true);
    let db = Builder {
        fasta: Some("dummy".into()), // required field, unused since we pass fasta directly
        variable_mods: cfg.variable_mods.clone(),
        max_variable_mods: cfg.max_variable_mods,
        max_combinations: cfg.max_combinations,
        generate_decoys: Some(true),
        ..Default::default()
    }
    .make_parameters()
    .build(fasta);

    (db.peptides.len(), db.fragments.len())
}

fn main() {
    let configs = vec![
        Config {
            label: "1. No variable mods",
            variable_mods: None,
            max_variable_mods: None,
            max_combinations: None,
        },
        Config {
            label: "2. M oxidation [15.9949], max_variable_mods=2, no cap",
            variable_mods: Some(HashMap::from([(
                "M".into(),
                vec![VarModEntry::Mass(15.9949)],
            )])),
            max_variable_mods: Some(2),
            max_combinations: None,
        },
        Config {
            label: "3. M oxidation MassWithLimit(1 per peptide), max_variable_mods=2, no cap",
            variable_mods: Some(HashMap::from([(
                "M".into(),
                vec![VarModEntry::MassWithLimit(15.9949, 1)],
            )])),
            max_variable_mods: Some(2),
            max_combinations: None,
        },
        Config {
            label: "4. M oxidation [15.9949], max_variable_mods=2, max_combinations=5",
            variable_mods: Some(HashMap::from([(
                "M".into(),
                vec![VarModEntry::Mass(15.9949)],
            )])),
            max_variable_mods: Some(2),
            max_combinations: Some(5),
        },
    ];

    println!("{:<70} {:>10} {:>13} {:>18}", "Config", "Peptides", "Fragments", "Rough Mem (bytes)");
    println!("{}", "-".repeat(115));

    for cfg in &configs {
        let (peptides, fragments) = run_config(cfg);
        let mem = peptides * 200 + fragments * 8;
        println!(
            "{:<70} {:>10} {:>13} {:>18}",
            cfg.label, peptides, fragments, mem
        );
    }
}
