//! Reproduction harness for the reported `prefilter_low_memory` bug.
//!
//! Claim under test (from the bug report): with `prefilter_low_memory: true`,
//! the prefilter keeps only the top-`report_psms` candidates per spectrum ranked
//! by a *cheap quick-score* that disagrees with the final hyperscore, so the
//! global hyperscore-best PSM can be pruned before full scoring.
//!
//! What the code actually does: the low-memory branch of `Scorer::quick_score`
//! calls `score_candidate` (the *full* hyperscore) and keeps the top-k by that
//! score. These tests check whether the highest-hyperscore candidate can ever be
//! dropped by the low-memory prefilter, and whether the true/false branches
//! diverge on the final identification through the real chunked pipeline.

use std::collections::HashMap;
use std::sync::atomic::{AtomicBool, Ordering};

use sage_core::database::{Builder, EnzymeBuilder, IndexedDatabase, Parameters, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::ion_series::{IonSeries, Kind};
use sage_core::mass::{Tolerance, PROTON};
use sage_core::modification::VariableMod;
use sage_core::scoring::{Scorer, ScoreType};
use sage_core::spectrum::{Peak, Precursor, ProcessedSpectrum};

const GLYGLY: f32 = 114.04293;

/// Variable GlyGly on K, like the diGly search in the report.
fn variable_mods() -> HashMap<String, Vec<VariableMod>> {
    let mut m = HashMap::new();
    m.insert(
        "K".to_string(),
        vec![VariableMod::Simple(
            serde_json::from_str("114.04293").unwrap(),
        )],
    );
    m
}

fn build_params(chunk_size: usize) -> Parameters {
    Builder {
        bucket_size: Some(8192),
        // Non-specific-ish: keep each protein as a single peptide so we control
        // the candidate set exactly. Trypsin would split on internal K/R; our
        // sequences only have a single terminal K.
        enzyme: Some(EnzymeBuilder {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            ..Default::default()
        }),
        variable_mods: Some(variable_mods()),
        max_variable_mods: Some(1),
        generate_decoys: Some(false),
        decoy_tag: Some("rev_".into()),
        fasta: Some("static".into()),
        prefilter_chunk_size: Some(chunk_size),
        ..Default::default()
    }
    .make_parameters()
}

/// Build a synthetic MS2 from the b/y ions (charge 1) of a fully-resolved
/// peptide, stored the way sage stores experimental peaks (neutral mass).
fn synthetic_spectrum(db: &IndexedDatabase, ix: PeptideIx) -> ProcessedSpectrum<Peak> {
    let peptide = &db[ix];
    let mut peaks: Vec<Peak> = Vec::new();
    for kind in [Kind::B, Kind::Y] {
        for ion in IonSeries::new(peptide, kind) {
            peaks.push(Peak {
                mass: ion.monoisotopic_mass,
                intensity: 1000.0,
            });
        }
    }
    peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
    let tic = peaks.iter().map(|p| p.intensity).sum();
    ProcessedSpectrum {
        level: 2,
        peaks,
        total_ion_current: tic,
        // charge-1 precursor at the peptide's neutral mass
        precursors: vec![Precursor {
            mz: peptide.monoisotopic + PROTON,
            charge: Some(1),
            ..Default::default()
        }],
        ..Default::default()
    }
}

fn scorer<'a>(db: &'a IndexedDatabase, report_psms: usize) -> Scorer<'a> {
    Scorer {
        db,
        // Wide precursor tolerance => every peptide in the db is a precursor
        // candidate, maximising window density (worst case for trim_hits).
        precursor_tol: Tolerance::Da(-5000.0, 5000.0),
        fragment_tol: Tolerance::Ppm(-10.0, 10.0),
        min_matched_peaks: 4,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 1,
        max_precursor_charge: 1,
        override_precursor_charge: false,
        max_fragment_charge: Some(1),
        chimera: false,
        report_psms,
        wide_window: false,
        annotate_matches: false,
        localize: false,
        mass_shift_ppm: 5.0,
        score_type: ScoreType::SageHyperScore,
    }
}

/// Find the db index of the GlyGly-modified form of `seq` (single mod on the K).
fn modified_ix(db: &IndexedDatabase, seq: &str) -> PeptideIx {
    db.peptides
        .iter()
        .position(|p| {
            p.sequence.as_ref() == seq.as_bytes()
                && p.modifications.iter().any(|m| (m - GLYGLY).abs() < 1e-3)
        })
        .map(|i| PeptideIx(i as u32))
        .expect("modified peptide form should exist in db")
}

/// Run one prefilter pass over `db` for `spectrum`, return the set of kept idxs.
fn prefilter_keep(
    db: &IndexedDatabase,
    spectrum: &ProcessedSpectrum<Peak>,
    report_psms: usize,
    low_memory: bool,
) -> Vec<bool> {
    let s = scorer(db, report_psms);
    let keep = (0..db.peptides.len())
        .map(|_| AtomicBool::new(false))
        .collect::<Vec<_>>();
    s.quick_score(spectrum, low_memory, &keep);
    keep.iter().map(|b| b.load(Ordering::Relaxed)).collect()
}

/// A FASTA where the GlyGly form of CORRECT is the true best match, surrounded by
/// many "filler" peptides that share fragment masses with it (so they earn high
/// *coarse* matched-peak counts and crowd the trim_hits / heap stages).
fn crowded_fasta() -> String {
    // Correct peptide: single terminal K, carries the diagnostic GlyGly site.
    let mut f = String::from(">CORRECT\nAAAGAAGAAGAADAAGAAGAAGK\n");
    // Fillers: near-anagrams reusing the same residue alphabet so many of their
    // b/y ions collide with CORRECT's, inflating their coarse matched counts.
    let fillers = [
        "GAAAGAAGAAGAADAAGAAGAAK",
        "AGAAGAAGAAAGAADAAGAAGAK",
        "AAGAAGAAGAADAGAAGAAAGAK",
        "AADAAGAAGAAGAAGAAAGAAGK",
        "GAAGAADAAGAAGAAAGAAGAAK",
        "AAGAADAAGAAGAAGAAAGAAGK",
        "AGAADAAGAAGAAAGAAGAAGAK",
        "AAGAAAGAAGAADAAGAAGAAGK",
    ];
    for (i, s) in fillers.iter().enumerate() {
        // Replicate each filler under distinct accessions to grow the window.
        for j in 0..8 {
            f.push_str(&format!(">FILL_{i}_{j}\n{s}\n"));
        }
    }
    f
}

#[test]
fn low_memory_keeps_the_hyperscore_best_candidate() {
    // Single chunk (chunk_size huge) => one db, the same one the final pass uses.
    let params = build_params(100_000);
    let fasta = Fasta::parse(crowded_fasta(), "rev_", false);
    let db = params.build(fasta);

    let correct = modified_ix(&db, "AAAGAAGAAGAADAAGAAGAAGK");
    let spectrum = synthetic_spectrum(&db, correct);

    // Ground truth: full scoring ranks `correct` #1 by hyperscore.
    let s = scorer(&db, 50);
    let features = s.score(&spectrum);
    assert_eq!(
        features[0].peptide_idx, correct,
        "sanity: full scoring should rank the GlyGly form first"
    );

    // The low-memory prefilter must keep the global hyperscore-best PSM, even at
    // report_psms=1. Regression guard for the derived-`PartialOrd` bug, where the
    // bounded min-heap selected by `PeptideIx` instead of `hyperscore`.
    let keep_true = prefilter_keep(&db, &spectrum, 1, true);
    let keep_false = prefilter_keep(&db, &spectrum, 1, false);

    assert!(
        keep_false[correct.0 as usize],
        "low_memory=false must keep the correct peptide"
    );
    assert!(
        keep_true[correct.0 as usize],
        "low_memory=true must keep the global hyperscore-best PSM (report_psms=1)"
    );
}

/// The low-memory prefilter must retain the *same* per-spectrum top-k set as the
/// final hyperscore ranking, for every reasonable `report_psms`.
#[test]
fn low_memory_keep_set_matches_hyperscore_ranking() {
    let params = build_params(100_000);
    let fasta = Fasta::parse(crowded_fasta(), "rev_", false);
    let db = params.build(fasta);

    let correct = modified_ix(&db, "AAAGAAGAAGAADAAGAAGAAGK");
    let spectrum = synthetic_spectrum(&db, correct);

    for report_psms in [1usize, 2, 3, 5] {
        // Expected: the top-`report_psms` PeptideIx by full hyperscore.
        let s = scorer(&db, report_psms);
        let mut expected: Vec<u32> = s
            .score(&spectrum)
            .iter()
            .take(report_psms)
            .map(|f| f.peptide_idx.0)
            .collect();
        expected.sort_unstable();

        let keep = prefilter_keep(&db, &spectrum, report_psms, true);
        let mut got: Vec<u32> = keep
            .iter()
            .enumerate()
            .filter(|(_, b)| **b)
            .map(|(i, _)| i as u32)
            .collect();
        got.sort_unstable();

        assert_eq!(
            got, expected,
            "low_memory keep-set != top-{report_psms} by hyperscore"
        );
    }
}
