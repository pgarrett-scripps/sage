use std::{
    collections::HashMap,
    fmt::{Display, Write},
    str::FromStr,
};

use serde::{Deserialize, Serialize};

/// A variable modification entry: either a bare mass or a (mass, max_count) pair.
/// When `max_count` is specified, at most that many instances of this modification
/// (at this exact mass) are allowed on a single peptide.
#[derive(Clone, Debug, Deserialize)]
#[serde(untagged)]
pub enum VarModEntry {
    Mass(f32),
    MassWithLimit(f32, usize),
}

impl VarModEntry {
    pub fn mass(&self) -> f32 {
        match self {
            VarModEntry::Mass(m) | VarModEntry::MassWithLimit(m, _) => *m,
        }
    }

    pub fn limit(&self) -> Option<usize> {
        match self {
            VarModEntry::Mass(_) => None,
            VarModEntry::MassWithLimit(_, l) => Some(*l),
        }
    }
}

use crate::mass::VALID_AA;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModificationSpecificity {
    PeptideN(Option<u8>),
    PeptideC(Option<u8>),
    ProteinN(Option<u8>),
    ProteinC(Option<u8>),
    Residue(u8),
}

impl Display for ModificationSpecificity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let r = match self {
            ModificationSpecificity::PeptideN(r) => {
                f.write_char('^')?;
                *r
            }
            ModificationSpecificity::PeptideC(r) => {
                f.write_char('$')?;
                *r
            }
            ModificationSpecificity::ProteinN(r) => {
                f.write_char('[')?;
                *r
            }
            ModificationSpecificity::ProteinC(r) => {
                f.write_char(']')?;
                *r
            }
            ModificationSpecificity::Residue(r) => Some(*r),
        };

        if let Some(r) = r {
            f.write_char(r as char)?;
        }

        Ok(())
    }
}

impl Serialize for ModificationSpecificity {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum InvalidModification {
    Empty,
    InvalidResidue(char),
    TooLong(String),
}

impl FromStr for ModificationSpecificity {
    type Err = InvalidModification;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.len() > 2 {
            return Err(InvalidModification::TooLong(s.into()));
        }
        if let Some(rest) = s.strip_prefix('^') {
            return Ok(ModificationSpecificity::PeptideN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('$') {
            return Ok(ModificationSpecificity::PeptideC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix('[') {
            return Ok(ModificationSpecificity::ProteinN(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        if let Some(rest) = s.strip_prefix(']') {
            return Ok(ModificationSpecificity::ProteinC(
                rest.chars().next().map(|ch| ch as u8),
            ));
        }
        match s.chars().next() {
            Some(c) => {
                if VALID_AA.contains(&(c as u8)) {
                    Ok(ModificationSpecificity::Residue(c as u8))
                } else {
                    Err(InvalidModification::InvalidResidue(c))
                }
            }
            None => Err(InvalidModification::Empty),
        }
    }
}

pub fn validate_mods(input: Option<HashMap<String, f32>>) -> HashMap<ModificationSpecificity, f32> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, mass) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(m, mass);
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!("Invalid modification string: unrecognized residue ({})", c)
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

pub fn validate_var_mods(
    input: Option<HashMap<String, Vec<VarModEntry>>>,
) -> HashMap<ModificationSpecificity, Vec<(f32, Option<usize>)>> {
    let mut output = HashMap::new();
    if let Some(input) = input {
        for (s, entries) in input {
            match ModificationSpecificity::from_str(&s) {
                Ok(m) => {
                    output.insert(
                        m,
                        entries.iter().map(|e| (e.mass(), e.limit())).collect(),
                    );
                }
                Err(InvalidModification::Empty) => {
                    log::error!("Skipping invalid modification string: empty")
                }
                Err(InvalidModification::InvalidResidue(c)) => {
                    log::error!(
                        "Skipping invalid modification string: unrecognized residue ({})",
                        c
                    )
                }
                Err(InvalidModification::TooLong(s)) => {
                    log::error!("Skipping invalid modification string: {} is too long", s)
                }
            }
        }
    }
    output
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn parse_modifications() {
        use InvalidModification::*;
        use ModificationSpecificity::*;
        assert_eq!("[".parse::<ModificationSpecificity>(), Ok(ProteinN(None)));
        assert_eq!(
            "[M".parse::<ModificationSpecificity>(),
            Ok(ProteinN(Some(b'M')))
        );
        assert_eq!(
            "]M".parse::<ModificationSpecificity>(),
            Ok(ProteinC(Some(b'M')))
        );
        assert_eq!("M".parse::<ModificationSpecificity>(), Ok(Residue(b'M')));
        assert_eq!(
            "Z".parse::<ModificationSpecificity>(),
            Err(InvalidResidue('Z'))
        );
    }

    #[test]
    fn var_mod_entry_bare_mass() {
        let entry = VarModEntry::Mass(15.9949);
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.limit(), None);
    }

    #[test]
    fn var_mod_entry_mass_with_limit() {
        let entry = VarModEntry::MassWithLimit(15.9949, 1);
        assert_eq!(entry.mass(), 15.9949);
        assert_eq!(entry.limit(), Some(1));
    }

    #[test]
    fn validate_var_mods_mixed() {
        use ModificationSpecificity::*;
        // Mix bare masses and MassWithLimit entries
        let mut raw = HashMap::new();
        raw.insert(
            "M".to_string(),
            vec![VarModEntry::Mass(15.9949), VarModEntry::MassWithLimit(15.9949, 1)],
        );
        raw.insert(
            "C".to_string(),
            vec![VarModEntry::MassWithLimit(57.0215, 2)],
        );
        let result = validate_var_mods(Some(raw));

        let m_entries = result.get(&Residue(b'M')).unwrap();
        assert_eq!(m_entries.len(), 2);
        assert!((m_entries[0].0 - 15.9949).abs() < 1e-4);
        assert_eq!(m_entries[0].1, None);
        assert!((m_entries[1].0 - 15.9949).abs() < 1e-4);
        assert_eq!(m_entries[1].1, Some(1));

        let c_entries = result.get(&Residue(b'C')).unwrap();
        assert_eq!(c_entries.len(), 1);
        assert!((c_entries[0].0 - 57.0215).abs() < 1e-4);
        assert_eq!(c_entries[0].1, Some(2));
    }

    #[test]
    fn validate_var_mods_invalid_residue_skipped() {
        let mut raw = HashMap::new();
        raw.insert("Z".to_string(), vec![VarModEntry::Mass(15.9949)]);
        raw.insert("M".to_string(), vec![VarModEntry::Mass(15.9949)]);
        let result = validate_var_mods(Some(raw));
        // Z is invalid — only M should survive
        assert_eq!(result.len(), 1);
        assert!(result.contains_key(&ModificationSpecificity::Residue(b'M')));
    }
}
