//! Environment catalogue — 3-step matching: hash → fingerprint → geometry.

use std::collections::HashMap;
use super::mechanism::Mechanism;
use super::environment::{Geometry, Fingerprint};

/// A catalogued environment with its mechanisms.
#[derive(Debug, Clone)]
pub struct CatalogueEntry {
    pub geometry: Geometry,
    pub fingerprint: Fingerprint,
    pub cat_index: usize,
    pub delta_max: f64,
    pub mechanisms: Vec<Mechanism>,
    pub freq: usize,
    pub false_pos: usize,
}

/// Per-atom runtime environment data.
pub struct AtomEnv {
    pub geo: Geometry,
    pub fingerprint: Fingerprint,
    pub graph_hash: u64,
    pub entry_idx: Option<(u64, usize)>, // (hash, bucket_index)
}

/// Environment-to-mechanism catalogue.
pub struct Catalogue {
    pub r_env: f64,
    pub r_edge: f64,
    pub delta_max: f64,
    pub overfuzz: f64,
    pub min_delta_max: f64,

    buckets: HashMap<u64, Vec<CatalogueEntry>>,
    size: usize,
    atoms: Vec<AtomEnv>,
}

impl Catalogue {
    pub fn new(r_env: f64, r_edge: f64, delta_max: f64) -> Self {
        Catalogue {
            r_env,
            r_edge,
            delta_max,
            overfuzz: 1.5,
            min_delta_max: 0.01,
            buckets: HashMap::new(),
            size: 0,
            atoms: Vec::new(),
        }
    }

    pub fn n_environments(&self) -> usize {
        self.size
    }

    pub fn n_mechanisms_total(&self) -> usize {
        self.buckets.values()
            .flat_map(|b| b.iter())
            .map(|e| e.mechanisms.len())
            .sum()
    }

    /// Rebuild catalogue. Returns indices of new (unknown) environments.
    pub fn rebuild(
        &mut self,
        geometries: Vec<Geometry>,
    ) -> Vec<usize> {
        self.atoms.clear();
        let mut new_indices = Vec::new();
        let mut new_entries: Vec<(u64, usize)> = Vec::new(); // (hash, bucket_idx)

        for mut geo in geometries {
            let fp = Fingerprint::new(&geo.positions);
            let gh = geo.graph_hash(self.r_edge);

            // Ensure bucket exists
            self.buckets.entry(gh).or_insert_with(Vec::new);

            let entry_idx = self.find_match(gh, &fp, &geo);

            if let Some(idx) = entry_idx {
                self.buckets.get_mut(&gh).unwrap()[idx].freq += 1;
                self.atoms.push(AtomEnv {
                    geo,
                    fingerprint: fp,
                    graph_hash: gh,
                    entry_idx: Some((gh, idx)),
                });
            } else {
                // Check against other new environments
                let mut found_dup = false;
                for &(nh, ni) in &new_entries {
                    if gh == nh {
                        let entry = &self.buckets[&nh][ni];
                        let delta = self.calc_delta(&fp, entry);
                        if entry.fingerprint.equiv(&fp, delta * self.overfuzz) {
                            if geo.permute_onto(&entry.geometry, delta).is_some() {
                                self.atoms.push(AtomEnv {
                                    geo: geo.clone(),
                                    fingerprint: fp.clone(),
                                    graph_hash: gh,
                                    entry_idx: Some((nh, ni)),
                                });
                                found_dup = true;
                                break;
                            }
                        }
                    }
                }

                if !found_dup {
                    let atom_idx = self.atoms.len();
                    let bucket_idx = self.insert(gh, &geo, &fp);
                    new_indices.push(atom_idx);
                    new_entries.push((gh, bucket_idx));

                    self.atoms.push(AtomEnv {
                        geo,
                        fingerprint: fp,
                        graph_hash: gh,
                        entry_idx: Some((gh, bucket_idx)),
                    });
                }
            }
        }

        new_indices
    }

    /// Get mechanisms for an atom.
    pub fn get_mechanisms(&self, atom_index: usize) -> &[Mechanism] {
        if let Some((h, idx)) = &self.atoms[atom_index].entry_idx {
            &self.buckets[h][*idx].mechanisms
        } else {
            &[]
        }
    }

    /// Get catalogue entry for an atom.
    pub fn get_entry(&self, atom_index: usize) -> Option<&CatalogueEntry> {
        if let Some((h, idx)) = &self.atoms[atom_index].entry_idx {
            Some(&self.buckets[h][*idx])
        } else {
            None
        }
    }

    /// Set mechanisms for a newly discovered environment.
    pub fn set_mechanisms(&mut self, atom_index: usize, mechanisms: Vec<Mechanism>) {
        if let Some((h, idx)) = self.atoms[atom_index].entry_idx {
            self.buckets.get_mut(&h).unwrap()[idx].mechanisms = mechanisms;
        }
    }

    // ---- Internal ----

    fn calc_delta(&self, fp: &Fingerprint, entry: &CatalogueEntry) -> f64 {
        (fp.r_min * 0.4).min(entry.delta_max)
    }

    fn find_match(&self, hash: u64, fp: &Fingerprint, geo: &Geometry) -> Option<usize> {
        let bucket = self.buckets.get(&hash)?;
        for (i, entry) in bucket.iter().enumerate() {
            let delta = self.calc_delta(fp, entry);
            if !entry.fingerprint.equiv(fp, delta * self.overfuzz) {
                continue;
            }
            if geo.permute_onto(&entry.geometry, delta).is_some() {
                return Some(i);
            }
        }
        None
    }

    fn insert(&mut self, hash: u64, geo: &Geometry, fp: &Fingerprint) -> usize {
        let delta_max = (fp.r_min * 0.4).min(self.delta_max);
        let entry = CatalogueEntry {
            geometry: geo.clone(),
            fingerprint: fp.clone(),
            cat_index: self.size,
            delta_max,
            mechanisms: Vec::new(),
            freq: 1,
            false_pos: 0,
        };
        self.size += 1;
        let bucket = self.buckets.entry(hash).or_insert_with(Vec::new);
        let idx = bucket.len();
        bucket.push(entry);
        idx
    }
}
