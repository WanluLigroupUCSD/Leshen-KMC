/// Data model definitions for KMC model construction.
///
/// Mirrors the Python mykmc API: Project, Species, Process, Condition, Action,
/// LateralInteraction, BEPRelation.

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A chemical species that can occupy a lattice site.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Species {
    pub name: String,
    pub id: usize,
    pub color: Option<String>,
}

/// A parameter that can be adjusted at runtime.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Parameter {
    pub name: String,
    pub value: f64,
    pub adjustable: bool,
    pub min: f64,
    pub max: f64,
}

/// Offset coordinate for conditions/actions (relative to process center site).
pub type Offset = Vec<i32>; // length = ndim

/// A condition that must be met for a process to fire.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Condition {
    pub offset: Offset,
    pub species: String,
}

/// An action that changes a site when a process fires.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Action {
    pub offset: Offset,
    pub species: String,
}

/// Rate expression: either a numeric constant or a symbolic string.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum RateExpr {
    Constant(f64),
    Expression(String),
}

/// An elementary process (reaction step).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Process {
    pub name: String,
    pub id: usize,
    pub conditions: Vec<Condition>,
    pub actions: Vec<Action>,
    pub rate_constant: RateExpr,
    /// TOF counters: observable_name -> coefficient
    pub tof_count: HashMap<String, f64>,
    pub enabled: bool,
    /// Required site type (None = any site)
    #[serde(default)]
    pub site_type: Option<i32>,
    /// Name of the reverse process (for BEP pairing)
    #[serde(default)]
    pub reverse_of: Option<String>,
}

/// Pairwise lateral interaction between neighboring adsorbates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LateralInteraction {
    pub species1: String,
    pub species2: String,
    /// Interaction energy in eV. Positive = repulsive.
    pub energy: f64,
}

/// Bronsted-Evans-Polanyi relation for a process.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BEPRelation {
    pub process_name: String,
    /// BEP slope (proximity factor), typically 0.0-1.0.
    pub alpha: f64,
}

/// Complete model definition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Project {
    pub model_name: String,
    pub model_dimension: usize,
    pub species_list: Vec<Species>,
    pub parameter_list: Vec<Parameter>,
    pub process_list: Vec<Process>,
    pub cell_size: Vec<f64>,
    #[serde(default)]
    pub lateral_interactions: Vec<LateralInteraction>,
    #[serde(default)]
    pub bep_relations: Vec<BEPRelation>,

    // Lookup maps (not serialized)
    #[serde(skip)]
    pub species_map: HashMap<String, usize>,
    #[serde(skip)]
    pub parameter_map: HashMap<String, usize>,
}

impl Project {
    pub fn new(name: &str, ndim: usize) -> Self {
        Self {
            model_name: name.to_string(),
            model_dimension: ndim,
            species_list: Vec::new(),
            parameter_list: Vec::new(),
            process_list: Vec::new(),
            cell_size: vec![1.0; ndim],
            lateral_interactions: Vec::new(),
            bep_relations: Vec::new(),
            species_map: HashMap::new(),
            parameter_map: HashMap::new(),
        }
    }

    /// Add a species. First added becomes default (id=0).
    pub fn add_species(&mut self, name: &str) -> usize {
        let id = self.species_list.len();
        self.species_list.push(Species {
            name: name.to_string(),
            id,
            color: None,
        });
        self.species_map.insert(name.to_string(), id);
        id
    }

    pub fn add_parameter(&mut self, name: &str, value: f64) -> usize {
        let idx = self.parameter_list.len();
        self.parameter_list.push(Parameter {
            name: name.to_string(),
            value,
            adjustable: false,
            min: 0.0,
            max: 0.0,
        });
        self.parameter_map.insert(name.to_string(), idx);
        idx
    }

    pub fn add_parameter_adjustable(
        &mut self, name: &str, value: f64, min: f64, max: f64,
    ) -> usize {
        let idx = self.parameter_list.len();
        self.parameter_list.push(Parameter {
            name: name.to_string(),
            value,
            adjustable: true,
            min,
            max,
        });
        self.parameter_map.insert(name.to_string(), idx);
        idx
    }

    pub fn add_process(
        &mut self,
        name: &str,
        conditions: Vec<Condition>,
        actions: Vec<Action>,
        rate_constant: RateExpr,
        tof_count: HashMap<String, f64>,
    ) -> usize {
        let id = self.process_list.len();
        self.process_list.push(Process {
            name: name.to_string(),
            id,
            conditions,
            actions,
            rate_constant,
            tof_count,
            enabled: true,
            site_type: None,
            reverse_of: None,
        });
        id
    }

    /// Add a process with site type and reverse_of.
    pub fn add_process_ext(
        &mut self,
        name: &str,
        conditions: Vec<Condition>,
        actions: Vec<Action>,
        rate_constant: RateExpr,
        tof_count: HashMap<String, f64>,
        site_type: Option<i32>,
        reverse_of: Option<&str>,
    ) -> usize {
        let id = self.process_list.len();
        self.process_list.push(Process {
            name: name.to_string(),
            id,
            conditions,
            actions,
            rate_constant,
            tof_count,
            enabled: true,
            site_type,
            reverse_of: reverse_of.map(|s| s.to_string()),
        });
        id
    }

    /// Add a pairwise lateral interaction.
    pub fn add_lateral_interaction(&mut self, sp1: &str, sp2: &str, energy: f64) {
        self.lateral_interactions.push(LateralInteraction {
            species1: sp1.to_string(),
            species2: sp2.to_string(),
            energy,
        });
    }

    /// Add a BEP relation for a process.
    pub fn add_bep_relation(&mut self, process_name: &str, alpha: f64) {
        self.bep_relations.push(BEPRelation {
            process_name: process_name.to_string(),
            alpha,
        });
    }

    /// Add diffusion processes for a species in all NN directions.
    pub fn add_diffusion(
        &mut self,
        species: &str,
        empty: &str,
        rate_constant: RateExpr,
    ) {
        let ndim = self.model_dimension;
        let offsets: Vec<(Vec<i32>, &str)> = match ndim {
            1 => vec![(vec![1], "right"), (vec![-1], "left")],
            2 => vec![
                (vec![1, 0], "right"), (vec![-1, 0], "left"),
                (vec![0, 1], "up"), (vec![0, -1], "down"),
            ],
            3 => vec![
                (vec![1, 0, 0], "right"), (vec![-1, 0, 0], "left"),
                (vec![0, 1, 0], "up"), (vec![0, -1, 0], "down"),
                (vec![0, 0, 1], "front"), (vec![0, 0, -1], "back"),
            ],
            _ => vec![],
        };
        let zero: Vec<i32> = vec![0; ndim];
        for (off, label) in offsets {
            let name = format!("{}_diff_{}", species, label);
            self.add_process(
                &name,
                vec![
                    Condition::new(&zero, species),
                    Condition::new(&off, empty),
                ],
                vec![
                    Action::new(&zero, empty),
                    Action::new(&off, species),
                ],
                rate_constant.clone(),
                HashMap::new(),
            );
        }
    }

    pub fn species_id(&self, name: &str) -> usize {
        self.species_map[name]
    }

    pub fn get_parameter(&self, name: &str) -> f64 {
        let idx = self.parameter_map[name];
        self.parameter_list[idx].value
    }

    pub fn set_parameter(&mut self, name: &str, value: f64) {
        if let Some(&idx) = self.parameter_map.get(name) {
            self.parameter_list[idx].value = value;
        }
    }

    /// Rebuild lookup maps after deserialization.
    pub fn rebuild_maps(&mut self) {
        self.species_map.clear();
        for sp in &self.species_list {
            self.species_map.insert(sp.name.clone(), sp.id);
        }
        self.parameter_map.clear();
        for (i, p) in self.parameter_list.iter().enumerate() {
            self.parameter_map.insert(p.name.clone(), i);
        }
    }

    pub fn summary(&self) {
        println!("Model: {}", self.model_name);
        println!("  Dimension: {}", self.model_dimension);
        println!(
            "  Species ({}): {}",
            self.species_list.len(),
            self.species_list.iter().map(|s| s.name.as_str()).collect::<Vec<_>>().join(", ")
        );
        println!(
            "  Parameters ({}): {}",
            self.parameter_list.len(),
            self.parameter_list.iter().map(|p| format!("{}={}", p.name, p.value)).collect::<Vec<_>>().join(", ")
        );
        println!("  Processes ({}):", self.process_list.len());
        for p in &self.process_list {
            let conds: Vec<String> = p.conditions.iter()
                .map(|c| format!("{}@{:?}", c.species, c.offset))
                .collect();
            let acts: Vec<String> = p.actions.iter()
                .map(|a| format!("{}@{:?}", a.species, a.offset))
                .collect();
            println!("    {}: [{}] -> [{}]", p.name, conds.join(", "), acts.join(", "));
        }
        if !self.lateral_interactions.is_empty() {
            println!("  Lateral interactions ({}):", self.lateral_interactions.len());
            for li in &self.lateral_interactions {
                let kind = if li.energy > 0.0 { "repulsive" } else { "attractive" };
                println!("    {}-{}: {:+.4} eV ({})", li.species1, li.species2, li.energy, kind);
            }
        }
        if !self.bep_relations.is_empty() {
            println!("  BEP relations ({}):", self.bep_relations.len());
            for b in &self.bep_relations {
                println!("    {}: alpha={}", b.process_name, b.alpha);
            }
        }
    }

    /// Save model to JSON file.
    pub fn save(&self, path: &str) -> std::io::Result<()> {
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        std::fs::write(path, json)?;
        println!("Model saved to {}", path);
        Ok(())
    }

    /// Load model from JSON file.
    pub fn load(path: &str) -> std::io::Result<Self> {
        let json = std::fs::read_to_string(path)?;
        let mut project: Self = serde_json::from_str(&json)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        project.rebuild_maps();
        Ok(project)
    }
}

// Convenience constructors
impl Condition {
    pub fn new(offset: &[i32], species: &str) -> Self {
        Self {
            offset: offset.to_vec(),
            species: species.to_string(),
        }
    }
}

impl Action {
    pub fn new(offset: &[i32], species: &str) -> Self {
        Self {
            offset: offset.to_vec(),
            species: species.to_string(),
        }
    }
}
