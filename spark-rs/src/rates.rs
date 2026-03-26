/// Rate constant calculations: Arrhenius, TST, Hertz-Knudsen, Butler-Volmer.

use crate::units::*;
use crate::model::{RateExpr, Parameter};
use std::collections::HashMap;

/// Evaluate a rate expression given current parameters.
pub fn evaluate_rate_expression(expr: &RateExpr, params: &[Parameter]) -> f64 {
    match expr {
        RateExpr::Constant(v) => *v,
        RateExpr::Expression(s) => evaluate_string_expr(s, params),
    }
}

/// Lazily-initialized base variable map (constants + molecular masses).
/// Avoids rebuilding ~30 HashMap entries on every rate evaluation.
use std::sync::OnceLock;

static BASE_VARS: OnceLock<HashMap<String, f64>> = OnceLock::new();

fn get_base_vars() -> &'static HashMap<String, f64> {
    BASE_VARS.get_or_init(|| {
        let mut vars = HashMap::with_capacity(32);
        // Physical constants
        vars.insert("kB".into(), KB);
        vars.insert("h".into(), H_PLANCK);
        vars.insert("eV".into(), EV);
        vars.insert("bar".into(), BAR);
        vars.insert("angstrom".into(), ANGSTROM);
        vars.insert("umass".into(), UMASS);
        vars.insert("pi".into(), PI);
        vars.insert("NA".into(), NA);
        vars.insert("e".into(), E_CHARGE);
        // Molecular masses
        for name in &["H", "H2", "N", "N2", "O2", "CO", "CO2", "H2O", "NH3",
                       "N2H4", "NNH", "Mo", "W"] {
            vars.insert(format!("m_{}", name), molecular_mass(name));
        }
        vars
    })
}

/// Evaluate a string-based rate expression.
///
/// Uses cached base constants; only merges user parameters per call.
fn evaluate_string_expr(expr: &str, params: &[Parameter]) -> f64 {
    let base = get_base_vars();
    let mut vars = base.clone();

    // User parameters (typically 3-5 entries)
    for p in params {
        vars.insert(p.name.clone(), p.value);
    }

    // Compute beta = 1/(kB*T) if T is defined
    if let Some(&t) = vars.get("T") {
        if t > 0.0 {
            vars.insert("beta".into(), 1.0 / (KB * t));
        }
    }

    match eval_expr(expr, &vars) {
        Ok(v) if v.is_finite() => v.max(0.0),
        Ok(v) => {
            eprintln!("Warning: rate '{}' evaluated to {}", expr, v);
            0.0
        }
        Err(e) => {
            eprintln!("Error evaluating '{}': {}", expr, e);
            0.0
        }
    }
}

// ─── Simple recursive-descent expression parser ───────────────────────────
//
// Supports: +, -, *, /, ** (power), unary -, parentheses,
//           exp(), log(), sqrt(), sin(), cos(), abs(), pow(a,b)
//           and variable lookup.

struct Parser<'a> {
    input: &'a [u8],
    pos: usize,
    vars: &'a HashMap<String, f64>,
}

impl<'a> Parser<'a> {
    fn new(input: &'a str, vars: &'a HashMap<String, f64>) -> Self {
        Self { input: input.as_bytes(), pos: 0, vars }
    }

    fn skip_ws(&mut self) {
        while self.pos < self.input.len() && self.input[self.pos].is_ascii_whitespace() {
            self.pos += 1;
        }
    }

    fn peek(&mut self) -> Option<u8> {
        self.skip_ws();
        self.input.get(self.pos).copied()
    }

    fn consume(&mut self, ch: u8) -> bool {
        self.skip_ws();
        if self.pos < self.input.len() && self.input[self.pos] == ch {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    // expr = term (('+' | '-') term)*
    fn parse_expr(&mut self) -> Result<f64, String> {
        let mut result = self.parse_term()?;
        loop {
            if self.consume(b'+') {
                result += self.parse_term()?;
            } else if self.consume(b'-') {
                result -= self.parse_term()?;
            } else {
                break;
            }
        }
        Ok(result)
    }

    // term = power (('*' | '/') power)*
    fn parse_term(&mut self) -> Result<f64, String> {
        let mut result = self.parse_power()?;
        loop {
            if self.consume(b'*') {
                if self.consume(b'*') {
                    // ** operator
                    let exp = self.parse_unary()?;
                    result = result.powf(exp);
                } else {
                    result *= self.parse_power()?;
                }
            } else if self.consume(b'/') {
                result /= self.parse_power()?;
            } else {
                break;
            }
        }
        Ok(result)
    }

    // power = unary ('**' unary)?
    fn parse_power(&mut self) -> Result<f64, String> {
        let base = self.parse_unary()?;
        self.skip_ws();
        if self.pos + 1 < self.input.len()
            && self.input[self.pos] == b'*'
            && self.input[self.pos + 1] == b'*'
        {
            self.pos += 2;
            let exp = self.parse_unary()?;
            Ok(base.powf(exp))
        } else {
            Ok(base)
        }
    }

    // unary = '-' unary | '+' unary | atom
    fn parse_unary(&mut self) -> Result<f64, String> {
        if self.consume(b'-') {
            Ok(-self.parse_unary()?)
        } else if self.consume(b'+') {
            self.parse_unary()
        } else {
            self.parse_atom()
        }
    }

    // atom = number | '(' expr ')' | func '(' args ')' | variable
    fn parse_atom(&mut self) -> Result<f64, String> {
        self.skip_ws();
        if self.pos >= self.input.len() {
            return Err("Unexpected end of expression".into());
        }

        // Parenthesized expression
        if self.input[self.pos] == b'(' {
            self.pos += 1;
            let v = self.parse_expr()?;
            if !self.consume(b')') {
                return Err("Missing ')'".into());
            }
            return Ok(v);
        }

        // Number
        if self.input[self.pos].is_ascii_digit() || self.input[self.pos] == b'.' {
            return self.parse_number();
        }

        // Identifier (variable or function)
        let ident = self.parse_ident()?;

        // Check if it's a function call
        if self.consume(b'(') {
            let result = match ident.as_str() {
                "exp" => {
                    let arg = self.parse_expr()?;
                    Ok(arg.exp())
                }
                "log" | "ln" => {
                    let arg = self.parse_expr()?;
                    Ok(arg.ln())
                }
                "sqrt" => {
                    let arg = self.parse_expr()?;
                    Ok(arg.sqrt())
                }
                "sin" => {
                    let arg = self.parse_expr()?;
                    Ok(arg.sin())
                }
                "cos" => {
                    let arg = self.parse_expr()?;
                    Ok(arg.cos())
                }
                "abs" => {
                    let arg = self.parse_expr()?;
                    Ok(arg.abs())
                }
                "pow" => {
                    let base = self.parse_expr()?;
                    if !self.consume(b',') {
                        return Err("pow() needs two arguments".into());
                    }
                    let exp = self.parse_expr()?;
                    Ok(base.powf(exp))
                }
                _ => Err(format!("Unknown function: {}", ident)),
            }?;
            if !self.consume(b')') {
                return Err(format!("Missing ')' after {}()", ident));
            }
            return Ok(result);
        }

        // Variable lookup
        if let Some(&val) = self.vars.get(&ident) {
            Ok(val)
        } else {
            Err(format!("Unknown variable: {}", ident))
        }
    }

    fn parse_number(&mut self) -> Result<f64, String> {
        let start = self.pos;
        while self.pos < self.input.len()
            && (self.input[self.pos].is_ascii_digit()
                || self.input[self.pos] == b'.'
                || self.input[self.pos] == b'e'
                || self.input[self.pos] == b'E'
                || ((self.input[self.pos] == b'+' || self.input[self.pos] == b'-')
                    && self.pos > start
                    && (self.input[self.pos - 1] == b'e' || self.input[self.pos - 1] == b'E')))
        {
            self.pos += 1;
        }
        let s = std::str::from_utf8(&self.input[start..self.pos])
            .map_err(|_| "Invalid UTF-8 in number")?;
        s.parse::<f64>().map_err(|e| format!("Bad number '{}': {}", s, e))
    }

    fn parse_ident(&mut self) -> Result<String, String> {
        let start = self.pos;
        while self.pos < self.input.len()
            && (self.input[self.pos].is_ascii_alphanumeric() || self.input[self.pos] == b'_')
        {
            self.pos += 1;
        }
        if self.pos == start {
            return Err(format!(
                "Expected identifier at position {}: '{}'",
                self.pos,
                String::from_utf8_lossy(&self.input[self.pos..])
            ));
        }
        Ok(String::from_utf8_lossy(&self.input[start..self.pos]).to_string())
    }
}

fn eval_expr(expr: &str, vars: &HashMap<String, f64>) -> Result<f64, String> {
    let mut parser = Parser::new(expr, vars);
    let result = parser.parse_expr()?;
    parser.skip_ws();
    if parser.pos < parser.input.len() {
        return Err(format!(
            "Unexpected character at position {}: '{}'",
            parser.pos,
            parser.input[parser.pos] as char
        ));
    }
    Ok(result)
}

// ─── Convenience rate functions ─────────────────────────────────────────

/// Arrhenius: k = A * exp(-Ea_eV * eV / (kB * T))
pub fn arrhenius(a: f64, ea_ev: f64, t: f64) -> f64 {
    a * (-ea_ev * EV / (KB * t)).exp()
}

/// Transition state theory: k = (kB*T/h) * exp(-Ea_eV * eV / (kB*T))
pub fn tst_rate(ea_ev: f64, t: f64) -> f64 {
    (KB * t / H_PLANCK) * (-ea_ev * EV / (KB * t)).exp()
}

/// Hertz-Knudsen: k = p * A_site / sqrt(2π * m_kg * kB * T)
pub fn hertz_knudsen(p_pa: f64, t: f64, m_kg: f64, a_site_m2: f64) -> f64 {
    p_pa * a_site_m2 / (2.0 * PI * m_kg * KB * t).sqrt()
}

/// BEP-modified rate: k = k_base * exp(-alpha * delta_delta_H * eV / (kB*T))
pub fn bep_modified_rate(base_rate: f64, delta_delta_h: f64, alpha: f64, t: f64) -> f64 {
    if t <= 0.0 || base_rate <= 0.0 { return 0.0; }
    let delta_ea = alpha * delta_delta_h;
    base_rate * (-delta_ea * EV / (KB * t)).exp()
}

/// Lateral-modified rate: k = k_base * exp(+E_lateral * eV / (kB*T))
/// Positive E_lateral (repulsive) -> higher rate (destabilized adsorbate).
pub fn lateral_modified_rate(base_rate: f64, e_lateral: f64, t: f64) -> f64 {
    if t <= 0.0 || base_rate <= 0.0 { return 0.0; }
    base_rate * (e_lateral * EV / (KB * t)).exp()
}

/// Electrochemical Butler-Volmer PCET rate:
/// k = (kB*T/h) * exp(-(Ea + beta_bv * (U - U0)) * eV / (kB*T))
pub fn electrochemical_rate(ea_ev: f64, t: f64, u: f64, u0: f64, beta_bv: f64) -> f64 {
    let ea_eff = (ea_ev + beta_bv * (u - u0)).max(0.0);
    (KB * t / H_PLANCK) * (-ea_eff * EV / (KB * t)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eval_basic() {
        let vars = HashMap::new();
        assert!((eval_expr("2+3", &vars).unwrap() - 5.0).abs() < 1e-10);
        assert!((eval_expr("2*3+1", &vars).unwrap() - 7.0).abs() < 1e-10);
        assert!((eval_expr("2*(3+1)", &vars).unwrap() - 8.0).abs() < 1e-10);
        assert!((eval_expr("2**3", &vars).unwrap() - 8.0).abs() < 1e-10);
        assert!((eval_expr("-5+3", &vars).unwrap() - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_eval_functions() {
        let vars = HashMap::new();
        assert!((eval_expr("exp(0)", &vars).unwrap() - 1.0).abs() < 1e-10);
        assert!((eval_expr("sqrt(4)", &vars).unwrap() - 2.0).abs() < 1e-10);
        assert!((eval_expr("log(1)", &vars).unwrap()).abs() < 1e-10);
    }

    #[test]
    fn test_eval_variables() {
        let mut vars = HashMap::new();
        vars.insert("T".into(), 300.0);
        vars.insert("kB".into(), KB);
        assert!((eval_expr("kB*T", &vars).unwrap() - KB * 300.0).abs() < 1e-30);
    }

    #[test]
    fn test_tst() {
        let k = tst_rate(0.5, 600.0);
        let expected = (KB * 600.0 / H_PLANCK) * (-0.5 * EV / (KB * 600.0)).exp();
        assert!((k - expected).abs() / expected < 1e-10);
    }
}
