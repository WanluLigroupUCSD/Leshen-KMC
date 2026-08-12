// Clean-room implementation from OFFLATTICE_OTF_KMC_SPEC_V1 + ERRATA_1.
// Independently authored; no implementation source consulted.
use serde::de::{self, MapAccess, Visitor};
use serde::ser::SerializeStruct;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use serde_json::{Map, Value};
use std::fmt;

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Serialize, Deserialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
pub enum StatusCode {
    Ok,
    DiscoveryConvergedHeuristic,
    DuplicateEvent,
    SaddleNotFound,
    InvalidSaddle,
    SaddleWrongBasin,
    EndpointCollapsed,
    EnvironmentAmbiguous,
    BasinDisabled,
    DiscoveryIncomplete,
    RelaxNotConverged,
    EventApplicationFailed,
    CalculatorFailure,
    NonfiniteResult,
    InvalidInput,
    SchemaUnsupported,
    InvalidState,
    RateInvalid,
    DetailedBalanceViolation,
    CatalogConflict,
    CatalogIncompatible,
    AtomCountChangeUnsupported,
    NoEnabledEvent,
    ResourceLimit,
    OutputExists,
    CheckpointCorrupt,
    CheckpointIncompatible,
    Cancelled,
    InternalError,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Serialize, Deserialize)]
pub enum Severity {
    // E2-STATUS-002 fixes the exact severity vocabulary. It is NOT uniformly kebab-case:
    // `success-with-qualification` is hyphenated while `candidate reject`,
    // `pause/qualified`, `transaction fail`, `fatal in strict mode`, and
    // `terminal-success if requested, else fatal` are not. A blanket
    // `rename_all = "kebab-case"` silently rewrote 12 of the 29 table rows, which the
    // cross-language parity harness caught at byte offset 527 of the status table.
    #[serde(rename = "success")]
    Success,
    #[serde(rename = "success-with-qualification")]
    SuccessWithQualification,
    #[serde(rename = "candidate reject")]
    CandidateReject,
    #[serde(rename = "recoverable")]
    Recoverable,
    #[serde(rename = "pause/qualified")]
    PauseQualified,
    #[serde(rename = "transaction fail")]
    TransactionFail,
    #[serde(rename = "fatal")]
    Fatal,
    #[serde(rename = "fatal in strict mode")]
    FatalInStrictMode,
    #[serde(rename = "terminal-success if requested, else fatal")]
    TerminalSuccessIfRequestedElseFatal,
    #[serde(rename = "pause")]
    Pause,
}

impl StatusCode {
    pub const fn severity(self) -> Severity {
        use Severity::*;
        use StatusCode::*;
        match self {
            Ok => Success,
            DiscoveryConvergedHeuristic => SuccessWithQualification,
            DuplicateEvent | SaddleNotFound | InvalidSaddle | SaddleWrongBasin
            | EndpointCollapsed => CandidateReject,
            EnvironmentAmbiguous | BasinDisabled => Recoverable,
            DiscoveryIncomplete => PauseQualified,
            RelaxNotConverged | EventApplicationFailed | CalculatorFailure => TransactionFail,
            RateInvalid | DetailedBalanceViolation => FatalInStrictMode,
            NoEnabledEvent => TerminalSuccessIfRequestedElseFatal,
            ResourceLimit | Cancelled => Pause,
            NonfiniteResult | InvalidInput | SchemaUnsupported | InvalidState
            | CatalogConflict | CatalogIncompatible | AtomCountChangeUnsupported
            | OutputExists | CheckpointCorrupt | CheckpointIncompatible | InternalError => Fatal,
        }
    }

    pub const fn exit_code(self, requested_absorbing: bool, exploratory: bool) -> Option<i32> {
        use StatusCode::*;
        match self {
            Ok | DiscoveryConvergedHeuristic => Some(0),
            DiscoveryIncomplete if exploratory => Some(0),
            NoEnabledEvent if requested_absorbing => Some(0),
            InvalidInput | SchemaUnsupported | OutputExists => Some(64),
            RelaxNotConverged | EventApplicationFailed | NonfiniteResult | InvalidState
            | RateInvalid | DetailedBalanceViolation | CatalogConflict | CatalogIncompatible
            | AtomCountChangeUnsupported | NoEnabledEvent => Some(65),
            CalculatorFailure => Some(69),
            InternalError => Some(70),
            CheckpointCorrupt | CheckpointIncompatible => Some(74),
            DiscoveryIncomplete | ResourceLimit | Cancelled => Some(75),
            DuplicateEvent | SaddleNotFound | InvalidSaddle | SaddleWrongBasin | EndpointCollapsed
            | EnvironmentAmbiguous | BasinDisabled => None,
        }
    }
    pub const fn message(self)->&'static str{use StatusCode::*;match self{Ok=>"transaction committed",DiscoveryConvergedHeuristic=>"heuristic discovery criterion passed",DuplicateEvent=>"duplicate event rejected",SaddleNotFound=>"saddle not found",InvalidSaddle=>"saddle validation failed",SaddleWrongBasin=>"neither endpoint matches origin",EndpointCollapsed=>"both endpoints match origin",EnvironmentAmbiguous=>"environment identity ambiguous",BasinDisabled=>"basin acceleration disabled",DiscoveryIncomplete=>"discovery budget exhausted",RelaxNotConverged=>"relaxation not converged",EventApplicationFailed=>"event application failed",CalculatorFailure=>"calculator callback failed",NonfiniteResult=>"nonfinite value rejected",InvalidInput=>"input invalid",SchemaUnsupported=>"schema unsupported",InvalidState=>"state invalid",RateInvalid=>"rate invalid",DetailedBalanceViolation=>"detailed balance violated",CatalogConflict=>"catalog conflict",CatalogIncompatible=>"catalog incompatible",AtomCountChangeUnsupported=>"atom count change unsupported",NoEnabledEvent=>"no enabled event",ResourceLimit=>"resource limit reached",OutputExists=>"output exists",CheckpointCorrupt=>"checkpoint corrupt",CheckpointIncompatible=>"checkpoint incompatible",Cancelled=>"run cancelled",InternalError=>"internal error"}}
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ErrorContext {
    pub component: String,
    pub requirement_id: String,
    pub state_id: Option<String>,
    pub search_or_event_id: Option<String>,
    pub retryable: bool,
    #[serde(default)]
    pub details: Map<String, Value>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Status {
    pub status: StatusCode,
    pub severity: Severity,
    pub message: String,
    pub context: ErrorContext,
    pub cause: Option<Box<Status>>,
}

impl Status {
    pub fn new(status: StatusCode, message: impl Into<String>, context: ErrorContext) -> Self {
        Self { status, severity: status.severity(), message: message.into(), context, cause: None }
    }

    pub fn caused_by(mut self, cause: Status) -> Self {
        self.cause = Some(Box::new(cause));
        self
    }

    pub fn simple(status: StatusCode, component: &str, requirement_id: &str, message: &str) -> Self {
        Self::new(status, message, ErrorContext {
            component: component.to_owned(),
            requirement_id: requirement_id.to_owned(),
            state_id: None,
            search_or_event_id: None,
            retryable: false,
            details: Map::new(),
        })
    }
}

impl Serialize for Status {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut s = serializer.serialize_struct("Status", 6)?;
        s.serialize_field("causal_status", &self.cause.as_ref().map(|x|x.status))?;
        s.serialize_field("context", &self.context)?;
        s.serialize_field("exit_code", &self.status.exit_code(false,false))?;
        s.serialize_field("message", self.status.message())?;
        s.serialize_field("severity", &self.severity)?;
        s.serialize_field("status", &self.status)?;
        s.end()
    }
}

impl<'de> Deserialize<'de> for Status {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        enum Field { CausalStatus,Context,ExitCode,Message,Severity,Status }
        impl<'de> Deserialize<'de> for Field {
            fn deserialize<D2: Deserializer<'de>>(d: D2) -> Result<Self, D2::Error> {
                struct F;
                impl<'de> Visitor<'de> for F {
                    type Value = Field;
                    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { f.write_str("a status field") }
                    fn visit_str<E: de::Error>(self, v: &str) -> Result<Field, E> {
                        match v {
                            "causal_status"=>Ok(Field::CausalStatus),"context"=>Ok(Field::Context),"exit_code"=>Ok(Field::ExitCode),"message"=>Ok(Field::Message),"severity"=>Ok(Field::Severity),"status"=>Ok(Field::Status),_ => Err(E::unknown_field(v, FIELDS)),
                        }
                    }
                }
                d.deserialize_identifier(F)
            }
        }
        struct V;
        impl<'de> Visitor<'de> for V {
            type Value = Status;
            fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { f.write_str("strict Status") }
            fn visit_map<A: MapAccess<'de>>(self, mut a: A) -> Result<Status, A::Error> {
                let (mut causal,mut ctx,mut exit,mut msg,mut sev,mut code)=(None,None,None,None,None,None);
                while let Some(k) = a.next_key()? {
                    match k {
                        Field::CausalStatus=>set_once(&mut causal,a.next_value::<Option<StatusCode>>()?,"causal_status")?,
                        Field::Context => set_once(&mut ctx, a.next_value()?, "context")?,
                        Field::ExitCode=>set_once(&mut exit,a.next_value::<Option<i32>>()?,"exit_code")?,
                        Field::Message=>set_once(&mut msg,a.next_value::<String>()?,"message")?,
                        Field::Severity=>set_once(&mut sev,a.next_value()?,"severity")?,
                        Field::Status=>set_once(&mut code,a.next_value::<StatusCode>()?,"status")?,
                    }
                }
                let status = code.ok_or_else(|| de::Error::missing_field("status"))?;
                let severity = sev.ok_or_else(|| de::Error::missing_field("severity"))?;
                if severity != status.severity() { return Err(de::Error::custom("severity does not match status")); }
                let message=msg.ok_or_else(||de::Error::missing_field("message"))?;if message!=status.message(){return Err(de::Error::custom("message does not match status"))}let exit_code=exit.ok_or_else(||de::Error::missing_field("exit_code"))?;if exit_code!=status.exit_code(false,false){return Err(de::Error::custom("exit code does not match status"))}let cause=causal.flatten().map(|c|Box::new(Status::simple(c,"causal-status","E2-STATUS-004",c.message())));Ok(Status{status,severity,message,context:ctx.ok_or_else(||de::Error::missing_field("context"))?,cause})
            }
        }
        const FIELDS:&[&str]=&["causal_status","context","exit_code","message","severity","status"];
        deserializer.deserialize_struct("Status", FIELDS, V)
    }
}

fn set_once<T, E: de::Error>(slot: &mut Option<T>, value: T, name: &'static str) -> Result<(), E> {
    if slot.is_some() { return Err(E::duplicate_field(name)); }
    *slot = Some(value);
    Ok(())
}
