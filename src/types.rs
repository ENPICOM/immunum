#![allow(clippy::upper_case_acronyms)]

use clap::ValueEnum;

/// Numbering schemes for immunoglobulin sequences
#[derive(Clone, Debug, ValueEnum)]
pub enum Scheme {
    /// IMGT numbering scheme
    IMGT,
    /// Kabat numbering scheme
    KABAT,
}


/// Immunoglobulin and T-cell receptor chain types
#[derive(Clone, Debug, ValueEnum)]
pub enum Chain {
    // IG Heavy chain variants
    #[value(alias = "Heavy", alias = "H")]
    IGH,

    // IG Kappa chain variants
    #[value(alias = "Kappa", alias = "K")]
    IGK,

    // IG Lambda chain variants
    #[value(alias = "Lambda", alias = "L")]
    IGL,

    // T-cell receptor Alpha chain variants
    #[value(alias = "Alpha", alias = "A")]
    TRA,

    // T-cell receptor Beta chain variants
    #[value(alias = "Beta", alias = "B")]
    TRB,

    // T-cell receptor Gamma chain variants
    #[value(alias = "Gamma", alias = "G")]
    TRG,

    // T-cell receptor Delta chain variants
    #[value(alias = "Delta", alias = "D")]
    TRD,
}