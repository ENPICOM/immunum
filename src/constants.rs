use phf::{phf_map, Map};
use std::collections::HashMap;

pub const GAP_PEN_START: f64 = 1.0;
pub const GAP_PEN_END: f64 = 1.0;
pub const MATCH_CP_MULTIPLIER: f64 = 8.0; // Multiplier of match score

#[derive(Clone)]
pub struct ScoringParams {
    pub gap_pen_cp: f64,
    pub gap_pen_fr: f64,
    pub gap_pen_ip: f64,
    pub gap_pen_op: f64,
    pub gap_pen_cdr: f64,
    pub gap_pen_other: f64,
    pub cdr_increase: f64, // increase per position away from insertion position CDR

    pub pen_leap_insertion_point_imgt: f64,  // better name
    pub pen_leap_insertion_point_kabat: f64, // better name
}

pub fn get_scoring_params() -> ScoringParams {
    ScoringParams::default()
    // use below if you want to chance params
    // ScoringParams {gap_pen_start: 10.0, ..Default::default()
}

impl Default for ScoringParams {
    fn default() -> ScoringParams {
        ScoringParams {
            gap_pen_cp: 55.0,
            gap_pen_fr: 26.0,
            gap_pen_ip: 1.5,
            gap_pen_op: 1.0,
            gap_pen_cdr: 2.5,
            gap_pen_other: 11.0,
            cdr_increase: 0.5, // increase per position away from insertion position CDR

            pen_leap_insertion_point_imgt: 1.0, // better name 6 best until now
            pen_leap_insertion_point_kabat: 10.0, // better name
        }
    }
}
// For prefiltering, how close the identity must be to the highest found in order to run the scheme
pub const WITHIN_IDENTITY_RANGE: f64 = 0.25;

pub const MINIMAL_CHAIN_IDENTITY: f64 = 0.7;

// Minimal chain length, minimal length for sequence to continue search for chain
pub const MINIMAL_CHAIN_LENGTH: i32 = 60;

pub const PRE_FILTER_TERMINAL_LENGTH: u8 = 10;

// Indexes for query and consensus gap columns in scoring matrix
pub const QUERY_GAP_COLUMN: usize = 23;
pub const CONSENSUS_GAP_COLUMN: usize = 24;
pub mod traceback_directions {
    pub const FROM_DIAG: u8 = 0;
    pub const FROM_LEFT: u8 = 1;
    pub const FROM_TOP: u8 = 2;
    pub const PERFECT_MATCH: u8 = 3;
}

// IMGT cdr insertion positions;
pub mod insertion_points {
    pub const CDR1_IMGT: u32 = 33;
    pub const CDR2_IMGT: u32 = 61;
    pub const CDR3_IMGT: u32 = 111;

    // KABAT cdr insertion positions
    pub const CDR1_KABAT_HEAVY: u32 = 35;
    pub const CDR1_KABAT_LIGHT: u32 = 28;
    pub const CDR2_KABAT_HEAVY: u32 = 52;
    pub const CDR2_KABAT_LIGHT: u32 = 52;
    pub const CDR3_KABAT_HEAVY: u32 = 100;
    pub const CDR3_KABAT_LIGHT: u32 = 95;
}
// ALPHABET for numbering
pub const ALPHABET: [char; 26] = [
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
    'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
];

// TODO get allowed amino acids
//pub const ALLOWED_AMINO_ACIDS: Vec<&str> = BLOSUM62.keys().collect();

pub const ACCEPTED_RESIDUES: [u8; 23] = [
    b'A', b'B', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R',
    b'S', b'T', b'V', b'W', b'X', b'Y', b'Z',
];

pub static ENCODED_RESIDUES_MAP: Map<u8, u8> = phf_map! {
        b'A'=> 0,
        b'B'=> 1,
        b'C'=> 2,
        b'D'=> 3,
        b'E'=> 4,
        b'F'=> 5,
        b'G'=> 6,
        b'H'=> 7,
        b'I'=> 8,
        b'K'=> 9,
        b'L'=> 10,
        b'M'=> 11,
        b'N'=> 12,
        b'P'=> 13,
        b'Q'=> 14,
        b'R'=> 15,
        b'S'=> 16,
        b'T'=> 17,
        b'V'=> 18,
        b'W'=> 19,
        b'X'=> 20,
        b'Y'=> 21,
        b'Z'=> 22,
};
pub static BLOSUM62: Map<&[u8; 2], i32> = phf_map! {
        b"BN" => 3,
        b"LW" => -2,
        b"GG" => 6,
        b"SX" => 0,
        b"DX" => -1,
        b"GK" => -2,
        b"ES" => 0,
        b"MX" => -1,
        b"EY" => -2,
        b"RW" => -3,
        b"IR" => -3,
        b"XZ" => -1,
        b"EH" => 0,
        b"MV" => 1,
        b"NR" => 0,
        b"DI" => -3,
        b"DF" => -3,
        b"CW" => -2,
        b"AN" => -2,
        b"QW" => -2,
        b"LQ" => -2,
        b"NS" => 1,
        b"KZ" => 1,
        b"NV" => -3,
        b"NQ" => 0,
        b"KM" => -1,
        b"HV" => -3,
        b"EG" => -2,
        b"LS" => -2,
        b"PR" => -2,
        b"AD" => -2,
        b"CS" => -1,
        b"DE" => 2,
        b"GY" => -3,
        b"PW" => -4,
        b"XX" => -1,
        b"LZ" => -3,
        b"AQ" => -1,
        b"VY" => -1,
        b"AW" => -3,
        b"DG" => -1,
        b"PX" => -2,
        b"DK" => -1,
        b"NT" => 0,
        b"FY" => 3,
        b"WW" => 11,
        b"MZ" => -1,
        b"DL" => -4,
        b"MR" => -1,
        b"KY" => -2,
        b"EF" => -3,
        b"EM" => -2,
        b"SS" => 4,
        b"CX" => -2,
        b"LY" => -1,
        b"HR" => 0,
        b"PP" => 7,
        b"CK" => -3,
        b"AS" => 1,
        b"IP" => -3,
        b"QQ" => 5,
        b"IL" => 2,
        b"FP" => -4,
        b"AB" => -2,
        b"NZ" => 0,
        b"MQ" => 0,
        b"IV" => 3,
        b"CQ" => -3,
        b"HI" => -3,
        b"DZ" => 1,
        b"PZ" => -1,
        b"WY" => 2,
        b"GT" => -2,
        b"BP" => -2,
        b"AP" => -1,
        b"CD" => -3,
        b"HY" => 2,
        b"VX" => -1,
        b"BB" => 4,
        b"FZ" => -3,
        b"LM" => 2,
        b"FG" => -3,
        b"MS" => -1,
        b"GM" => -3,
        b"QZ" => 3,
        b"QS" => 0,
        b"AX" => 0,
        b"TV" => 0,
        b"FW" => 1,
        b"HS" => -1,
        b"NX" => -1,
        b"BQ" => 0,
        b"AK" => -1,
        b"IQ" => -3,
        b"WX" => -2,
        b"NN" => 6,
        b"TW" => -2,
        b"DP" => -1,
        b"BC" => -3,
        b"CI" => -1,
        b"KV" => -2,
        b"XY" => -1,
        b"KR" => 2,
        b"RZ" => 0,
        b"EW" => -3,
        b"ET" => -1,
        b"BR" => -1,
        b"LR" => -2,
        b"QR" => 1,
        b"FX" => -1,
        b"ST" => 1,
        b"BD" => 4,
        b"AZ" => -1,
        b"MN" => -2,
        b"DV" => -3,
        b"AF" => -2,
        b"EX" => -1,
        b"FH" => -1,
        b"AM" => -1,
        b"KQ" => 1,
        b"SZ" => 0,
        b"GX" => -1,
        b"VV" => 4,
        b"DW" => -4,
        b"HX" => -1,
        b"FS" => -2,
        b"LX" => -1,
        b"BS" => 0,
        b"GS" => 0,
        b"MP" => -2,
        b"MY" => -1,
        b"DH" => -1,
        b"BE" => 1,
        b"BZ" => 1,
        b"EI" => -3,
        b"EV" => -2,
        b"TX" => 0,
        b"RX" => -1,
        b"RR" => 5,
        b"TZ" => -1,
        b"DY" => -3,
        b"VW" => -3,
        b"FL" => 0,
        b"CT" => -1,
        b"QX" => -1,
        b"BT" => -1,
        b"KN" => 0,
        b"HT" => -2,
        b"IY" => -1,
        b"FQ" => -3,
        b"IT" => -1,
        b"QT" => -1,
        b"LP" => -3,
        b"AR" => -1,
        b"BF" => -3,
        b"CZ" => -3,
        b"HM" => -2,
        b"FV" => -1,
        b"CF" => -2,
        b"LL" => 4,
        b"CM" => -1,
        b"CR" => -3,
        b"DD" => 6,
        b"ER" => 0,
        b"PV" => -2,
        b"DS" => 0,
        b"EE" => 5,
        b"GW" => -2,
        b"CP" => -3,
        b"FR" => -3,
        b"BG" => -1,
        b"CC" => 9,
        b"GI" => -4,
        b"GV" => -3,
        b"KW" => -3,
        b"GN" => 0,
        b"IN" => -3,
        b"VZ" => -2,
        b"AA" => 4,
        b"QV" => -2,
        b"FK" => -3,
        b"AT" => 0,
        b"BV" => -3,
        b"KL" => -2,
        b"LN" => -3,
        b"NY" => -2,
        b"FF" => 6,
        b"GL" => -4,
        b"BH" => 0,
        b"EZ" => 4,
        b"DQ" => 0,
        b"BX" => -1,
        b"WZ" => -3,
        b"KS" => 0,
        b"KX" => -1,
        b"RV" => -3,
        b"EK" => 1,
        b"AI" => -1,
        b"HP" => -2,
        b"BW" => -4,
        b"KK" => 5,
        b"CH" => -3,
        b"EN" => 0,
        b"QY" => -1,
        b"HH" => 8,
        b"BI" => -3,
        b"AC" => 0,
        b"II" => 4,
        b"AV" => 0,
        b"IW" => -3,
        b"FT" => -2,
        b"SV" => -2,
        b"TT" => 5,
        b"FM" => 0,
        b"EL" => -3,
        b"MM" => 5,
        b"GZ" => -2,
        b"DR" => -2,
        b"DM" => -3,
        b"HW" => -2,
        b"CG" => -3,
        b"RS" => -1,
        b"IS" => -2,
        b"PQ" => -1,
        b"AY" => -2,
        b"IX" => -1,
        b"AE" => -1,
        b"BY" => -3,
        b"IK" => -3,
        b"AH" => -2,
        b"GP" => -2,
        b"FN" => -3,
        b"HN" => 1,
        b"BK" => 0,
        b"CV" => -1,
        b"LT" => -1,
        b"KP" => -1,
        b"SW" => -3,
        b"DT" => -1,
        b"MT" => -1,
        b"NP" => -2,
        b"HK" => -1,
        b"RT" => -1,
        b"RY" => -2,
        b"CL" => -1,
        b"BL" => -4,
        b"YZ" => -2,
        b"NW" => -4,
        b"AG" => 0,
        b"PS" => -1,
        b"EQ" => 2,
        b"CN" => -3,
        b"HQ" => 0,
        b"DN" => 1,
        b"CY" => -2,
        b"HL" => -3,
        b"CE" => -4,
        b"HZ" => 0,
        b"GH" => -2,
        b"EP" => -1,
        b"SY" => -2,
        b"GR" => -2,
        b"BM" => -3,
        b"ZZ" => 4,
        b"MW" => -1,
        b"TY" => -2,
        b"PY" => -3,
        b"YY" => 7,
        b"KT" => -1,
        b"IZ" => -3,
        b"PT" => -1,
        b"LV" => 1,
        b"FI" => 0,
        b"GQ" => -2,
        b"AL" => -1,
        b"IM" => 1
};

// Embedded consensus data at compile time
const IMGT_CONSENSUS_H_DATA: &str = include_str!("../resources/consensus/IMGT_CONSENSUS_H.txt");
const IMGT_CONSENSUS_K_DATA: &str = include_str!("../resources/consensus/IMGT_CONSENSUS_K.txt");
const IMGT_CONSENSUS_L_DATA: &str = include_str!("../resources/consensus/IMGT_CONSENSUS_L.txt");
const KABAT_CONSENSUS_H_DATA: &str = include_str!("../resources/consensus/KABAT_CONSENSUS_H.txt");
const KABAT_CONSENSUS_K_DATA: &str = include_str!("../resources/consensus/KABAT_CONSENSUS_K.txt");
const KABAT_CONSENSUS_L_DATA: &str = include_str!("../resources/consensus/KABAT_CONSENSUS_L.txt");

/// Parse consensus data from embedded string at compile time
pub fn parse_consensus_data(data: &str) -> HashMap<u32, Vec<u8>> {
    let mut consensus_aas: HashMap<u32, Vec<u8>> = HashMap::new();
    let total_lines = data.lines().count();
    
    // Skip first and last line
    for line in data.lines().skip(1).take(total_lines - 2) {
        let split_line: Vec<&str> = line.split(',').collect();
        if let Ok(position) = split_line[0].parse::<u32>() {
            let amino_acids = split_line[1..].join("").into_bytes();
            consensus_aas.insert(position, amino_acids);
        }
    }
    consensus_aas
}

/// Get IMGT Heavy consensus data
pub fn get_imgt_heavy_consensus() -> HashMap<u32, Vec<u8>> {
    parse_consensus_data(IMGT_CONSENSUS_H_DATA)
}

/// Get IMGT Kappa consensus data
pub fn get_imgt_kappa_consensus() -> HashMap<u32, Vec<u8>> {
    parse_consensus_data(IMGT_CONSENSUS_K_DATA)
}

/// Get IMGT Lambda consensus data
pub fn get_imgt_lambda_consensus() -> HashMap<u32, Vec<u8>> {
    parse_consensus_data(IMGT_CONSENSUS_L_DATA)
}

/// Get KABAT Heavy consensus data
pub fn get_kabat_heavy_consensus() -> HashMap<u32, Vec<u8>> {
    parse_consensus_data(KABAT_CONSENSUS_H_DATA)
}

/// Get KABAT Kappa consensus data
pub fn get_kabat_kappa_consensus() -> HashMap<u32, Vec<u8>> {
    parse_consensus_data(KABAT_CONSENSUS_K_DATA)
}

/// Get KABAT Lambda consensus data
pub fn get_kabat_lambda_consensus() -> HashMap<u32, Vec<u8>> {
    parse_consensus_data(KABAT_CONSENSUS_L_DATA)
}
