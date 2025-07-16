use crate::constants::{insertion_points, ALPHABET};
use crate::types::Scheme;

pub(crate) fn name_insertions(numbering: &mut Vec<String>, scheme: &Scheme) {
    // Function that assigns correct labels to insertions
    let original_numbering_length = numbering.len();

    // For IMGT reverse numbering
    let reverse_numbering_positions = [
        insertion_points::CDR1_IMGT,
        insertion_points::CDR2_IMGT,
        insertion_points::CDR3_IMGT,
    ];

    let mut latest_number: Option<String> = None;
    let mut first_number: Option<String> = None;
    let mut gaps_dict: std::collections::HashMap<String, Vec<usize>> =
        std::collections::HashMap::new();

    // find first and last numbered position, also all gaps
    for (i, item) in numbering.iter().enumerate() {
        if item != "-" {
            // number
            if first_number.is_none() {
                first_number = Some(item.to_string());
            }
            latest_number = Some(item.to_string()); // set to latest
        } else if item == "-" {
            // gap
            if first_number.is_some() && latest_number.is_some() {
                // inside antibody
                let latest = latest_number.as_ref().unwrap();
                gaps_dict.entry(latest.to_string()).or_default().push(i); // add i
            }
        }
    }

    // remove last number (end of sequence) from gaps
    if let Some(latest) = latest_number {
        gaps_dict.remove(&latest);
    }

    // number non special gaps
    for (gap_position, gaps) in gaps_dict.iter() {
        let mut numbered_gap: Vec<String> = Vec::new();

        if *scheme != Scheme::IMGT
            || !reverse_numbering_positions.contains(&gap_position.parse::<u32>().unwrap_or(0))
        {
            // create list of correct named gap positions
            let mut alphabet_index: usize = 0;
            let mut base: String = String::new();
            let mut base_index: usize = 0;

            for _ in 0..gaps.len() {
                let addition = format!("{}{}", base, ALPHABET[alphabet_index]);
                let insertion_name = format!("{}{}", gap_position, addition);
                numbered_gap.push(insertion_name);
                alphabet_index += 1;

                if alphabet_index == 26 {
                    // add extra letter if gap longer than 26
                    base = ALPHABET[base_index].to_string();
                    base_index += 1;
                    alphabet_index = 0;
                }
            }
        } else {
            // special reverse numbering
            numbered_gap =
                imgt_reverse_numbering(gap_position.parse::<u32>().unwrap(), gaps.len(), true);
        }

        // insert this list in place of gap
        let start = gaps[0];
        let end = gaps[gaps.len() - 1] + 1;

        numbering.splice(start..end, numbered_gap);
    }

    if numbering.len() != original_numbering_length {
        panic!("Placing insertions caused numbering length to change");
    }
}

fn imgt_reverse_numbering(position: u32, insertion_length: usize, letters: bool) -> Vec<String> {
    // Function to get a list of insertion names according to imgt scheme
    if position != insertion_points::CDR1_IMGT
        && position != insertion_points::CDR2_IMGT
        && position != insertion_points::CDR3_IMGT
    {
        panic!("Trying to reverse number a wrong position {}", position);
    }

    let mut numbered_gap: Vec<String> = Vec::new();
    let mut alphabet_index: usize = 0;
    let mut plus_one: bool = true;
    let mut base: String = String::new();
    let mut base_index: usize = 0;

    for i in 0..insertion_length {
        let insertion_name = if letters {
            let addition = format!("{}{}", base, ALPHABET[alphabet_index]);
            format!("{}{}", position + if plus_one { 1 } else { 0 }, addition)
        } else {
            // use original .x format instead
            format!(
                "{}.{}",
                position + if plus_one { 1 } else { 0 },
                alphabet_index + 1
            )
        };

        numbered_gap.insert(numbered_gap.len() / 2, insertion_name);
        plus_one = !plus_one; // goes from e.g. 111 to 112 and back

        if i % 2 == 1 {
            // every 2 insertions go to next alphabet character
            alphabet_index += 1;
        }

        if alphabet_index == 26 {
            // handle very long insertions
            base = ALPHABET[base_index].to_string();
            base_index += 1;
            alphabet_index = 0;
        }
    }

    numbered_gap
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_insertion_numbering() {
        let mut numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            '2'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '3'.to_string(),
            '4'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        let correct_numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            '2'.to_string(),
            "2A".to_string(),
            "2B".to_string(),
            "2C".to_string(),
            '3'.to_string(),
            '4'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        let scheme = Scheme::IMGT;
        name_insertions(&mut numbering, &scheme);
        assert_eq!(numbering, correct_numbering);
        let mut numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            '2'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '3'.to_string(),
            '4'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        let scheme = Scheme::KABAT;
        name_insertions(&mut numbering, &scheme);
        assert_eq!(numbering, correct_numbering);
    }

    #[test]
    fn test_reverse_numbering() {
        let mut numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "33".to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            "34".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        let correct_numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "33".to_string(),
            "33A".to_string(),
            "33B".to_string(),
            "34B".to_string(),
            "34A".to_string(),
            "34".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        name_insertions(&mut numbering, &Scheme::IMGT);
        assert_eq!(numbering, correct_numbering);

        let mut numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "61".to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            "70".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        let correct_numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "61".to_string(),
            "61A".to_string(),
            "61B".to_string(),
            "62B".to_string(),
            "62A".to_string(),
            "70".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        name_insertions(&mut numbering, &Scheme::IMGT);
        assert_eq!(numbering, correct_numbering);

        let mut numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "111".to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            "114".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        let correct_numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "111".to_string(),
            "111A".to_string(),
            "111B".to_string(),
            "112B".to_string(),
            "112A".to_string(),
            "114".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        name_insertions(&mut numbering, &Scheme::IMGT);
        assert_eq!(numbering, correct_numbering);

        let mut numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "33".to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            '-'.to_string(),
            "34".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        let correct_numbering = vec![
            '-'.to_string(),
            '-'.to_string(),
            '1'.to_string(),
            "33".to_string(),
            "33A".to_string(),
            "33B".to_string(),
            "33C".to_string(),
            "33D".to_string(),
            "34".to_string(),
            '-'.to_string(),
            '-'.to_string(),
        ];
        name_insertions(&mut numbering, &Scheme::KABAT);
        assert_eq!(numbering, correct_numbering);
    }

    #[test]
    fn test_long_reverse_insertion() {
        let mut test_vec: Vec<String> = vec!["-".to_string(); 33]; // 10 dashes
                                                                   // Add numbers 1 to 10
        for i in 10..=33 {
            test_vec.push(i.to_string());
        }
        // Add 50 dashes
        test_vec.extend(vec!["-".to_string(); 100]);
        // Add numbers 11 to 15
        for i in 34..=50 {
            test_vec.push(i.to_string());
        }
        // Add more dashes
        test_vec.extend(vec!["-".to_string(); 10]);
        name_insertions(&mut test_vec, &Scheme::IMGT);
        // check reversal point to see if numbering goes well
        assert!(test_vec.contains(&String::from("33AX")));
        assert!(test_vec.contains(&String::from("34AX")));
    }

    #[test]
    fn test_long_insertion() {
        let mut test_vec: Vec<String> = vec!["-".to_string(); 33]; // 10 dashes
                                                                   // Add numbers 1 to 10
        for i in 1..=10 {
            test_vec.push(i.to_string());
        }
        // Add 50 dashes
        test_vec.extend(vec!["-".to_string(); 100]);
        // Add numbers 11 to 15
        for i in 11..=15 {
            test_vec.push(i.to_string());
        }
        // Add more dashes
        test_vec.extend(vec!["-".to_string(); 20]);
        name_insertions(&mut test_vec, &Scheme::IMGT);
        assert!(test_vec.contains(&String::from("10CV")));
    }
}
