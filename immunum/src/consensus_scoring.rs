use crate::constants::BLOSUM62;
use std::collections::HashMap;
use std::io;
use std::fs;

fn read_consensus_file(path: &str) -> Result<HashMap<u32, Vec<char>>, io::Error> {
    // TODO BUFREADER instead???
    // TODO convert to char??
    let content = fs::read_to_string(path);
    println!("{:?}", content);
    let content = content.unwrap_or("".to_string());
    let mut consensus_aas: HashMap<u32, Vec<char>> = HashMap::new();
    // Loop over every line of content
    let total_lines = content.lines().count();
    // Skip first and last line
    for line in content.lines().skip(1).take(total_lines - 2) {
        let split_line: Vec<&str> = line.split(',').collect();
        consensus_aas.insert(
            split_line[0].parse::<u32>().unwrap(),
            split_line[1..].iter()
                .flat_map(|s| s.chars())
                .collect()
        );
    }
    Ok(consensus_aas)
}

fn best_score_consensus(position:u32, residue: char, consensus:&HashMap<u32, Vec<char>>) -> i32{
    // Initialize to minimal score
    let to_check_residues: &Vec<char> =
        consensus.get(&position).expect("Position outside of consensus");

    // when any is allowed, best score is perfect match
    if to_check_residues.iter().all(|&c| c == '-'){
        return blosum_lookup(&residue, &residue)
    }
    // iterate over residues to get max score
    blosum_lookup(&residue,
                  to_check_residues
                      .iter()
                      .max_by_key(|&c| blosum_lookup(&residue, c))
                      .unwrap())
}

fn blosum_lookup(residue1:&char, residue2:&char) -> i32 {
    if residue1 < residue2 {
        let lookup: String = format!("{}{}", residue1, residue2);
        *BLOSUM62.get(&lookup).unwrap()
    }else{
        let lookup: String = format!("{}{}", residue2, residue1);
        *BLOSUM62.get(&lookup).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    //TODO tests for this part

    #[test]
    fn test_lookup(){
        // working
        assert_eq!(blosum_lookup(&'B', &'N'), 3);
        assert_eq!(blosum_lookup(&'W', &'Y'), 2);
        assert_eq!(blosum_lookup(&'X', &'E'), -1);
        // same either way of inputting arguments
        assert_eq!(blosum_lookup(&'C', &'S'),
                   blosum_lookup(&'S', &'C'));
    }

    #[test]
    fn finding_best_scores(){
        let consensus = read_consensus_file(
            r"C:\Anti_Num\numbering\consensus\IMGT_CONSENSUS_H.txt");
        assert!(consensus.is_ok());
        let consensus = consensus.unwrap();

        // correct scores
        assert_eq!(best_score_consensus(1, 'A', &consensus), -1);
        assert_eq!(best_score_consensus(1, 'B', &consensus), 4);
        assert_eq!(best_score_consensus(1, 'C', &consensus), -3);

        assert_eq!(best_score_consensus(128, 'X', &consensus), 0);
        assert_eq!(best_score_consensus(128, 'Y', &consensus), -2);
        assert_eq!(best_score_consensus(128, 'Z', &consensus), 0);

    }
    #[test]
    #[should_panic(expected = "Position outside of consensus")]
    fn look_outside_consensus(){
        let consensus = read_consensus_file(
            r"C:\Anti_Num\numbering\consensus\IMGT_CONSENSUS_H.txt");
        assert!(consensus.is_ok());
        let consensus = consensus.unwrap();
        // test wrong positions
        best_score_consensus(200, 'A', &consensus);
    }
}