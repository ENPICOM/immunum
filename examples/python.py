from immunum import Annotator, Chain, Scheme, default_scoring_params

my_sequence = "AAAAAAAAAAAA"


"""
The Annotator class is the starting point for numbering. It needs to be initialized with a scheme, chains and can take option parameters like custom scoring parameters.
Then the annotator will supply methods for numbering. We can number a single sequence, many in parallel or use files as input.
"""


# Numbering a sequence using a specific scheme and chain
annotator = Annotator(scheme=Scheme.IMGT, chains=Chain.IGH)
result = annotator.number_sequence(sequence=my_sequence)

# Numbering of sequence in bulk, checking multiple chains for each
my_sequences = [my_sequence] * 1000
annotator = Annotator(scheme=Scheme.IMGT, chains=[Chain.IGK, Chain.IGL])
result = annotator.number_sequences(sequences=my_sequences, parallel=True)

# Numbering a sequence with a custom penalties
scoring_params = default_scoring_params()
scoring_params.gap_pen_op = 1.5
annotator = Annotator(
    scheme=Scheme.IMGT, chains=Chain.IGH, scoring_params=scoring_params
)
result = annotator.number_sequence(sequence=my_sequence)

# Output structure
# result.numbers: list[str] # Number for each index position of the sequence
# result.sequence: str
# result.scheme: Scheme
# result.chain: Chain
# result.identity: int # Confidence alignment score
# result.regions: dict[str, (int, int)] # Start and end position of each region
