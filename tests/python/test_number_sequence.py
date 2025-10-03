"""
Basic tests for the immunum Python API.

This module contains tests for the new Annotator-based API functionality.
"""

import pytest
import immunum


class TestImmunumAPI:
    """Tests for the main immunum API."""

    def test_scheme_enum_class(self):
        """Test that Scheme enum class is accessible."""
        assert hasattr(immunum, "Scheme")
        assert immunum.Scheme is not None

    def test_chain_enum_class(self):
        """Test that Chain enum class is accessible."""
        assert hasattr(immunum, "Chain")
        assert immunum.Chain is not None

    def test_annotator_creation_default(self):
        """Test Annotator creation with default parameters."""
        annotator = immunum.Annotator()
        assert annotator is not None

    def test_annotator_creation_custom(self):
        """Test Annotator creation with custom parameters."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT,
            chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
            disable_prefiltering=False,
            threads=2,
            min_confidence=0.8,
        )
        assert annotator is not None

    def test_annotator_creation_kabat(self):
        """Test Annotator creation with KABAT scheme."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.KABAT,
            chains=[immunum.Chain.IGH],
            disable_prefiltering=True,
            min_confidence=0.5,
        )
        assert annotator is not None

    def test_api_structure(self):
        """Test that the expected API is available."""
        expected_api = [
            "Annotator",
            "Chain",
            "Scheme",
        ]
        available_api = [attr for attr in dir(immunum) if not attr.startswith("_")]

        for api in expected_api:
            assert api in available_api, (
                f"Expected API '{api}' not found. Available: {available_api}"
            )


class TestAnnotatorMethods:
    """Test Annotator methods with the new API."""

    def test_number_sequences_basic(self):
        """Test basic number_sequences functionality."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT,
            chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
        )

        # Test with realistic antibody sequences
        heavy_chain = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS"
        light_chain = "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNSKNYLAWYQDKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYSFGQGTKLEIK"

        sequences = [heavy_chain, light_chain]
        results = annotator.number_sequences(sequences, max_chains=2)

        assert isinstance(results, list)
        assert len(results) == 2  # Two input sequences

        # Check heavy chain result
        heavy_result = results[0]
        assert len(heavy_result) >= 1  # Should find at least one chain
        numbers, confidence, chain_type = heavy_result[0]
        print(dir(chain_type))
        assert isinstance(numbers, list)
        assert isinstance(confidence, float)
        assert confidence >= 0.0 and confidence <= 1.0
        assert chain_type == immunum.Chain.IGH  # Should be heavy chain

        # Check light chain result
        light_result = results[1]
        assert len(light_result) >= 1  # Should find at least one chain
        numbers, confidence, chain_type = light_result[0]
        assert isinstance(numbers, list)
        assert isinstance(confidence, float)
        assert confidence >= 0.0 and confidence <= 1.0
        assert chain_type in [
            immunum.Chain.IGK,
            immunum.Chain.IGL,
        ]  # Should be light chain

    def test_number_sequences_with_tuples(self):
        """Test number_sequences with tuple input (id, sequence)."""
        annotator = immunum.Annotator()

        sequences = [
            (
                "heavy_seq",
                "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS",
            ),
            (
                "light_seq",
                "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNSKNYLAWYQDKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYSFGQGTKLEIK",
            ),
        ]

        results = annotator.number_sequences(sequences)
        assert isinstance(results, list)
        assert len(results) == 2

    def test_number_sequences_mixed_input(self):
        """Test number_sequences with mixed string and tuple input."""
        annotator = immunum.Annotator(min_confidence=0.5)

        sequences = [
            "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS",
            (
                "light_seq",
                "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNSKNYLAWYQDKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYSFGQGTKLEIK",
            ),
        ]

        results = annotator.number_sequences(sequences, max_chains=1)
        assert isinstance(results, list)
        assert len(results) == 2

    def test_number_sequences_empty(self):
        """Test number_sequences with empty input."""
        annotator = immunum.Annotator()

        results = annotator.number_sequences([])
        assert isinstance(results, list)
        assert len(results) == 0

    def test_number_sequences_no_chains_found(self):
        """Test number_sequences with sequences that likely won't match."""
        annotator = immunum.Annotator(min_confidence=0.9)  # High confidence threshold

        # Very short, non-antibody sequence
        sequences = ["ATCGATCG"]
        results = annotator.number_sequences(sequences)

        assert isinstance(results, list)
        assert len(results) == 1
        assert len(results[0]) == 0  # No chains found

    def test_number_sequences_max_chains(self):
        """Test max_chains parameter."""
        annotator = immunum.Annotator()

        # Test with a longer sequence that might contain multiple chains
        long_sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSSDIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNSKNYLAWYQDKPGQPPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYSTPYSFGQGTKLEIK"

        results_max1 = annotator.number_sequences([long_sequence], max_chains=1)
        results_max2 = annotator.number_sequences([long_sequence], max_chains=2)

        assert len(results_max1[0]) <= 1
        assert len(results_max2[0]) <= 2

    def test_confidence_filtering(self):
        """Test that min_confidence parameter works."""
        # Create annotators with different confidence thresholds
        low_confidence = immunum.Annotator(min_confidence=0.1)
        high_confidence = immunum.Annotator(min_confidence=0.95)

        sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS"

        low_results = low_confidence.number_sequences([sequence])
        high_results = high_confidence.number_sequences([sequence])

        # Low confidence should find chains (if any exist)
        # High confidence might filter some out
        assert len(low_results[0]) >= len(high_results[0])

    def test_different_schemes(self):
        """Test different numbering schemes."""
        imgt_annotator = immunum.Annotator(scheme=immunum.Scheme.IMGT)
        kabat_annotator = immunum.Annotator(scheme=immunum.Scheme.KABAT)

        sequence = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS"

        imgt_results = imgt_annotator.number_sequences([sequence])
        kabat_results = kabat_annotator.number_sequences([sequence])

        # Both should process without errors
        assert isinstance(imgt_results, list)
        assert isinstance(kabat_results, list)

        # If chains are found, the numbering might be different
        if len(imgt_results[0]) > 0 and len(kabat_results[0]) > 0:
            imgt_numbers = imgt_results[0][0][0]  # First chain's numbers
            kabat_numbers = kabat_results[0][0][0]  # First chain's numbers
            # Numbers might be different between schemes
            assert isinstance(imgt_numbers, list)
            assert isinstance(kabat_numbers, list)


class TestStringBasedAPI:
    """Tests for string-based scheme and chain arguments."""

    def test_annotator_with_string_scheme(self):
        """Test Annotator creation with string scheme."""
        annotator = immunum.Annotator(scheme="IMGT")
        assert annotator is not None

        annotator = immunum.Annotator(scheme="KABAT")
        assert annotator is not None

    def test_annotator_with_scheme_alias(self):
        """Test Annotator creation with scheme aliases."""
        annotator = immunum.Annotator(scheme="I")
        assert annotator is not None

        annotator = immunum.Annotator(scheme="K")
        assert annotator is not None

    def test_annotator_with_string_chains(self):
        """Test Annotator creation with string chains."""
        annotator = immunum.Annotator(chains=["IGH", "IGK", "IGL"])
        assert annotator is not None

    def test_annotator_with_chain_aliases(self):
        """Test Annotator creation with chain aliases."""
        annotator = immunum.Annotator(chains=["H", "K", "L"])
        assert annotator is not None

        annotator = immunum.Annotator(chains=["Heavy", "Kappa", "Lambda"])
        assert annotator is not None

    def test_annotator_string_and_enum_mixed(self):
        """Test that string and enum arguments work together."""
        annotator = immunum.Annotator(
            scheme="IMGT", chains=[immunum.Chain.IGH, immunum.Chain.IGK]
        )
        assert annotator is not None

    def test_number_sequences_with_string_scheme(self):
        """Test number_sequences with string scheme."""
        annotator = immunum.Annotator(scheme="IMGT", chains=["IGH", "IGK", "IGL"])

        heavy_chain = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS"
        sequences = [heavy_chain]
        results = annotator.number_sequences(sequences)

        assert isinstance(results, list)
        assert len(results) == 1
        assert len(results[0]) >= 1

    def test_number_sequences_with_kabat_string(self):
        """Test number_sequences with KABAT scheme as string."""
        annotator = immunum.Annotator(scheme="KABAT", chains=["IGH"])

        heavy_chain = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS"
        results = annotator.number_sequences([heavy_chain])

        assert isinstance(results, list)
        assert len(results) == 1

    def test_invalid_scheme_string(self):
        """Test that invalid scheme strings raise ValueError."""
        with pytest.raises(ValueError) as excinfo:
            immunum.Annotator(scheme="INVALID")
        assert "Invalid scheme" in str(excinfo.value)

    def test_invalid_chain_string(self):
        """Test that invalid chain strings raise ValueError."""
        with pytest.raises(ValueError) as excinfo:
            immunum.Annotator(chains=["INVALID"])
        assert "Invalid chain" in str(excinfo.value)

    def test_tcr_chains_with_strings(self):
        """Test T-cell receptor chains with string arguments."""
        annotator = immunum.Annotator(
            scheme="IMGT", chains=["TRA", "TRB", "TRG", "TRD"]
        )
        assert annotator is not None

    def test_tcr_chain_aliases(self):
        """Test T-cell receptor chain aliases."""
        annotator = immunum.Annotator(scheme="IMGT", chains=["A", "B", "G", "D"])
        assert annotator is not None

        annotator = immunum.Annotator(
            scheme="IMGT", chains=["Alpha", "Beta", "Gamma", "Delta"]
        )
        assert annotator is not None

    def test_backward_compatibility_enums(self):
        """Test that old enum-based API still works."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT,
            chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
        )
        assert annotator is not None

        heavy_chain = "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARGQHTLQGVVVVPATDYWGQGTLVTVSS"
        results = annotator.number_sequences([heavy_chain])

        assert isinstance(results, list)
        assert len(results) == 1
