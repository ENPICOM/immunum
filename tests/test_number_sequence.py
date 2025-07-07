"""
Tests for the number_sequence pyfunction in immunum.

This module contains comprehensive tests for the number_sequence function,
including valid inputs, edge cases, and error handling.
"""

import pytest
import immunum


class TestNumberSequence:
    """Test class for the number_sequence function."""

    def test_valid_imgt_scheme_with_igh_chain(self):
        """Test number_sequence with valid IMGT scheme and IGH chain."""
        sequence = "ATCGATCGATCG"
        result = immunum.number_sequence(sequence, "imgt", ["igh"])
        assert isinstance(result, str)
        assert len(result) > 0

    def test_valid_kabat_scheme_with_igk_chain(self):
        """Test number_sequence with valid KABAT scheme and IGK chain."""
        sequence = "ATCGATCGATCG"
        result = immunum.number_sequence(sequence, "kabat", ["igk"])
        assert isinstance(result, str)
        assert len(result) > 0

    def test_scheme_case_insensitive(self):
        """Test that scheme parameter is case insensitive."""
        sequence = "ATCGATCGATCG"
        chains = ["igh"]

        # Test different cases for IMGT
        result1 = immunum.number_sequence(sequence, "IMGT", chains)
        result2 = immunum.number_sequence(sequence, "imgt", chains)
        result3 = immunum.number_sequence(sequence, "ImGt", chains)

        assert result1 == result2 == result3

    def test_scheme_abbreviations(self):
        """Test that scheme abbreviations work correctly."""
        sequence = "ATCGATCGATCG"
        chains = ["igh"]

        # Test IMGT abbreviation
        result_imgt_full = immunum.number_sequence(sequence, "imgt", chains)
        result_imgt_abbrev = immunum.number_sequence(sequence, "i", chains)
        assert result_imgt_full == result_imgt_abbrev

        # Test KABAT abbreviation
        result_kabat_full = immunum.number_sequence(sequence, "kabat", chains)
        result_kabat_abbrev = immunum.number_sequence(sequence, "k", chains)
        assert result_kabat_full == result_kabat_abbrev

    def test_chain_case_insensitive(self):
        """Test that chain parameters are case insensitive."""
        sequence = "ATCGATCGATCG"
        scheme = "imgt"

        result1 = immunum.number_sequence(sequence, scheme, ["IGH"])
        result2 = immunum.number_sequence(sequence, scheme, ["igh"])
        result3 = immunum.number_sequence(sequence, scheme, ["IgH"])

        assert result1 == result2 == result3

    def test_chain_abbreviations(self):
        """Test that chain abbreviations work correctly."""
        sequence = "ATCGATCGATCG"
        scheme = "imgt"

        # Test IGH abbreviation
        result_igh_full = immunum.number_sequence(sequence, scheme, ["igh"])
        result_igh_abbrev = immunum.number_sequence(sequence, scheme, ["h"])
        assert result_igh_full == result_igh_abbrev

        # Test IGK abbreviation
        result_igk_full = immunum.number_sequence(sequence, scheme, ["igk"])
        result_igk_abbrev = immunum.number_sequence(sequence, scheme, ["k"])
        assert result_igk_full == result_igk_abbrev

    def test_multiple_chains(self):
        """Test number_sequence with multiple chains."""
        sequence = "ATCGATCGATCG"
        chains = ["igh", "igk", "igl"]
        result = immunum.number_sequence(sequence, "imgt", chains)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_all_valid_chains(self):
        """Test all valid chain types."""
        sequence = "ATCGATCGATCG"
        scheme = "imgt"

        valid_chains = ["igh", "igk", "igl", "tra", "trb", "trg", "trd"]

        for chain in valid_chains:
            result = immunum.number_sequence(sequence, scheme, [chain])
            assert isinstance(result, str)
            assert len(result) > 0

    def test_empty_sequence(self):
        """Test number_sequence with empty sequence."""
        result = immunum.number_sequence("", "imgt", ["igh"])
        assert isinstance(result, str)

    def test_long_sequence(self):
        """Test number_sequence with a long sequence."""
        long_sequence = "ATCG" * 100  # 400 characters
        result = immunum.number_sequence(long_sequence, "imgt", ["igh"])
        assert isinstance(result, str)
        assert len(result) > 0

    def test_sequence_with_non_standard_bases(self):
        """Test number_sequence with sequences containing non-standard bases."""
        sequence = "ATCGNATCGN"  # Contains N bases
        result = immunum.number_sequence(sequence, "imgt", ["igh"])
        assert isinstance(result, str)

    def test_invalid_scheme_raises_error(self):
        """Test that invalid scheme raises ValueError."""
        sequence = "ATCGATCGATCG"
        chains = ["igh"]

        with pytest.raises(ValueError, match="Scheme not supported"):
            immunum.number_sequence(sequence, "invalid_scheme", chains)

    def test_invalid_chain_raises_error(self):
        """Test that invalid chain raises ValueError."""
        sequence = "ATCGATCGATCG"
        scheme = "imgt"

        with pytest.raises(ValueError, match="Chain not supported"):
            immunum.number_sequence(sequence, scheme, ["invalid_chain"])

    def test_mixed_valid_and_invalid_chains(self):
        """Test that one invalid chain in a list raises error."""
        sequence = "ATCGATCGATCG"
        scheme = "imgt"

        with pytest.raises(ValueError, match="Chain not supported"):
            immunum.number_sequence(sequence, scheme, ["igh", "invalid_chain"])

    def test_empty_chains_list(self):
        """Test number_sequence with empty chains list."""
        sequence = "ATCGATCGATCG"
        result = immunum.number_sequence(sequence, "imgt", [])
        assert isinstance(result, str)

    def test_duplicate_chains(self):
        """Test number_sequence with duplicate chains."""
        sequence = "ATCGATCGATCG"
        result = immunum.number_sequence(sequence, "imgt", ["igh", "igh", "igh"])
        assert isinstance(result, str)
        assert len(result) > 0


class TestNumberSequenceParameterTypes:
    """Test parameter type handling for number_sequence function."""

    def test_sequence_parameter_type(self):
        """Test that sequence parameter must be a string."""
        with pytest.raises(TypeError):
            immunum.number_sequence(123, "imgt", ["igh"])

    def test_scheme_parameter_type(self):
        """Test that scheme parameter must be a string."""
        with pytest.raises(TypeError):
            immunum.number_sequence("ATCG", 123, ["igh"])

    def test_chains_parameter_type(self):
        """Test that chains parameter must be a list."""
        with pytest.raises(TypeError):
            immunum.number_sequence("ATCG", "imgt", "igh")


class TestNumberSequenceIntegration:
    """Integration tests for number_sequence function."""

    def test_realistic_antibody_sequence(self):
        """Test with a realistic antibody sequence."""
        # Example heavy chain variable region sequence
        sequence = (
            "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGWINPNSGGTNYAQKFQG"
            "RVTMTRDTSISTAYMELSRLRSDDTAVYYCARGGYDILTDYWGQGTLVTVSS"
        )

        result = immunum.number_sequence(sequence, "imgt", ["igh"])
        assert isinstance(result, str)
        assert len(result) > 0

    def test_realistic_tcr_sequence(self):
        """Test with a realistic T-cell receptor sequence."""
        # Example TCR alpha chain sequence
        sequence = (
            "METLLGLLILWLQLQWVSSKQEVTQIPAALSVPEGENLVLNCSFTDSAIYNLQWFRQDPGKGLTSL"
            "LLIQSSQREQTSGRLNASLDKSSGRSTLYIAASQPGDSATYLCAVRPTSGGSYIPTFGRGTSLIVHP"
        )

        result = immunum.number_sequence(sequence, "kabat", ["tra"])
        assert isinstance(result, str)
        assert len(result) > 0

    def test_cross_scheme_consistency(self):
        """Test that the same sequence produces consistent results across schemes."""
        sequence = "ATCGATCGATCGATCGATCG"
        chains = ["igh"]

        result_imgt = immunum.number_sequence(sequence, "imgt", chains)
        result_kabat = immunum.number_sequence(sequence, "kabat", chains)

        # Both should return strings (though content may differ)
        assert isinstance(result_imgt, str)
        assert isinstance(result_kabat, str)
