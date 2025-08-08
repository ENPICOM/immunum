"""
Basic tests for the immunum Python API.

This module contains tests for the new Annotator-based API functionality.
"""

import pytest
import immunum


class TestImmunumAPI:
    """Tests for the main immunum API."""

    def test_scheme_enums(self):
        """Test that scheme enums are accessible."""
        assert immunum.Scheme.IMGT is not None
        assert immunum.Scheme.KABAT is not None

    def test_chain_enums(self):
        """Test that chain enums are accessible."""
        assert immunum.Chain.IGH is not None
        assert immunum.Chain.IGK is not None
        assert immunum.Chain.IGL is not None
        assert immunum.Chain.TRA is not None
        assert immunum.Chain.TRB is not None
        assert immunum.Chain.TRG is not None
        assert immunum.Chain.TRD is not None

    def test_scoring_params(self):
        """Test ScoringParams creation."""
        params = immunum.ScoringParams()
        assert params is not None
        assert params.gap_pen_cp > 0

        custom_params = immunum.ScoringParams(gap_pen_cp=60.0)
        assert custom_params.gap_pen_cp == 60.0

    def test_default_scoring_params(self):
        """Test default scoring parameters function."""
        params = immunum.default_scoring_params()
        assert params is not None
        assert params.gap_pen_cp > 0

    def test_annotator_creation(self):
        """Test Annotator creation."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT, chains=[immunum.Chain.IGH]
        )
        assert annotator is not None

    def test_annotator_with_custom_params(self):
        """Test Annotator creation with custom scoring parameters."""
        custom_params = immunum.ScoringParams(gap_pen_cp=60.0)
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT,
            chains=[immunum.Chain.IGH],
            scoring_params=custom_params,
        )
        assert annotator is not None

    def test_annotator_with_prefiltering(self):
        """Test Annotator creation with pre-filtering enabled."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT,
            chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
            use_prefiltering=True,
        )
        assert annotator is not None

    def test_api_structure(self):
        """Test that the expected API is available."""
        expected_api = [
            "Annotator",
            "AnnotationResult",
            "Chain",
            "Scheme",
            "ScoringParams",
            "default_scoring_params",
        ]
        available_api = [attr for attr in dir(immunum) if not attr.startswith("_")]

        for api in expected_api:
            assert api in available_api, (
                f"Expected API '{api}' not found. Available: {available_api}"
            )


class TestAnnotationResult:
    """Test AnnotationResult functionality."""

    def test_annotation_result_properties(self):
        """Test that AnnotationResult has expected properties."""
        # We'll need to create a mock or use a real result
        # For now, just test the class exists
        assert hasattr(immunum, "AnnotationResult")

    def test_annotation_result_methods(self):
        """Test that AnnotationResult has expected methods."""
        # This would require creating an actual result to test
        # The methods exist based on our implementation:
        # - sequence, numbers, scheme, chain, identity (getters)
        # - regions, start, end (getters)
        # - get_region_sequence, get_cdr_sequences, get_framework_sequences
        # - is_high_confidence, summary
        pass


class TestAnnotatorMethods:
    """Test Annotator methods."""

    def test_annotator_number_sequence(self):
        """Test Annotator.number_sequence method."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT, chains=[immunum.Chain.IGH]
        )

        # Test with a short sequence - might fail but shouldn't crash
        test_sequence = "ATCGATCGATCG"
        results = annotator.number_sequence(test_sequence)
        assert isinstance(results, list)

    def test_annotator_number_sequences(self):
        """Test Annotator.number_sequences method."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT, chains=[immunum.Chain.IGH]
        )

        sequences = ["ATCGATCGATCG", "GCTAGCTAGCTA"]
        results = annotator.number_sequences(sequences)
        assert isinstance(results, list)
        assert len(results) == len(sequences)

    def test_annotator_paired_flag(self):
        """Test paired-mode via parameterized API."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT, chains=[immunum.Chain.IGH, immunum.Chain.IGK]
        )

        # Single sequence paired
        results = annotator.number_sequence("ATCGATCGATCG", all_chains=True)
        assert isinstance(results, list)

        # Multiple sequences paired
        batch = annotator.number_sequences(["ATCGATCGATCG", "GCTAGCTAGCTA"], all_chains=True)
        assert isinstance(batch, list)
        for item in batch:
            assert isinstance(item, list)

        # File paired
        results_file = annotator.number_file("fixtures/test.fasta", all_chains=True)
        assert isinstance(results_file, list)

    def test_annotator_number_file(self):
        """Test Annotator.number_file method."""
        annotator = immunum.Annotator(
            scheme=immunum.Scheme.IMGT,
            chains=[immunum.Chain.IGH, immunum.Chain.IGK, immunum.Chain.IGL],
        )

        # Test with the fixture file (sequential processing)
        results = annotator.number_file("fixtures/test.fasta")
        assert isinstance(results, list)
        assert len(results) == 4  # Expected 4 sequences in test.fasta

        # Check that results are tuples of (name, AnnotationResult)
        for name, result in results:
            assert isinstance(name, str)
            assert hasattr(result, "sequence")
            assert hasattr(result, "scheme")
            assert hasattr(result, "chain")
            assert hasattr(result, "identity")

        # Test with parallel processing
        parallel_results = annotator.number_file("fixtures/test.fasta", parallel=True)
        assert isinstance(parallel_results, list)
        assert len(parallel_results) == 4  # Same number of results

        # Test with non-existent file
        try:
            annotator.number_file("fixtures/nonexistent.fasta")
            assert False, "Should have raised an exception for non-existent file"
        except Exception as e:
            assert "not found" in str(e).lower() or "error" in str(e).lower()
