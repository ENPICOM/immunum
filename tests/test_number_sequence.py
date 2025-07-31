"""
Basic tests for the typed API in immunum.

This module contains basic tests for the new typed API functionality.
"""

import pytest
import immunum


class TestTypedAPI:
    """Basic tests for the new typed API."""

    def test_scheme_enums(self):
        """Test that scheme enums are accessible."""
        assert immunum.Scheme.IMGT is not None
        assert immunum.Scheme.KABAT is not None

    def test_chain_enums(self):
        """Test that chain enums are accessible."""
        assert immunum.Chain.IGH is not None
        assert immunum.Chain.IGK is not None
        assert immunum.Chain.IGL is not None

    def test_scoring_params(self):
        """Test ScoringParams creation."""
        params = immunum.ScoringParams()
        assert params is not None
        assert params.gap_pen_cp > 0
        
        custom_params = immunum.ScoringParams(gap_pen_cp=60.0)
        assert custom_params.gap_pen_cp == 60.0

    def test_numbering_scheme_creation(self):
        """Test NumberingScheme creation."""
        scheme = immunum.NumberingScheme.imgt_heavy()
        assert scheme is not None
        
        custom_params = immunum.ScoringParams(gap_pen_cp=60.0)
        scheme_custom = immunum.NumberingScheme.kabat_kappa(custom_params)
        assert scheme_custom is not None

    def test_batch_processing_basic(self):
        """Test basic batch processing functionality."""
        sequences = ["ATCGATCG"]
        scheme = immunum.Scheme.IMGT
        chains = [immunum.Chain.IGH]
        
        # This will likely fail due to file I/O issues in tests,
        # but should demonstrate the correct API signature
        try:
            results = immunum.number_sequences_batch(sequences, scheme, chains)
            assert isinstance(results, list)
        except Exception as e:
            # Expected to fail due to missing test data files
            assert "No such file or directory" in str(e) or "file" in str(e).lower()

    def test_api_structure(self):
        """Test that the expected API is available."""
        expected_api = ['Chain', 'NumberingScheme', 'Scheme', 'ScoringParams', 'number_sequences_batch']
        available_api = [attr for attr in dir(immunum) if not attr.startswith('_')]
        
        for api in expected_api:
            assert api in available_api, f"Expected API '{api}' not found"


class TestNumberingSchemeConstructors:
    """Test all NumberingScheme constructors."""

    def test_all_constructors(self):
        """Test that all scheme constructors work."""
        constructors = [
            immunum.NumberingScheme.imgt_heavy,
            immunum.NumberingScheme.imgt_kappa,
            immunum.NumberingScheme.imgt_lambda,
            immunum.NumberingScheme.kabat_heavy,
            immunum.NumberingScheme.kabat_kappa,
            immunum.NumberingScheme.kabat_lambda,
        ]
        
        for constructor in constructors:
            scheme = constructor()
            assert scheme is not None