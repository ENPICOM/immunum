import immunum
import pytest
import pickle


ALL_CHAINS = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
AB_CHAINS = ["IGH", "IGK", "IGL"]
TCR_CHAINS = ["TRA", "TRB", "TRG", "TRD"]
SEQ = "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA"
IGH_SEQ = (
    "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN"
)


@pytest.fixture(
    params=[
        # (chains, scheme) — all valid initialization scenarios
        pytest.param((["IGH"], "IMGT"), id="single_IGH_imgt"),
        pytest.param((["IGK"], "IMGT"), id="single_IGK_imgt"),
        pytest.param((["IGL"], "IMGT"), id="single_IGL_imgt"),
        pytest.param((["TRA"], "IMGT"), id="single_TRA_imgt"),
        pytest.param((["TRB"], "IMGT"), id="single_TRB_imgt"),
        pytest.param((["TRG"], "IMGT"), id="single_TRG_imgt"),
        pytest.param((["TRD"], "IMGT"), id="single_TRD_imgt"),
        pytest.param((AB_CHAINS, "IMGT"), id="all_ab_chains_imgt"),
        pytest.param((TCR_CHAINS, "IMGT"), id="all_tcr_chains_imgt"),
        pytest.param((ALL_CHAINS, "IMGT"), id="all_chains_imgt"),
        pytest.param((["IGH"], "Kabat"), id="single_IGH_kabat"),
        pytest.param((AB_CHAINS, "Kabat"), id="all_ab_chains_kabat"),
        # Case-insensitive chain/scheme aliases
        pytest.param((["igh"], "imgt"), id="lowercase_chain_and_scheme"),
        pytest.param((["H"], "IMGT"), id="short_chain_alias_H"),
        pytest.param((["K"], "IMGT"), id="short_chain_alias_K"),
        pytest.param((["L"], "IMGT"), id="short_chain_alias_L"),
        pytest.param((["heavy"], "IMGT"), id="named_chain_alias_heavy"),
        pytest.param((["kappa"], "IMGT"), id="named_chain_alias_kappa"),
        pytest.param((["lambda"], "IMGT"), id="named_chain_alias_lambda"),
        pytest.param((["alpha"], "IMGT"), id="named_chain_alias_alpha"),
        pytest.param((["beta"], "IMGT"), id="named_chain_alias_beta"),
    ]
)
def annotator(request):
    chains, scheme = request.param
    return immunum.Annotator(chains, scheme)


class TestAnnotatorInit:
    def test_valid_init(self, annotator):
        assert annotator is not None

    @pytest.mark.parametrize(
        "chains,scheme",
        [
            (["TRA"], "Kabat"),
            (["TRB"], "Kabat"),
            (["IGH", "TRA"], "Kabat"),
            (ALL_CHAINS, "Kabat"),
        ],
    )
    def test_kabat_tcr_raises(self, chains, scheme):
        with pytest.raises(
            ValueError, match="Kabat scheme only supported for antibody chains"
        ):
            immunum.Annotator(chains, scheme)

    @pytest.mark.parametrize(
        "chains,scheme",
        [
            (["INVALID"], "IMGT"),
            (["IGH"], "INVALID"),
            ([], "IMGT"),
        ],
    )
    def test_invalid_args_raise(self, chains, scheme):
        with pytest.raises(ValueError):
            immunum.Annotator(chains, scheme)

    def test_number_smoke(self, annotator, seq=SEQ):
        annotator.number(seq)

    def test_pickle(self, annotator, seq=SEQ):
        re_annotator = pickle.loads(pickle.dumps(annotator))
        re_annotator.number(seq)


class TestNumbering:
    def test_number_igh_sequence(self):
        annotator = immunum.Annotator(ALL_CHAINS, "IMGT")
        result = annotator.number(IGH_SEQ)
        assert result.chain == "H"
        assert result.scheme == "IMGT"

    def test_number_with_single_chain(self):
        annotator = immunum.Annotator(["IGH"], "IMGT")
        result = annotator.number(IGH_SEQ)
        assert result.chain == "H"

    def test_empty_sequence_raises(self):
        annotator = immunum.Annotator(ALL_CHAINS, "IMGT")
        with pytest.raises(ValueError):
            annotator.number("")

    def test_confidence_is_float(self, annotator):
        result = annotator.number(SEQ)
        assert isinstance(result.confidence, float)

    def test_segmentation(self, annotator):
        result = annotator.segment(SEQ)
        assert set(result.as_dict()) == {f"cdr{i}" for i in (1, 2, 3)} | {
            f"fr{i}" for i in (1, 2, 3, 4)
        } | {"prefix", "postfix"}
