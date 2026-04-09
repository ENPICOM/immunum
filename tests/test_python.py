import immunum
import pytest
import pickle


ALL_CHAINS = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
AB_CHAINS = ["IGH", "IGK", "IGL"]
TCR_CHAINS = ["TRA", "TRB", "TRG", "TRD"]

IGL_SEQ = "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA"
IGH_SEQ = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
TRA_SEQ = "DSVTQTEGQVALSEEDFLTIHCNYSASGYPALFWYVQYPGEGPQFLFRASRDKEKGSSRGFEATYNKEATSFHLQKASVQESDSAVYYCALSGGNNKLTFGAGTKLTIKP"
TRB_SEQ = "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE"
TRG_SEQ = "AGHLEQPQISSTKTLSKTARLECVVSGITISATSVYWYRERPGEVIQFLVSISYDGTVRKESGIPSGKFEVDRIPETSTSTLTIHNVEKQDIATYYCALWEAQQEGLKKIKVFGPGTKLIITD"
TRD_SEQ = "QKVTQAQSSVSMPVRKAVTLNCLYETSWWSYYIFWYKQLPSKEMIFLIRQGSDEQNAKSGRYSVNFKKAAKSVALTISALQLEDSAKYFCALGDPGGNLTDKLIFGKGTRVTVEP"


@pytest.fixture(
    params=[
        # (chains, scheme, seq) — all valid initialization scenarios
        pytest.param((["IGH"], "IMGT", IGH_SEQ), id="single_IGH_imgt"),
        pytest.param((["IGK"], "IMGT", IGL_SEQ), id="single_IGK_imgt"),
        pytest.param((["IGL"], "IMGT", IGL_SEQ), id="single_IGL_imgt"),
        pytest.param((["TRA"], "IMGT", TRA_SEQ), id="single_TRA_imgt"),
        pytest.param((["TRB"], "IMGT", TRB_SEQ), id="single_TRB_imgt"),
        pytest.param((["TRG"], "IMGT", TRG_SEQ), id="single_TRG_imgt"),
        pytest.param((["TRD"], "IMGT", TRD_SEQ), id="single_TRD_imgt"),
        pytest.param((AB_CHAINS, "IMGT", IGL_SEQ), id="all_ab_chains_imgt"),
        pytest.param((TCR_CHAINS, "IMGT", TRB_SEQ), id="all_tcr_chains_imgt"),
        pytest.param((ALL_CHAINS, "IMGT", IGL_SEQ), id="all_chains_imgt"),
        pytest.param((["IGH"], "Kabat", IGH_SEQ), id="single_IGH_kabat"),
        pytest.param((AB_CHAINS, "Kabat", IGL_SEQ), id="all_ab_chains_kabat"),
        # Case-insensitive chain/scheme aliases
        pytest.param((["igh"], "imgt", IGH_SEQ), id="lowercase_chain_and_scheme"),
        pytest.param((["H"], "IMGT", IGH_SEQ), id="short_chain_alias_H"),
        pytest.param((["K"], "IMGT", IGL_SEQ), id="short_chain_alias_K"),
        pytest.param((["L"], "IMGT", IGL_SEQ), id="short_chain_alias_L"),
        pytest.param((["heavy"], "IMGT", IGH_SEQ), id="named_chain_alias_heavy"),
        pytest.param((["kappa"], "IMGT", IGL_SEQ), id="named_chain_alias_kappa"),
        pytest.param((["lambda"], "IMGT", IGL_SEQ), id="named_chain_alias_lambda"),
        pytest.param((["alpha"], "IMGT", TRA_SEQ), id="named_chain_alias_alpha"),
        pytest.param((["beta"], "IMGT", TRB_SEQ), id="named_chain_alias_beta"),
        pytest.param((["A"], "IMGT", TRA_SEQ), id="short_chain_alias_A"),
        pytest.param((["B"], "IMGT", TRB_SEQ), id="short_chain_alias_B"),
        pytest.param((["G"], "IMGT", TRG_SEQ), id="short_chain_alias_G"),
        pytest.param((["D"], "IMGT", TRD_SEQ), id="short_chain_alias_D"),
        pytest.param((["gamma"], "IMGT", TRG_SEQ), id="named_chain_alias_gamma"),
        pytest.param((["delta"], "IMGT", TRD_SEQ), id="named_chain_alias_delta"),
        pytest.param((["IGH"], "i", IGH_SEQ), id="scheme_alias_i"),
        pytest.param((AB_CHAINS, "k", IGL_SEQ), id="scheme_alias_k"),
    ]
)
def annotator_and_seq(request):
    chains, scheme, seq = request.param
    return immunum.Annotator(chains, scheme), seq


class TestAnnotatorInit:
    def test_valid_init(self, annotator_and_seq):
        annotator, _ = annotator_and_seq
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

    def test_number_smoke(self, annotator_and_seq):
        annotator, seq = annotator_and_seq
        annotator.number(seq)

    def test_pickle(self, annotator_and_seq):
        annotator, seq = annotator_and_seq
        re_annotator = pickle.loads(pickle.dumps(annotator))
        re_annotator.number(seq)


class TestNumbering:
    def test_number_igh_sequence(self):
        annotator = immunum.Annotator(ALL_CHAINS, "IMGT")
        result = annotator.number(IGH_SEQ)
        assert result.chain == "H"
        assert result.scheme == "IMGT"
        assert result.error is None

    def test_number_with_single_chain(self):
        annotator = immunum.Annotator(["IGH"], "IMGT")
        result = annotator.number(IGH_SEQ)
        assert result.chain == "H"

    def test_query_start_end_on_success(self):
        annotator = immunum.Annotator(["IGH"], "IMGT")
        result = annotator.number(IGH_SEQ)
        assert isinstance(result.query_start, int)
        assert isinstance(result.query_end, int)
        assert 0 <= result.query_start <= result.query_end < len(IGH_SEQ)
        assert result.query_end - result.query_start + 1 == len(result.numbering)

    def test_empty_sequence_returns_error(self):
        annotator = immunum.Annotator(ALL_CHAINS, "IMGT")
        result = annotator.number("")
        assert result.error is not None
        assert result.chain is None
        assert result.query_start is None
        assert result.query_end is None

    def test_invalid_sequence_returns_error(self):
        annotator = immunum.Annotator(ALL_CHAINS, "IMGT")
        result = annotator.number("AAAAAAAAAAAAAAAA")
        assert result.error is not None
        assert result.chain is None

    def test_confidence_is_float(self, annotator_and_seq):
        annotator, seq = annotator_and_seq
        result = annotator.number(seq)
        assert isinstance(result.confidence, float)

    def test_segmentation(self, annotator_and_seq):
        annotator, seq = annotator_and_seq
        result = annotator.segment(seq)
        assert set(result.as_dict()) == {f"cdr{i}" for i in (1, 2, 3)} | {
            f"fr{i}" for i in (1, 2, 3, 4)
        } | {"prefix", "postfix"}
        assert result.error is None

    def test_segmentation_invalid_sequence_returns_error(self):
        annotator = immunum.Annotator(ALL_CHAINS, "IMGT")
        result = annotator.segment("AAAAAAAAAAAAAAAA")
        assert result.error is not None
        assert result.fr1 is None


class TestNormalization:
    @pytest.mark.parametrize(
        "alias_chains,canonical_chains,scheme,seq",
        [
            (["H"], ["IGH"], "IMGT", IGH_SEQ),
            (["K"], ["IGK"], "IMGT", IGL_SEQ),
            (["L"], ["IGL"], "IMGT", IGL_SEQ),
            (["A"], ["TRA"], "IMGT", TRA_SEQ),
            (["B"], ["TRB"], "IMGT", TRB_SEQ),
            (["G"], ["TRG"], "IMGT", TRG_SEQ),
            (["D"], ["TRD"], "IMGT", TRD_SEQ),
            (["heavy"], ["IGH"], "IMGT", IGH_SEQ),
            (["kappa"], ["IGK"], "IMGT", IGL_SEQ),
            (["lambda"], ["IGL"], "IMGT", IGL_SEQ),
            (["alpha"], ["TRA"], "IMGT", TRA_SEQ),
            (["beta"], ["TRB"], "IMGT", TRB_SEQ),
            (["gamma"], ["TRG"], "IMGT", TRG_SEQ),
            (["delta"], ["TRD"], "IMGT", TRD_SEQ),
            (["igh"], ["IGH"], "imgt", IGH_SEQ),
            (["IGH"], ["IGH"], "i", IGH_SEQ),
            (AB_CHAINS, AB_CHAINS, "k", IGL_SEQ),
        ],
    )
    def test_alias_produces_identical_result(
        self, alias_chains, canonical_chains, scheme, seq
    ):
        alias_result = immunum.Annotator(alias_chains, scheme).number(seq)
        canonical_result = immunum.Annotator(
            canonical_chains, "IMGT" if scheme in ("i", "imgt", "IMGT") else "Kabat"
        ).number(seq)
        assert alias_result.chain == canonical_result.chain
        assert alias_result.scheme == canonical_result.scheme
        assert alias_result.numbering == canonical_result.numbering

    @pytest.mark.parametrize(
        "chains,scheme",
        [
            (["INVALID"], "IMGT"),
            (["IGH"], "INVALID"),
            (["Z"], "IMGT"),
            (["IGX"], "IMGT"),
        ],
    )
    def test_unknown_alias_raises(self, chains, scheme):
        with pytest.raises(ValueError):
            immunum.Annotator(chains, scheme)
