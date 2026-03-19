import immunum
import pytest
import pickle


ALL_CHAINS = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD"]
AB_CHAINS = ["IGH", "IGK", "IGL"]
TCR_CHAINS = ["TRA", "TRB", "TRG", "TRD"]

IGL_SEQ = "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA"
IGH_SEQ = (
    "EVQLVESGGGLVKPGGSLKLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNAKN"
)
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

    def test_number_with_single_chain(self):
        annotator = immunum.Annotator(["IGH"], "IMGT")
        result = annotator.number(IGH_SEQ)
        assert result.chain == "H"

    def test_empty_sequence_raises(self):
        annotator = immunum.Annotator(ALL_CHAINS, "IMGT")
        with pytest.raises(ValueError):
            annotator.number("")

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
