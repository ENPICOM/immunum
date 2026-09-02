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
        pytest.param((["IGH"], "Chothia", IGH_SEQ), id="single_IGH_chothia"),
        pytest.param((AB_CHAINS, "Chothia", IGL_SEQ), id="all_ab_chains_chothia"),
        pytest.param((["IGH"], "Martin", IGH_SEQ), id="single_IGH_martin"),
        pytest.param((AB_CHAINS, "Martin", IGL_SEQ), id="all_ab_chains_martin"),
        pytest.param((["IGH"], "Aho", IGH_SEQ), id="single_IGH_aho"),
        pytest.param((AB_CHAINS, "Aho", IGL_SEQ), id="all_ab_chains_aho"),
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
        pytest.param((AB_CHAINS, "c", IGL_SEQ), id="scheme_alias_c"),
        pytest.param((AB_CHAINS, "m", IGL_SEQ), id="scheme_alias_m"),
        pytest.param((AB_CHAINS, "a", IGL_SEQ), id="scheme_alias_a"),
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
            (["TRA"], "Chothia"),
            (ALL_CHAINS, "Chothia"),
            (["TRA"], "Martin"),
            (ALL_CHAINS, "Martin"),
            (["TRA"], "Aho"),
            (ALL_CHAINS, "Aho"),
        ],
    )
    def test_non_imgt_tcr_raises(self, chains, scheme):
        with pytest.raises(
            ValueError, match=f"{scheme} scheme only supported for antibody chains"
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
        "alias_chains,canonical_chains,scheme,canonical_scheme,seq",
        [
            (["H"], ["IGH"], "IMGT", "IMGT", IGH_SEQ),
            (["K"], ["IGK"], "IMGT", "IMGT", IGL_SEQ),
            (["L"], ["IGL"], "IMGT", "IMGT", IGL_SEQ),
            (["A"], ["TRA"], "IMGT", "IMGT", TRA_SEQ),
            (["B"], ["TRB"], "IMGT", "IMGT", TRB_SEQ),
            (["G"], ["TRG"], "IMGT", "IMGT", TRG_SEQ),
            (["D"], ["TRD"], "IMGT", "IMGT", TRD_SEQ),
            (["heavy"], ["IGH"], "IMGT", "IMGT", IGH_SEQ),
            (["kappa"], ["IGK"], "IMGT", "IMGT", IGL_SEQ),
            (["lambda"], ["IGL"], "IMGT", "IMGT", IGL_SEQ),
            (["alpha"], ["TRA"], "IMGT", "IMGT", TRA_SEQ),
            (["beta"], ["TRB"], "IMGT", "IMGT", TRB_SEQ),
            (["gamma"], ["TRG"], "IMGT", "IMGT", TRG_SEQ),
            (["delta"], ["TRD"], "IMGT", "IMGT", TRD_SEQ),
            (["igh"], ["IGH"], "imgt", "IMGT", IGH_SEQ),
            (["IGH"], ["IGH"], "i", "IMGT", IGH_SEQ),
            (AB_CHAINS, AB_CHAINS, "k", "Kabat", IGL_SEQ),
            (AB_CHAINS, AB_CHAINS, "c", "Chothia", IGL_SEQ),
            (AB_CHAINS, AB_CHAINS, "m", "Martin", IGL_SEQ),
            (AB_CHAINS, AB_CHAINS, "a", "Aho", IGL_SEQ),
            (["igh"], ["IGH"], "chothia", "Chothia", IGH_SEQ),
            (["igh"], ["IGH"], "martin", "Martin", IGH_SEQ),
            (["igh"], ["IGH"], "aho", "Aho", IGH_SEQ),
        ],
    )
    def test_alias_produces_identical_result(
        self, alias_chains, canonical_chains, scheme, canonical_scheme, seq
    ):
        alias_result = immunum.Annotator(alias_chains, scheme).number(seq)
        canonical_result = immunum.Annotator(canonical_chains, canonical_scheme).number(
            seq
        )
        assert alias_result.chain == canonical_result.chain
        assert alias_result.scheme == canonical_scheme
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


# Region boundaries per scheme and chain, inclusive on both ends. Pinned literals on purpose:
# these are the numbers downstream code diffs its own FR/CDR definitions against, so they have to
# fail loudly when a table moves rather than track it.
IMGT_REGIONS = {
    "fr1": (1, 26),
    "cdr1": (27, 38),
    "fr2": (39, 55),
    "cdr2": (56, 65),
    "fr3": (66, 104),
    "cdr3": (105, 117),
    "fr4": (118, 128),
}
AHO_REGIONS = {
    "fr1": (1, 24),
    "cdr1": (25, 42),
    "fr2": (43, 56),
    "cdr2": (57, 77),
    "fr3": (78, 108),
    "cdr3": (109, 138),
    "fr4": (139, 149),
}
KABAT_HEAVY_REGIONS = {
    "fr1": (1, 30),
    "cdr1": (31, 35),
    "fr2": (36, 49),
    "cdr2": (50, 65),
    "fr3": (66, 94),
    "cdr3": (95, 102),
    "fr4": (103, 113),
}
KABAT_LIGHT_REGIONS = {
    "fr1": (1, 23),
    "cdr1": (24, 34),
    "fr2": (35, 49),
    "cdr2": (50, 56),
    "fr3": (57, 88),
    "cdr3": (89, 97),
    "fr4": (98, 107),
}
CHOTHIA_HEAVY_REGIONS = {
    "fr1": (1, 25),
    "cdr1": (26, 32),
    "fr2": (33, 51),
    "cdr2": (52, 56),
    "fr3": (57, 95),
    "cdr3": (96, 101),
    "fr4": (102, 113),
}
CHOTHIA_LIGHT_REGIONS = {
    "fr1": (1, 25),
    "cdr1": (26, 32),
    "fr2": (33, 49),
    "cdr2": (50, 52),
    "fr3": (53, 90),
    "cdr3": (91, 96),
    "fr4": (97, 107),
}
MARTIN_HEAVY_REGIONS = {
    "fr1": (1, 25),
    "cdr1": (26, 35),
    "fr2": (36, 49),
    "cdr2": (50, 58),
    "fr3": (59, 94),
    "cdr3": (95, 102),
    "fr4": (103, 113),
}
MARTIN_LIGHT_REGIONS = {
    "fr1": (1, 23),
    "cdr1": (24, 34),
    "fr2": (35, 49),
    "cdr2": (50, 56),
    "fr3": (57, 88),
    "cdr3": (89, 97),
    "fr4": (98, 107),
}

REGION_ORDER = ["fr1", "cdr1", "fr2", "cdr2", "fr3", "cdr3", "fr4"]

AB_SCHEMES = ["IMGT", "Kabat", "Chothia", "Martin", "Aho"]
NON_IMGT_SCHEMES = ["Kabat", "Chothia", "Martin", "Aho"]


class TestRegionsFor:
    @pytest.mark.parametrize(
        "scheme,chain,expected",
        [
            ("IMGT", "IGH", IMGT_REGIONS),
            ("IMGT", "IGK", IMGT_REGIONS),
            ("IMGT", "IGL", IMGT_REGIONS),
            ("Aho", "IGH", AHO_REGIONS),
            ("Aho", "IGK", AHO_REGIONS),
            ("Aho", "IGL", AHO_REGIONS),
            ("Kabat", "IGH", KABAT_HEAVY_REGIONS),
            ("Kabat", "IGK", KABAT_LIGHT_REGIONS),
            ("Kabat", "IGL", KABAT_LIGHT_REGIONS),
            ("Chothia", "IGH", CHOTHIA_HEAVY_REGIONS),
            ("Chothia", "IGK", CHOTHIA_LIGHT_REGIONS),
            ("Chothia", "IGL", CHOTHIA_LIGHT_REGIONS),
            ("Martin", "IGH", MARTIN_HEAVY_REGIONS),
            ("Martin", "IGK", MARTIN_LIGHT_REGIONS),
            ("Martin", "IGL", MARTIN_LIGHT_REGIONS),
        ],
    )
    def test_boundaries_match_pinned_table(self, scheme, chain, expected):
        assert immunum.regions_for(scheme, chain) == expected

    @pytest.mark.parametrize("scheme", AB_SCHEMES)
    @pytest.mark.parametrize("chain", AB_CHAINS)
    def test_regions_are_contiguous_from_position_one(self, scheme, chain):
        regions = immunum.regions_for(scheme, chain)
        assert list(regions) == REGION_ORDER
        starts = [start for start, _ in regions.values()]
        ends = [end for _, end in regions.values()]
        assert starts[0] == 1
        assert starts[1:] == [end + 1 for end in ends[:-1]]

    @pytest.mark.parametrize("scheme", ["IMGT", "Aho"])
    def test_chain_agnostic_schemes_share_one_table(self, scheme):
        heavy = immunum.regions_for(scheme, "IGH")
        assert immunum.regions_for(scheme, "IGK") == heavy
        assert immunum.regions_for(scheme, "IGL") == heavy

    @pytest.mark.parametrize("scheme", ["Kabat", "Chothia", "Martin"])
    def test_chain_specific_schemes_differ_between_heavy_and_light(self, scheme):
        assert immunum.regions_for(scheme, "IGH") != immunum.regions_for(scheme, "IGK")

    @pytest.mark.parametrize("scheme", NON_IMGT_SCHEMES)
    @pytest.mark.parametrize("chain", TCR_CHAINS)
    def test_non_imgt_tcr_raises(self, scheme, chain):
        with pytest.raises(
            ValueError, match=f"{scheme} scheme only supported for antibody chains"
        ):
            immunum.regions_for(scheme, chain)

    @pytest.mark.parametrize("chain", TCR_CHAINS)
    def test_imgt_covers_tcr_chains(self, chain):
        assert immunum.regions_for("IMGT", chain) == IMGT_REGIONS

    @pytest.mark.parametrize(
        "scheme,chain",
        [
            ("INVALID", "IGH"),
            ("Z", "IGH"),
            ("", "IGH"),
            ("IMGT", "INVALID"),
            ("IMGT", "IGX"),
            ("IMGT", ""),
        ],
    )
    def test_unknown_scheme_or_chain_raises(self, scheme, chain):
        with pytest.raises(ValueError):
            immunum.regions_for(scheme, chain)

    @pytest.mark.parametrize(
        "alias_scheme,alias_chain,scheme,chain",
        [
            ("k", "h", "Kabat", "IGH"),
            ("kabat", "heavy", "Kabat", "IGH"),
            ("c", "kappa", "Chothia", "IGK"),
            ("m", "l", "Martin", "IGL"),
            ("a", "lambda", "Aho", "IGL"),
            ("i", "b", "IMGT", "TRB"),
        ],
    )
    def test_aliases_resolve_to_the_same_table(
        self, alias_scheme, alias_chain, scheme, chain
    ):
        assert immunum.regions_for(alias_scheme, alias_chain) == immunum.regions_for(
            scheme, chain
        )
