# Sequence and chain-set fixtures shared across test files.
# Mirrors the constants in tests/test_python.py.

ALL_CHAINS_R <- c("IGH", "IGK", "IGL", "TRA", "TRB", "TRG", "TRD")
AB_CHAINS_R  <- c("IGH", "IGK", "IGL")
TCR_CHAINS_R <- c("TRA", "TRB", "TRG", "TRD")

IGL_SEQ <- "SALTQPPAVSGTPGQRVTISCSGSDIGRRSVNWYQQFPGTAPKLLIYSNDQRPSVVPDRFSGSKSGTSASLAISGLQSEDEAEYYCAAWDDSLAVFGGGTQLTVGQPKA"
IGH_SEQ <- "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS"
TRA_SEQ <- "DSVTQTEGQVALSEEDFLTIHCNYSASGYPALFWYVQYPGEGPQFLFRASRDKEKGSSRGFEATYNKEATSFHLQKASVQESDSAVYYCALSGGNNKLTFGAGTKLTIKP"
TRB_SEQ <- "GVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTE"
TRG_SEQ <- "AGHLEQPQISSTKTLSKTARLECVVSGITISATSVYWYRERPGEVIQFLVSISYDGTVRKESGIPSGKFEVDRIPETSTSTLTIHNVEKQDIATYYCALWEAQQEGLKKIKVFGPGTKLIITD"
TRD_SEQ <- "QKVTQAQSSVSMPVRKAVTLNCLYETSWWSYYIFWYKQLPSKEMIFLIRQGSDEQNAKSGRYSVNFKKAAKSVALTISALQLEDSAKYFCALGDPGGNLTDKLIFGKGTRVTVEP"

INIT_CASES <- list(
  list(id = "single_IGH_imgt",       chains = "IGH",          scheme = "IMGT",  seq = IGH_SEQ),
  list(id = "single_IGK_imgt",       chains = "IGK",          scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "single_IGL_imgt",       chains = "IGL",          scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "single_TRA_imgt",       chains = "TRA",          scheme = "IMGT",  seq = TRA_SEQ),
  list(id = "single_TRB_imgt",       chains = "TRB",          scheme = "IMGT",  seq = TRB_SEQ),
  list(id = "single_TRG_imgt",       chains = "TRG",          scheme = "IMGT",  seq = TRG_SEQ),
  list(id = "single_TRD_imgt",       chains = "TRD",          scheme = "IMGT",  seq = TRD_SEQ),
  list(id = "all_ab_chains_imgt",    chains = AB_CHAINS_R,    scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "all_tcr_chains_imgt",   chains = TCR_CHAINS_R,   scheme = "IMGT",  seq = TRB_SEQ),
  list(id = "all_chains_imgt",       chains = ALL_CHAINS_R,   scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "single_IGH_kabat",      chains = "IGH",          scheme = "Kabat", seq = IGH_SEQ),
  list(id = "all_ab_chains_kabat",   chains = AB_CHAINS_R,    scheme = "Kabat", seq = IGL_SEQ),
  list(id = "lowercase_chain_scheme", chains = "igh",         scheme = "imgt",  seq = IGH_SEQ),
  list(id = "short_chain_alias_H",   chains = "H",            scheme = "IMGT",  seq = IGH_SEQ),
  list(id = "short_chain_alias_K",   chains = "K",            scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "short_chain_alias_L",   chains = "L",            scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "named_chain_heavy",     chains = "heavy",        scheme = "IMGT",  seq = IGH_SEQ),
  list(id = "named_chain_kappa",     chains = "kappa",        scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "named_chain_lambda",    chains = "lambda",       scheme = "IMGT",  seq = IGL_SEQ),
  list(id = "named_chain_alpha",     chains = "alpha",        scheme = "IMGT",  seq = TRA_SEQ),
  list(id = "named_chain_beta",      chains = "beta",         scheme = "IMGT",  seq = TRB_SEQ),
  list(id = "short_chain_alias_A",   chains = "A",            scheme = "IMGT",  seq = TRA_SEQ),
  list(id = "short_chain_alias_B",   chains = "B",            scheme = "IMGT",  seq = TRB_SEQ),
  list(id = "short_chain_alias_G",   chains = "G",            scheme = "IMGT",  seq = TRG_SEQ),
  list(id = "short_chain_alias_D",   chains = "D",            scheme = "IMGT",  seq = TRD_SEQ),
  list(id = "named_chain_gamma",     chains = "gamma",        scheme = "IMGT",  seq = TRG_SEQ),
  list(id = "named_chain_delta",     chains = "delta",        scheme = "IMGT",  seq = TRD_SEQ),
  list(id = "scheme_alias_i",        chains = "IGH",          scheme = "i",     seq = IGH_SEQ),
  list(id = "scheme_alias_k",        chains = AB_CHAINS_R,    scheme = "k",     seq = IGL_SEQ)
)
