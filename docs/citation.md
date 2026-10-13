# Citation

If you use immunum in your research, please cite it, with the version you used:
`immunum --version`, `immunum.__version__` in Python or `version()` in JavaScript prints it.
Numbering can change between versions (see the [changelog](changelog.md)), so the version makes
your results reproducible.

The citation metadata is in
[CITATION.cff](https://github.com/ENPICOM/immunum/blob/main/CITATION.cff), which GitHub's "Cite
this repository" button formats. For this version of the documentation:

```bibtex
@software{immunum,
  author  = {{ENPICOM}},
  title   = {immunum: antibody and T-cell receptor sequence numbering},
  version = {<!-- immunum-version -->},
  url     = {https://immunum.enpicom.com},
  note    = {https://github.com/ENPICOM/immunum}
}
```

ENPICOM. *immunum: antibody and T-cell receptor sequence numbering* (version
<!-- immunum-version -->). https://github.com/ENPICOM/immunum

## Numbering schemes

If your work depends on a particular scheme, cite its definition too:

- **IMGT**: Lefranc M-P et al. IMGT unique numbering for immunoglobulin and T cell receptor
  variable domains and Ig superfamily V-like domains. *Dev Comp Immunol* 27:55–77 (2003).
  [doi:10.1016/S0145-305X(02)00039-3](https://doi.org/10.1016/S0145-305X(02)00039-3)
- **Kabat**: Kabat EA, Wu TT, Perry HM, Gottesman KS, Foeller C. *Sequences of Proteins of
  Immunological Interest*, 5th ed. NIH Publication 91-3242 (1991).
- **Chothia**: Chothia C, Lesk AM. Canonical structures for the hypervariable regions of
  immunoglobulins. *J Mol Biol* 196:901–917 (1987).
  [doi:10.1016/0022-2836(87)90412-8](https://doi.org/10.1016/0022-2836(87)90412-8)
- **Martin** (extended Chothia): Abhinandan KR, Martin ACR. Analysis and improvements to Kabat
  and structurally correct numbering of antibody variable domains. *Mol Immunol* 45:3832–3839
  (2008). [doi:10.1016/j.molimm.2008.05.022](https://doi.org/10.1016/j.molimm.2008.05.022)
- **AHo**: Honegger A, Plückthun A. Yet another numbering scheme for immunoglobulin variable
  domains: an automatic modeling and analysis tool. *J Mol Biol* 309:657–670 (2001).
  [doi:10.1006/jmbi.2001.4662](https://doi.org/10.1006/jmbi.2001.4662)
