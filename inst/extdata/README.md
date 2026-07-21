# Empirical artifact provenance

`empirical-artifacts.json` records the repository-relative path, SHA-256
checksum, cited source and version, transformation provenance, and known reuse
status for the empirical artifacts added with the Berv et al. (2026)
vignettes. Its scope is deliberately limited to these directories:

- `inst/extdata/avian-skeleton`
- `inst/extdata/simulation-study-cache`
- `vignettes/avian-skeleton`

In a repository checkout, validate the manifest from the checkout root with:

```sh
python3 tools/validate-empirical-artifacts.py
```

The validator is a repository-maintenance tool. The built R source package
installs this provenance manifest and the package extdata, but excludes
`tools/` and the repository-only vignette image artifacts via `.Rbuildignore`.

After intentionally replacing an artifact, review its provenance fields and
then regenerate its recorded checksum with:

```sh
python3 tools/validate-empirical-artifacts.py --update-checksums
```

The package declares `GPL (>= 2)` in `DESCRIPTION`. The Zenodo API record for
the v1.0.0 supplementary archive identifies its license as CC BY 4.0; the
manifest cites that API record and assigns the archive-derived artifacts to a
separate `zenodo_cc_by_4_0` license record. The API metadata were inspected
without downloading the 15.2 GB archive.

The published Nature Figure 1 is intentionally separate. Neither the
repository metadata nor the article metadata recorded in this manifest
establishes figure-specific reuse permission. Consult the article rights
statement and rights holder before redistributing that figure outside the
package.

The manifest does not cover generated Colab notebooks or rendered package
documentation because those are generated software documentation, not
empirical inputs or scientific regression artifacts.
