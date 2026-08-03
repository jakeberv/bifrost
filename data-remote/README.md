# Empirical example data and provenance

`empirical-artifacts.json` is the authoritative schema-2 downloader,
provenance, and integrity manifest for the repository's empirical artifacts.
It exposes exactly these eight public identifiers:

- `jaw-tree`
- `jaw-landmarks`
- `passerine-tree`
- `passerine-traits`
- `passerine-search`
- `passerine-sensitivity`
- `passerine-posthoc`
- `simulation-preview-tables`

`bifrost_example_file()` resolves these names to checksum- and size-verified
files. On the first request without a valid cache entry, the function retrieves
the manifest and artifact currently tracked on the repository's `main` branch.
Subsequent calls reuse that verified cache and work offline. `refresh = TRUE`
bypasses the cached manifest and checks `main` for the current entry. A verified
local mirror supplied through
`BIFROST_ARTIFACT_DIR` always takes precedence, including during refreshes.

Repository maintainers must update ordinary empirical artifacts on tracked
`main` together with their manifest metadata. Validate all 11 provenance
records from the checkout root with:

```sh
python3 tools/validate-empirical-artifacts.py
```

After intentionally replacing an artifact, update its recorded checksum with:

```sh
python3 tools/validate-empirical-artifacts.py --update-checksums
```

An artifact replacement and its updated manifest entry, byte size, checksum,
and provenance must be reviewed and committed together. This prevents a new
file from becoming downloadable under stale identity or provenance metadata.

The passerine records retain their source locations in the Zenodo v1.0.0
supplementary archive. Zenodo preserves the durable versioned scientific
source and citation record; these compact repository copies support the
documented examples without replacing that archive. The manifest records the
archive's CC BY 4.0 metadata separately from the repository's GPL declaration.

The two jaw records are repository copies derived from the `R_code.zip`
analysis materials in the Troyer et al. (2025) Dryad dataset
(`10.5061/dryad.z08kprrqf`). Dryad's terms place those materials under CC0
1.0/public-domain status; scholarly citation remains recommended.

The three `vignettes/avian-skeleton/*.png` entries are provenance-only image
records and deliberately have no downloader identifiers. Published Figure 1
also retains its separate rights-review warning in the manifest.

`data-remote/` is excluded from the R source package by `.Rbuildignore`.
Generated Colab notebooks and rendered package documentation are outside this
manifest because they are generated software documentation, not empirical
inputs or scientific regression artifacts.
