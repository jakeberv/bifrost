# CRAN Downloads Chart

This directory generates the CRAN download chart used for `bifrost`.

## Data Source

Download counts are fetched from the [cranlogs](https://cranlogs.r-pkg.org/) API, which reports daily package download counts from RStudio CRAN mirror logs.

Counts are download events, not unique users.

`data/bifrost-cran-downloads.json` is the repository's canonical stored copy of the daily series. The fetch script merges new API results into this file, preserving previously archived days and overwriting only matching days returned by cranlogs so upstream corrections are retained.

## Generated Series

The chart contains:

- Cumulative CRAN downloads over all available recorded days.
- A trailing 14-day moving average of daily downloads, stepped daily.

## Commands

```sh
pnpm install
pnpm run update
```

`pnpm run update` writes:

- `data/bifrost-cran-downloads.json`
- `output/cran-downloads.svg` as the SVG render source
- `output/cran-downloads.png` as the rendered PNG when `rsvg-convert` is available
- `../../man/figures/cran-downloads.png` as the README/pkgdown copy

Use the PNG in Markdown. The SVG is retained as a render source, but PNG output avoids renderer/font differences across GitHub, pkgdown, and local preview tools.

## Automation

`.github/workflows/update-cran-downloads.yml` updates the stored JSON and chart once per week.

Each run rebuilds the persistent `automation/cran-downloads` branch from the current default branch, regenerates the data and chart, and force-updates that branch only when the generated artifacts differ. This keeps a single standing automated pull request with:

- `data/bifrost-cran-downloads.json`
- `output/cran-downloads.svg`
- `output/cran-downloads.png`
- `../../man/figures/cran-downloads.png`

The workflow uses a pull request instead of pushing directly to `main`, so it works with branch protection rules. Merging the PR makes the stored JSON on `main` the canonical archive for the next run.

## Attribution

The chart renderer is adapted from [star-history/star-history](https://github.com/star-history/star-history), which is MIT licensed. The upstream license is included at `vendor/star-history/LICENSE`.
