# CRAN Downloads Dark Mode

## Goal

Render matching light and dark CRAN download charts and automatically show the appropriate version in both GitHub and pkgdown. Preserve the approved dark palette: `#0d1117` background, white labels and axes, coral cumulative series, and cyan rate series.

## Design

- Compute the download series once, then render light and dark theme configurations without duplicating chart logic.
- Keep the current light outputs and add `cran-downloads-dark.svg` and `cran-downloads-dark.png` under `output/`.
- Copy both PNG variants to `man/figures/`. Add `CRAN_DOWNLOADS_DARK_SVG`, `CRAN_DOWNLOADS_DARK_PNG`, and `CRAN_DOWNLOADS_DARK_README_PNG` for custom paths.
- Use GitHub's supported `<picture>` and `prefers-color-scheme` pattern in the README, with the light image as fallback.
- Enable pkgdown's Bootstrap 5 light switch. Extend `pkgdown/extra.js` to align the picture with pkgdown's current `data-bs-theme`, including manual overrides.
- Track both variants in the weekly workflow and update the tracker documentation.

## Verification

- Renderer tests verify both themes and their expected colors; existing invalid-window tests remain.
- The pkgdown browser smoke test covers light, dark, and automatic selection.
- A render or conversion failure prevents the workflow from committing partial output.

## Acceptance

GitHub and pkgdown select the appropriate chart; pkgdown theme changes update it without reloading; weekly automation updates both variants together.
