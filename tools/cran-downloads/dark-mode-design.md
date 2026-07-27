# CRAN Downloads Dark Mode

## Goal

Render matching light and dark CRAN download charts and automatically show the appropriate version in GitHub. Keep pkgdown in its default light mode without exposing a theme control, while preserving the existing switching mechanism for future use. Preserve the approved dark palette: `#0d1117` background, white labels and axes, coral cumulative series, and cyan rate series.

## Design

- Compute the download series once, then render light and dark theme configurations without duplicating chart logic.
- Keep the current light outputs and add `cran-downloads-dark.svg` and `cran-downloads-dark.png` under `output/`.
- Copy both PNG variants to `man/figures/`. Add `CRAN_DOWNLOADS_DARK_SVG`, `CRAN_DOWNLOADS_DARK_PNG`, and `CRAN_DOWNLOADS_DARK_README_PNG` for custom paths.
- Use GitHub's supported `<picture>` and `prefers-color-scheme` pattern in the README, with the light image as fallback.
- Keep pkgdown's Bootstrap light switch disabled so the site remains in its default light mode without a visible selector. Default the existing `pkgdown/extra.js` synchronizer to light while preserving explicit light, dark, and auto handling.
- Track both variants in the weekly workflow and update the tracker documentation.

## Verification

- Renderer tests verify both themes and their expected colors; existing invalid-window tests remain.
- The pkgdown browser smoke test verifies that no theme selector is rendered, light is the default, and explicit light, dark, and automatic source behavior remains available.
- A render or conversion failure prevents the workflow from committing partial output.

## Acceptance

GitHub selects the chart from the system preference; pkgdown exposes no theme selector and uses its default light mode; weekly automation updates both variants together.
