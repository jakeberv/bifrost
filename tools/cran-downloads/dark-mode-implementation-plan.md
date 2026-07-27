# CRAN Downloads Dark Mode Implementation Plan

> Use the checkboxes below to track implementation and verification.

**Goal:** Generate and automatically select matching light and dark CRAN download charts on GitHub and pkgdown.

**Architecture:** Compute the chart data once and render it through two explicit theme configurations. GitHub selects with `<picture>` media sources; pkgdown adds its native light switch and synchronizes those sources with Bootstrap's active theme.

**Tech Stack:** Node.js, tsx, JSDOM/D3 SVG rendering, rsvg-convert, pkgdown/Bootstrap 5, Playwright.

## Global Constraints

- Preserve the current light chart and the approved `#0d1117`/white/coral/cyan dark palette.
- Add no dependencies and keep PNG as the displayed format.
- A failed render or conversion must prevent partial automation commits.

---

### Task 1: Render Both Themes

**Files:**
- Modify: `tools/cran-downloads/tests/render-cran-downloads.test.mjs`
- Modify: `tools/cran-downloads/scripts/render-cran-downloads.mjs`
- Modify: `tools/cran-downloads/README.md`
- Modify: `.github/workflows/update-cran-downloads.yml`
- Create: `tools/cran-downloads/output/cran-downloads-dark.svg`
- Create: `tools/cran-downloads/output/cran-downloads-dark.png`
- Create: `man/figures/cran-downloads-dark.png`

**Interfaces:**
- Produces `CRAN_DOWNLOADS_DARK_SVG`, `CRAN_DOWNLOADS_DARK_PNG`, and `CRAN_DOWNLOADS_DARK_README_PNG` path overrides.

- [ ] Add a fixture-based test that runs the real renderer with light and dark temporary paths and asserts both SVGs exist with `#fff` and `#0d1117` outer backgrounds.
- [ ] Run `pnpm test` in `tools/cran-downloads`; expect the new test to fail because no dark output exists.
- [ ] Replace the single hard-coded render with light/dark theme objects and a shared `renderComposite(theme)` path:

```js
const themes = [
  { name: "light", background: "white", stroke: "#111827", colors: ["#dd4528", "#28a3dd"] },
  { name: "dark", background: "#0d1117", stroke: "white", colors: ["#ff6b6b", "#48dbfb"] },
];
```

- [ ] Render and convert both variants, copying each PNG to `man/figures/`; throw on either conversion failure.
- [ ] Run `pnpm test`; expect all renderer tests to pass.
- [ ] Generate the dark artifacts, document them in the tracker README, and add them to the workflow's `paths` array.
- [ ] Commit the renderer slice with `git commit -m "Add dark CRAN downloads chart"`.

### Task 2: Switch Themes on GitHub and pkgdown

**Files:**
- Modify: `README.md`
- Modify: `_pkgdown.yml`
- Modify: `pkgdown/extra.js`
- Modify: `tools/pkgdown-browser-smoke/pkgdown-smoke.spec.js`

**Interfaces:**
- Consumes `.cran-downloads-picture` with `source[data-theme="light|dark"]`.
- Produces source `media="all"` for the active pkgdown theme and `media="not all"` for the inactive theme.

- [ ] Add a Playwright test that asserts pkgdown exposes no theme control while retaining light, dark, and automatic source behavior.
- [ ] Build/serve pkgdown and run the focused smoke test; expect failure because the picture and synchronizer do not exist.
- [ ] Replace the README image with GitHub's supported markup:

```html
<picture class="cran-downloads-picture">
  <source data-theme="dark" media="(prefers-color-scheme: dark)" srcset="man/figures/cran-downloads-dark.png">
  <source data-theme="light" media="(prefers-color-scheme: light)" srcset="man/figures/cran-downloads.png">
  <img src="man/figures/cran-downloads.png" alt="Cumulative CRAN downloads for bifrost" width="560">
</picture>
```

- [ ] Leave pkgdown's `template.light-switch` disabled so the site uses its default light mode.
- [ ] Add a DOM-ready synchronizer and `MutationObserver` in `pkgdown/extra.js` that updates the two source media attributes from the root `data-bs-theme` value.
- [ ] Rebuild pkgdown and run the focused smoke test; expect no visible selector and preserved light, dark, and auto behavior.
- [ ] Commit the display slice with `git commit -m "Switch CRAN chart with site theme"`.

### Task 3: Verify and Publish

**Files:** All files above.

- [ ] Run `pnpm test` in `tools/cran-downloads`.
- [ ] Run the pkgdown browser smoke suite.
- [ ] Run `git diff --check` and confirm only intended artifacts and source files changed.
- [ ] Push `automation/cran-downloads` and confirm PR #209 points at the new head.
