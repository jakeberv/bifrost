const { test, expect } = require("@playwright/test");

async function stubExternalServices(page) {
  await page.route(
    /https:\/\/(?:cdn\.jsdelivr\.net|gc\.zgo\.at|(?:[^/]+\.)?goatcounter\.com|github\.com|raw\.githubusercontent\.com|codecov\.io|www\.r-pkg\.org|cranlogs\.r-pkg\.org|img\.shields\.io|colab\.research\.google\.com|jakeberv\.com|www\.rse\.ox\.ac\.uk)\//,
    (route) => route.abort()
  );
  await page.route("**/mermaid.esm.min.mjs*", (route) =>
    route.fulfill({
      contentType: "application/javascript",
      body: `const calls = { initialize: 0, run: 0 };
        window.mermaidCalls = calls;
        export default {
          initialize() { calls.initialize += 1; },
          run() { calls.run += 1; }
        };`
    })
  );
}

async function removeSourceLabels(page, url) {
  await stubExternalServices(page);
  await page.route(url, async (route) => {
    const response = await route.fetch();
    const html = (await response.text())
      .replace(/<([a-z][\w-]*)\b(?=[^>]*\bclass=(["'])[^"']*\bdont-index\b[^"']*\2)[^>]*>[\s\S]*?<\/\1>/gi, "")
      .replace(/<([a-z][\w-]*)\b(?=[^>]*\bclass=(["'])[^"']*\bname\b[^"']*\2)[^>]*>[\s\S]*?<\/\1>/gi, "");
    await route.fulfill({ response, body: html });
  });
}

async function serveUnderPrefix(page, prefix) {
  await page.route(`**${prefix}**`, async (route) => {
    const url = new URL(route.request().url());
    const sourceUrl = new URL(url.pathname.slice(prefix.length), url.origin);
    sourceUrl.search = url.search;
    const response = await route.fetch({ url: sourceUrl.href });
    await route.fulfill({ response });
  });
}

test("built home page includes the pkgdown extensions and analytics marker", async ({ page }) => {
  await stubExternalServices(page);
  await page.goto("/");

  await expect(page.locator('link[href$="extra.css"]')).toHaveCount(1);
  await expect(page.locator('script[src$="extra.js"]')).toHaveCount(1);
  await expect(page.locator("script[data-goatcounter]")).toHaveAttribute(
    "data-goatcounter",
    "https://bifrost.goatcounter.com/count"
  );
});

test("homepage provenance guide targets the tracked main documentation", async ({ page }) => {
  await stubExternalServices(page);
  await page.goto("/");

  await expect(page.getByRole("link", { name: "example-data guide" })).toHaveAttribute(
    "href",
    "https://github.com/jakeberv/bifrost/blob/main/data-remote/README.md"
  );
});

test("CRAN downloads chart defaults to light without exposing color-mode controls", async ({ page }) => {
  await stubExternalServices(page);
  await page.goto("/");

  const picture = page.locator("picture.cran-downloads-picture");
  const darkSource = picture.locator('source[data-theme="dark"]');
  const lightSource = picture.locator('source[data-theme="light"]');
  await expect(picture).toHaveCount(1);
  await expect(page.locator("#dropdown-lightswitch")).toHaveCount(0);
  await expect(lightSource).toHaveAttribute("media", "all");
  await expect(darkSource).toHaveAttribute("media", "not all");

  await page.evaluate(() => document.documentElement.setAttribute("data-bs-theme", "dark"));
  await expect(darkSource).toHaveAttribute("media", "all");
  await expect(lightSource).toHaveAttribute("media", "not all");

  await page.evaluate(() => document.documentElement.setAttribute("data-bs-theme", "light"));
  await expect(lightSource).toHaveAttribute("media", "all");
  await expect(darkSource).toHaveAttribute("media", "not all");

  await page.evaluate(() => document.documentElement.setAttribute("data-bs-theme", "auto"));
  await expect(darkSource).toHaveAttribute("media", "(prefers-color-scheme: dark)");
  await expect(lightSource).toHaveAttribute("media", "(prefers-color-scheme: light)");

  await page.evaluate(() => document.documentElement.removeAttribute("data-bs-theme"));
  await expect(lightSource).toHaveAttribute("media", "all");
  await expect(darkSource).toHaveAttribute("media", "not all");
});

test("vignette pages expose artifact actions and initialize Mermaid", async ({ page }) => {
  await stubExternalServices(page);
  await page.goto("/articles/quick-start-vignette.html");

  await expect(page.locator(".article-pdf-badge")).toHaveAttribute(
    "href",
    "./quick-start-vignette.pdf"
  );
  await expect(page.locator(".article-colab-badge")).toHaveAttribute(
    "href",
    "https://colab.research.google.com/github/jakeberv/bifrost/blob/main/vignettes/colab/quick-start-vignette.ipynb"
  );
  await expect.poll(() => page.evaluate(() => window.mermaidCalls?.initialize || 0)).toBe(1);
  await expect.poll(() => page.evaluate(() => window.mermaidCalls?.run || 0)).toBe(1);
});

test("vignette pages under the deployed base path expose artifact actions", async ({ page }) => {
  await stubExternalServices(page);
  await serveUnderPrefix(page, "/bifrost/");
  await page.goto("/bifrost/articles/quick-start-vignette.html");

  await expect(page.locator(".article-pdf-badge")).toHaveAttribute(
    "href",
    "./quick-start-vignette.pdf"
  );
  await expect(page.locator(".article-colab-badge")).toHaveAttribute(
    "href",
    "https://colab.research.google.com/github/jakeberv/bifrost/blob/main/vignettes/colab/quick-start-vignette.ipynb"
  );
});

test("Part 2 maximum-age slider remains a normal-flow row below the plot", async ({ page }) => {
  await stubExternalServices(page);
  await page.goto("/articles/avian-skeleton-part-2.html");

  const layout = await page.evaluate(() => {
    const plotShell = document.querySelector(
      "#lineage-decay-widget-part2 .ldw-shell > .ldw-plot-shell"
    );
    const slider = document.querySelector(
      "#lineage-decay-widget-part2 .ldw-shell > .ldw-axis-slider"
    );
    if (!plotShell || !slider) return null;

    const plotBox = plotShell.getBoundingClientRect();
    const sliderBox = slider.getBoundingClientRect();
    return {
      immediatelyAfterPlot: plotShell.nextElementSibling === slider,
      position: window.getComputedStyle(slider).position,
      plotBottom: plotBox.bottom,
      sliderTop: sliderBox.top
    };
  });

  expect(layout).not.toBeNull();
  expect(layout.immediatelyAfterPlot).toBe(true);
  expect(layout.position).toBe("static");
  expect(layout.sliderTop).toBeGreaterThanOrEqual(layout.plotBottom - 1);
});

test("website-only articles remain excluded without pkgdown source labels", async ({ page }) => {
  await removeSourceLabels(page, "/articles/development-status.html");
  const response = await page.goto("/articles/development-status.html");
  expect(response && response.status()).toBe(200);
  await expect(page.locator(".article-artifact-actions")).toHaveCount(0);
});
