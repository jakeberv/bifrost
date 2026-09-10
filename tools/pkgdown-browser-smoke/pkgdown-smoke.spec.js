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

for (const javaScriptEnabled of [false, true]) {
  test.describe(`light homepage images with JavaScript ${javaScriptEnabled ? "enabled" : "disabled"}`, () => {
    test.use({ colorScheme: "dark", javaScriptEnabled });

    test("never requests dark chart or logo variants", async ({ page }) => {
      const imageRequests = [];
      page.on("request", request => {
        if (request.resourceType() === "image") imageRequests.push(request.url());
      });
      await stubExternalServices(page);
      await page.goto("/");

      for (const [selector, imageUrl] of [
        ["picture.cran-downloads-picture", "https://raw.githubusercontent.com/jakeberv/bifrost/main/man/figures/cran-downloads.png"],
        ["picture.schmidt-sciences-picture", "https://jakeberv.com/images/SchmidtSciencesLogo.png"]
      ]) {
        const picture = page.locator(selector);
        await expect(picture).toHaveCount(1);
        await expect(picture.locator('source[data-theme="light"]')).toHaveAttribute("media", "all");
        await expect(picture.locator('source[data-theme="dark"]')).toHaveAttribute("media", "not all");
        await expect(picture.locator("img")).toHaveAttribute("src", imageUrl);
        expect(imageRequests).toContain(imageUrl);
      }
      expect(imageRequests.some(url => /(?:cran-downloads-dark|schmidt-sciences-dark)\.png/.test(url))).toBe(false);
      await expect(page.locator("#dropdown-lightswitch")).toHaveCount(0);
    });
  });
}

test("homepage images still follow page-theme changes", async ({ page }) => {
  await stubExternalServices(page);
  await page.goto("/");
  for (const theme of ["dark", "light", "auto", null]) {
    await page.evaluate(theme => {
      if (theme === null) document.documentElement.removeAttribute("data-bs-theme");
      else document.documentElement.setAttribute("data-bs-theme", theme);
    }, theme);
    for (const selector of ["picture.cran-downloads-picture", "picture.schmidt-sciences-picture"]) {
      for (const imageTheme of ["light", "dark"]) {
        const media = theme === "auto"
          ? `(prefers-color-scheme: ${imageTheme})`
          : imageTheme === (theme || "light") ? "all" : "not all";
        await expect(page.locator(`${selector} source[data-theme="${imageTheme}"]`)).toHaveAttribute("media", media);
      }
    }
  }
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
