import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { fileURLToPath } from "node:url";

const projectRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");

function createFixture(t) {
  const tempDir = fs.mkdtempSync(path.join(os.tmpdir(), "bifrost-cran-downloads-"));
  t.after(() => fs.rmSync(tempDir, { recursive: true, force: true }));
  const paths = {
    data: path.join(tempDir, "downloads.json"),
    lightSvg: path.join(tempDir, "chart.svg"),
    lightPng: path.join(tempDir, "chart.png"),
    lightReadmePng: path.join(tempDir, "readme.png"),
    darkSvg: path.join(tempDir, "chart-dark.svg"),
    darkPng: path.join(tempDir, "chart-dark.png"),
    darkReadmePng: path.join(tempDir, "readme-dark.png"),
  };
  fs.writeFileSync(
    paths.data,
    JSON.stringify({
      downloads: [
        { day: "2026-01-01", downloads: 10 },
        { day: "2026-01-02", downloads: 20 },
        { day: "2026-01-03", downloads: 30 },
      ],
    }),
  );
  return paths;
}

function runRenderer(paths, extraEnv = {}) {
  return spawnSync("pnpm", ["exec", "tsx", "scripts/render-cran-downloads.mjs"], {
    cwd: projectRoot,
    encoding: "utf8",
    env: {
      ...process.env,
      CRAN_DOWNLOADS_DATA: paths.data,
      CRAN_DOWNLOADS_SVG: paths.lightSvg,
      CRAN_DOWNLOADS_PNG: paths.lightPng,
      CRAN_DOWNLOADS_README_PNG: paths.lightReadmePng,
      CRAN_DOWNLOADS_DARK_SVG: paths.darkSvg,
      CRAN_DOWNLOADS_DARK_PNG: paths.darkPng,
      CRAN_DOWNLOADS_DARK_README_PNG: paths.darkReadmePng,
      ...extraEnv,
    },
  });
}

function embeddedSvgContent(composite) {
  return [...composite.matchAll(/href="data:image\/svg\+xml;base64,([^"]+)"/g)]
    .map((match) => Buffer.from(match[1], "base64").toString("utf8"))
    .join("\n");
}

test("rejects invalid rolling-average window sizes", async (t) => {
  for (const value of ["0", "-1", "not-a-number", "1.5"]) {
    await t.test(value, (subtest) => {
      const paths = createFixture(subtest);
      const result = runRenderer(paths, { CRAN_DOWNLOADS_RATE_WINDOW_DAYS: value });

      assert.notEqual(result.status, 0, `expected ${value} to fail`);
      assert.match(
        result.stderr,
        /CRAN_DOWNLOADS_RATE_WINDOW_DAYS must be a positive integer/,
      );
    });
  }
});

test("renders light and dark chart variants", (t) => {
  const paths = createFixture(t);
  const result = runRenderer(paths);

  assert.equal(result.status, 0, result.stderr);
  assert.ok(fs.existsSync(paths.darkSvg), "expected dark SVG output");

  const lightSvg = fs.readFileSync(paths.lightSvg, "utf8");
  const darkSvg = fs.readFileSync(paths.darkSvg, "utf8");
  assert.match(lightSvg, /fill="#fff"/);
  assert.match(darkSvg, /fill="#0d1117"/);
  assert.match(embeddedSvgContent(lightSvg), /#dd4528/);
  assert.match(embeddedSvgContent(lightSvg), /#28a3dd/);
  assert.match(embeddedSvgContent(darkSvg), /#ff6b6b/);
  assert.match(embeddedSvgContent(darkSvg), /#48dbfb/);
});
