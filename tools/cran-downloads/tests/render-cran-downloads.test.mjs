import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { spawnSync } from "node:child_process";
import test from "node:test";
import { fileURLToPath } from "node:url";

const projectRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");

test("rejects invalid rolling-average window sizes", async (t) => {
  for (const value of ["0", "-1", "not-a-number", "1.5"]) {
    await t.test(value, (subtest) => {
      const tempDir = fs.mkdtempSync(path.join(os.tmpdir(), "bifrost-cran-downloads-"));
      subtest.after(() => fs.rmSync(tempDir, { recursive: true, force: true }));
      const dataPath = path.join(tempDir, "downloads.json");
      fs.writeFileSync(
        dataPath,
        JSON.stringify({
          downloads: [
            { day: "2026-01-01", downloads: 10 },
            { day: "2026-01-02", downloads: 20 },
            { day: "2026-01-03", downloads: 30 },
          ],
        }),
      );

      const result = spawnSync("pnpm", ["exec", "tsx", "scripts/render-cran-downloads.mjs"], {
        cwd: projectRoot,
        encoding: "utf8",
        env: {
          ...process.env,
          CRAN_DOWNLOADS_DATA: dataPath,
          CRAN_DOWNLOADS_RATE_WINDOW_DAYS: value,
          CRAN_DOWNLOADS_SVG: path.join(tempDir, "chart.svg"),
          CRAN_DOWNLOADS_PNG: path.join(tempDir, "chart.png"),
          CRAN_DOWNLOADS_README_PNG: path.join(tempDir, "readme.png"),
        },
      });

      assert.notEqual(result.status, 0, `expected ${value} to fail`);
      assert.match(
        result.stderr,
        /CRAN_DOWNLOADS_RATE_WINDOW_DAYS must be a positive integer/,
      );
    });
  }
});
