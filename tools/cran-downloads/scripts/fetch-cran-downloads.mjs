import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const projectRoot = path.resolve(__dirname, "..");

const packageName = process.env.CRAN_PACKAGE || "bifrost";
const startDate = process.env.CRAN_DOWNLOADS_START || "1900-01-01";
const endDate = process.env.CRAN_DOWNLOADS_END || yesterdayUtc();
const outputPath =
  process.env.CRAN_DOWNLOADS_DATA ||
  path.join(projectRoot, "data", `${packageName}-cran-downloads.json`);

const endpoint = `https://cranlogs.r-pkg.org/downloads/daily/${startDate}:${endDate}/${encodeURIComponent(
  packageName,
)}`;

const payload = await fetchJson(endpoint);
const result = Array.isArray(payload) ? payload[0] : payload;
const rows = Array.isArray(result?.downloads) ? result.downloads : [];

if (rows.length === 0) {
  throw new Error(`cranlogs returned no daily download rows for ${packageName}`);
}

const fetchedDownloads = normalizeRows(rows);
const existingDownloads = readExistingDownloads(outputPath);
const downloadsByDay = new Map(existingDownloads.map((row) => [row.day, row]));

// Preserve previously archived days, while allowing cranlogs to correct any
// day that is present in the newest response.
for (const row of fetchedDownloads) {
  downloadsByDay.set(row.day, row);
}

const downloads = [...downloadsByDay.values()].sort((a, b) => a.day.localeCompare(b.day));

const body = {
  package: packageName,
  source: "https://cranlogs.r-pkg.org/",
  source_endpoint_template: "https://cranlogs.r-pkg.org/downloads/daily/{start}:{end}/{package}",
  source_note: "Daily package download counts from RStudio CRAN mirror logs; counts are download events, not unique users.",
  archive_note: "This file is the repository's canonical CRAN downloads archive. Existing days are preserved; fetched rows overwrite matching days so upstream corrections are retained.",
  query_start: startDate,
  start: downloads[0].day,
  end: downloads.at(-1).day,
  total_downloads: downloads.reduce((sum, row) => sum + row.downloads, 0),
  downloads,
};

fs.mkdirSync(path.dirname(outputPath), { recursive: true });
fs.writeFileSync(outputPath, `${JSON.stringify(body, null, 2)}\n`);
console.log(outputPath);

function normalizeRows(sourceRows) {
  return sourceRows
  .map((row) => ({
    day: row.day,
    downloads: Number(row.downloads || 0),
  }))
  .filter((row) => Number.isFinite(row.downloads) && !Number.isNaN(new Date(row.day).getTime()))
  .sort((a, b) => a.day.localeCompare(b.day));
}

function readExistingDownloads(filePath) {
  if (!fs.existsSync(filePath)) {
    return [];
  }

  const existing = JSON.parse(fs.readFileSync(filePath, "utf8"));
  return normalizeRows(existing.downloads || existing[0]?.downloads || []);
}

function yesterdayUtc() {
  const now = new Date();
  const yesterday = new Date(Date.UTC(now.getUTCFullYear(), now.getUTCMonth(), now.getUTCDate() - 1));
  return yesterday.toISOString().slice(0, 10);
}

async function fetchJson(url) {
  try {
    const response = await fetch(url);
    if (!response.ok) {
      throw new Error(`cranlogs request failed: ${response.status} ${response.statusText}`);
    }
    return await response.json();
  } catch (error) {
    const curl = spawnSync("curl", ["-fsSL", url], { encoding: "utf8" });
    if (curl.status === 0) {
      return JSON.parse(curl.stdout);
    }
    throw new Error(`cranlogs request failed via fetch and curl fallback: ${error.message}`);
  }
}
