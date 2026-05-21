import fs from "node:fs";
import path from "node:path";
import { spawnSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import { JSDOM } from "jsdom";
import { optimize } from "svgo";
import XYChart from "../vendor/star-history/shared/packages/xy-chart.tsx";
import { fixJsdomSvgCasing } from "./svg-utils.ts";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const projectRoot = path.resolve(__dirname, "..");

const packageName = process.env.CRAN_PACKAGE || "bifrost";
const rateWindowDays = Number(process.env.CRAN_DOWNLOADS_RATE_WINDOW_DAYS || 14);
const dataPath =
  process.env.CRAN_DOWNLOADS_DATA ||
  path.join(projectRoot, "data", `${packageName}-cran-downloads.json`);
const outputSvgPath =
  process.env.CRAN_DOWNLOADS_SVG ||
  path.join(projectRoot, "output", "cran-downloads.svg");
const outputPngPath =
  process.env.CRAN_DOWNLOADS_PNG ||
  path.join(projectRoot, "output", "cran-downloads.png");
const readmePngPath =
  process.env.CRAN_DOWNLOADS_README_PNG ||
  path.resolve(projectRoot, "..", "..", "man", "figures", "cran-downloads.png");

const canvasWidth = 1200;
const canvasHeight = 500;
const chartWidth = 560;
const chartHeight = 430;
const chartY = 28;

const payload = JSON.parse(fs.readFileSync(dataPath, "utf8"));
const rows = (payload.downloads || payload[0]?.downloads || [])
  .map((row) => ({
    day: row.day,
    downloads: Number(row.downloads || 0),
  }))
  .filter((row) => Number.isFinite(row.downloads) && !Number.isNaN(new Date(row.day).getTime()))
  .sort((a, b) => a.day.localeCompare(b.day));

if (rows.length === 0) {
  throw new Error(`No CRAN download rows found in ${dataPath}`);
}

let cumulative = 0;
const cumulativeData = rows.map((row) => {
  cumulative += row.downloads;
  return {
    x: new Date(row.day),
    y: cumulative,
  };
});

const rollingAverageData = rows.map((row, index) => {
  const window = rows.slice(Math.max(0, index - rateWindowDays + 1), index + 1);
  const average = window.reduce((sum, item) => sum + item.downloads, 0) / window.length;
  return {
    x: new Date(row.day),
    y: Math.round(average * 10) / 10,
  };
});

const leftChart = renderChart({
  title: "Cumulative CRAN downloads",
  yLabel: "Downloads",
  colors: ["#dd4528"],
  datasets: [
    {
      label: `total ${cumulative.toLocaleString("en-US")}`,
      logo: "",
      data: cumulativeData,
    },
  ],
});

const rightChart = renderChart({
  title: "Download rate",
  yLabel: "Downloads/day",
  colors: ["#28a3dd"],
  legendPosition: "bottom-right",
  datasets: [
    {
      label: `${rateWindowDays}-day avg`,
      logo: "",
      data: rollingAverageData,
    },
  ],
});

const composite = `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${canvasWidth}" height="${canvasHeight}" viewBox="0 0 ${canvasWidth} ${canvasHeight}">
  <rect width="100%" height="100%" fill="white"/>
  <image href="${escapeDataUri(leftChart)}" x="24" y="${chartY}" width="${chartWidth}" height="${chartHeight}"/>
  <image href="${escapeDataUri(rightChart)}" x="616" y="${chartY}" width="${chartWidth}" height="${chartHeight}"/>
</svg>
`;

const optimized = optimize(composite, { multipass: true }).data;
fs.mkdirSync(path.dirname(outputSvgPath), { recursive: true });
fs.writeFileSync(outputSvgPath, `${optimized}\n`);
console.log(outputSvgPath);

if (commandExists("rsvg-convert")) {
  const png = spawnSync("rsvg-convert", ["-b", "white", "-f", "png", "-o", outputPngPath, outputSvgPath], {
    stdio: "inherit",
  });
  if (png.status !== 0) {
    throw new Error("rsvg-convert failed");
  }
  console.log(outputPngPath);

  fs.mkdirSync(path.dirname(readmePngPath), { recursive: true });
  fs.copyFileSync(outputPngPath, readmePngPath);
  console.log(readmePngPath);
}

function renderChart({ title, yLabel, datasets, colors, legendPosition = "top-left" }) {
  const dom = new JSDOM(`<!DOCTYPE html><body></body>`);
  const svg = dom.window.document.createElement("svg");

  svg.setAttribute("width", String(chartWidth));
  svg.setAttribute("height", String(chartHeight));
  svg.setAttribute("viewBox", `0 0 ${chartWidth} ${chartHeight}`);
  svg.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  dom.window.document.body.append(svg);

  XYChart(
    svg,
    {
      title,
      xLabel: "Date",
      yLabel,
      data: { datasets },
      showDots: false,
      transparent: false,
      theme: "light",
    },
    {
      envType: "node",
      xTickLabelType: "Date",
      chartWidth,
      chartHeight,
      xTickCount: 4,
      yTickCount: 4,
      dataColors: colors,
      backgroundColor: "white",
      strokeColor: "#111827",
      legendPosition,
    },
  );

  const background = dom.window.document.createElement("rect");
  background.setAttribute("width", "100%");
  background.setAttribute("height", "100%");
  background.setAttribute("fill", "white");
  svg.insertBefore(background, svg.firstChild);

  return fixJsdomSvgCasing(svg.outerHTML);
}

function escapeDataUri(content) {
  return `data:image/svg+xml;base64,${Buffer.from(content, "utf8").toString("base64")}`;
}

function commandExists(command) {
  return spawnSync("sh", ["-lc", `command -v ${command}`], { stdio: "ignore" }).status === 0;
}
