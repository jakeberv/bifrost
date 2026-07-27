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
if (!Number.isInteger(rateWindowDays) || rateWindowDays < 1) {
  throw new Error("CRAN_DOWNLOADS_RATE_WINDOW_DAYS must be a positive integer");
}
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
const outputDarkSvgPath =
  process.env.CRAN_DOWNLOADS_DARK_SVG ||
  path.join(projectRoot, "output", "cran-downloads-dark.svg");
const outputDarkPngPath =
  process.env.CRAN_DOWNLOADS_DARK_PNG ||
  path.join(projectRoot, "output", "cran-downloads-dark.png");
const readmeDarkPngPath =
  process.env.CRAN_DOWNLOADS_DARK_README_PNG ||
  path.resolve(projectRoot, "..", "..", "man", "figures", "cran-downloads-dark.png");

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

const themes = [
  {
    name: "light",
    background: "#fff",
    stroke: "#111827",
    colors: ["#dd4528", "#28a3dd"],
    svgPath: outputSvgPath,
    pngPath: outputPngPath,
    readmePngPath,
  },
  {
    name: "dark",
    background: "#0d1117",
    stroke: "white",
    colors: ["#ff6b6b", "#48dbfb"],
    svgPath: outputDarkSvgPath,
    pngPath: outputDarkPngPath,
    readmePngPath: readmeDarkPngPath,
  },
];

const hasRsvgConvert = commandExists("rsvg-convert");
for (const theme of themes) {
  const optimized = optimize(renderComposite(theme), { multipass: true }).data;
  fs.mkdirSync(path.dirname(theme.svgPath), { recursive: true });
  fs.writeFileSync(theme.svgPath, `${optimized}\n`);
  console.log(theme.svgPath);

  if (hasRsvgConvert) {
    fs.mkdirSync(path.dirname(theme.pngPath), { recursive: true });
    const png = spawnSync(
      "rsvg-convert",
      ["-b", theme.background, "-f", "png", "-o", theme.pngPath, theme.svgPath],
      { stdio: "inherit" },
    );
    if (png.status !== 0) {
      throw new Error(`rsvg-convert failed for ${theme.name} chart`);
    }
    console.log(theme.pngPath);

    fs.mkdirSync(path.dirname(theme.readmePngPath), { recursive: true });
    fs.copyFileSync(theme.pngPath, theme.readmePngPath);
    console.log(theme.readmePngPath);
  }
}

function renderComposite(theme) {
  const leftChart = renderChart({
    title: "Cumulative CRAN downloads",
    yLabel: "Downloads",
    colors: [theme.colors[0]],
    datasets: [
      {
        label: `total ${cumulative.toLocaleString("en-US")}`,
        logo: "",
        data: cumulativeData,
      },
    ],
    theme,
  });

  const rightChart = renderChart({
    title: "Download rate",
    yLabel: "Downloads/day",
    colors: [theme.colors[1]],
    legendPosition: "bottom-right",
    datasets: [
      {
        label: `${rateWindowDays}-day avg`,
        logo: "",
        data: rollingAverageData,
      },
    ],
    theme,
  });

  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${canvasWidth}" height="${canvasHeight}" viewBox="0 0 ${canvasWidth} ${canvasHeight}">
  <rect width="100%" height="100%" fill="${theme.background}"/>
  <image href="${escapeDataUri(leftChart)}" x="24" y="${chartY}" width="${chartWidth}" height="${chartHeight}"/>
  <image href="${escapeDataUri(rightChart)}" x="616" y="${chartY}" width="${chartWidth}" height="${chartHeight}"/>
</svg>
`;
}

function renderChart({ title, yLabel, datasets, colors, theme, legendPosition = "top-left" }) {
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
      theme: theme.name,
    },
    {
      envType: "node",
      xTickLabelType: "Date",
      chartWidth,
      chartHeight,
      xTickCount: 4,
      yTickCount: 4,
      dataColors: colors,
      backgroundColor: theme.background,
      strokeColor: theme.stroke,
      legendPosition,
    },
  );

  const background = dom.window.document.createElement("rect");
  background.setAttribute("width", "100%");
  background.setAttribute("height", "100%");
  background.setAttribute("fill", theme.background);
  svg.insertBefore(background, svg.firstChild);

  return fixJsdomSvgCasing(svg.outerHTML);
}

function escapeDataUri(content) {
  return `data:image/svg+xml;base64,${Buffer.from(content, "utf8").toString("base64")}`;
}

function commandExists(command) {
  return spawnSync("sh", ["-lc", `command -v ${command}`], { stdio: "ignore" }).status === 0;
}
