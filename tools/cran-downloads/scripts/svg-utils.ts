/**
 * JSDOM lowercases camelCase SVG attribute names and element names.
 * Fix the known ones used by the adapted Star History D3 filter generation.
 */
export const fixJsdomSvgCasing = (svgContent: string): string => {
  return svgContent
    .replace(/feturbulence/g, "feTurbulence")
    .replace(/fedisplacementmap/g, "feDisplacementMap")
    .replace(/filterunits/g, "filterUnits")
    .replace(/basefrequency/g, "baseFrequency")
    .replace(/xchannelselector/g, "xChannelSelector")
    .replace(/ychannelselector/g, "yChannelSelector");
};
