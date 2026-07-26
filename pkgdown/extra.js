const ARTICLE_ARTIFACT_SLUGS = new Set([
  "avian-skeleton-part-1",
  "avian-skeleton-part-2",
  "avian-skeleton-part-3",
  "avian-skeleton-part-4",
  "avian-skeleton-part-5",
  "jaw-shape-vignette",
  "pca-model-selection-and-bifrost-vignette",
  "quick-start-vignette",
  "rate-map-jaw-shape-part-2-comparisons",
  "rate-map-jaw-shape-vignette",
  "simulation-study-part-1",
  "simulation-study-part-2",
  "theoretical-background-vignette"
]);

function getArticleArtifactSlug(pathname) {
  const match = /^(?:\/[^/]+)*\/articles\/([^/]+)\.html$/.exec(pathname);
  if (!match) return null;

  let slug;
  try {
    slug = decodeURIComponent(match[1]);
  } catch {
    return null;
  }
  return ARTICLE_ARTIFACT_SLUGS.has(slug) ? slug : null;
}

// Load a pinned, stable Mermaid (v10) and force it to be the one we use.
import('https://cdn.jsdelivr.net/npm/mermaid@10.9.1/dist/mermaid.esm.min.mjs?v=1091')
  .then(({ default: mermaid }) => {
    window.mermaid = mermaid;

    (function initMermaid() {
      const run = () => {
        // Find mermaid code blocks produced by pkgdown/pandoc (covers multiple shapes)
        const sels = [
          'pre code.language-mermaid',
          'pre code.mermaid',
          'pre.mermaid',
          'code.mermaid',
          'div.sourceCode pre code.language-mermaid',
          'div.sourceCode pre code.mermaid'
        ];
        const nodes = new Set();
        sels.forEach(sel => document.querySelectorAll(sel).forEach(el => nodes.add(el)));

        nodes.forEach(code => {
          const pre = code.closest('pre') || code;
          const div = document.createElement('div');
          div.className = 'mermaid';
          div.textContent = (code.textContent || code.innerText || '').trim(); // important
          pre.replaceWith(div);
        });

        mermaid.initialize({ startOnLoad: false, theme: "neutral" });
        mermaid.run({ querySelector: '.mermaid' });
      };

      if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', run);
      } else {
        run();
      }
    })();
  });

(function addArticleArtifacts() {
  const run = () => {
    const slug = getArticleArtifactSlug(window.location.pathname);
    if (!slug) return;

    const article = document.querySelector('.template-article');
    const header = article && article.querySelector('main .page-header');
    if (!header || header.querySelector('.article-artifact-actions')) return;

    const encodedSlug = encodeURIComponent(slug);

    const actions = document.createElement('div');
    actions.className = 'article-artifact-actions';
    actions.setAttribute('role', 'group');
    actions.setAttribute('aria-label', 'Article downloads and notebooks');

    const pdf = document.createElement('a');
    pdf.className = 'article-pdf-badge';
    pdf.href = './' + encodedSlug + '.pdf';
    pdf.title = 'Download PDF';
    pdf.setAttribute('aria-label', 'Download PDF');

    const pdfImage = document.createElement('img');
    pdfImage.src = 'https://img.shields.io/badge/Download%20as%20PDF-EF3939?style=flat-square&logo=adobeacrobatreader&logoColor=white&labelColor=ec1c24';
    pdfImage.alt = 'Download as PDF';
    pdf.appendChild(pdfImage);

    const colab = document.createElement('a');
    colab.className = 'article-colab-badge';
    colab.href = 'https://colab.research.google.com/github/jakeberv/bifrost/blob/main/vignettes/colab/' +
      encodedSlug + '.ipynb';
    colab.target = '_blank';
    colab.rel = 'noopener noreferrer';
    colab.title = 'Open in Colab';
    colab.setAttribute('aria-label', 'Open in Colab');

    const colabImage = document.createElement('img');
    colabImage.src = 'https://colab.research.google.com/assets/colab-badge.svg';
    colabImage.alt = 'Open In Colab';
    colab.appendChild(colabImage);

    actions.append(pdf, colab);
    header.appendChild(actions);
  };

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', run);
  } else {
    run();
  }
})();
