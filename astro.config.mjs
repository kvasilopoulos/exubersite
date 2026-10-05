import { defineConfig } from "astro/config";
import tailwindcss from "@tailwindcss/vite";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";

export default defineConfig({
  site: "https://exuber.kvasilopoulos.com",
  // Methodology/settings/critical-values moved under /guide/ when the nav was
  // consolidated from 6 top-level links to 4 -- keep any inbound/bookmarked
  // links to the old paths working.
  redirects: {
    "/methodology": "/guide/methodology",
    "/settings": "/guide/critical-values",
    "/guide/settings": "/guide/critical-values",
    "/critical-values": "/guide/critical-values",
    // The suite section moved onto the landing page itself.
    "/suite": "/#suite",
  },
  markdown: {
    remarkPlugins: [remarkMath],
    rehypePlugins: [rehypeKatex],
    // Astro ships Shiki; a light theme to sit on the white ground (<Code> in
    // .astro files defaults to github-dark, so each passes theme too). The panel
    // background is overridden in global.css so it matches the palette.
    shikiConfig: { theme: "github-light", wrap: false },
  },
  vite: {
    plugins: [tailwindcss()],
  },
});
