// Post-build check over dist/: every internal link resolves to a built page,
// and every fragment resolves to an id on that page. The replication corpus is
// vendored from docs/enhancements/ with its links rewritten, so this is what
// catches an anchor that drifted upstream.
//
//   npm run build && node scripts/check-links.mjs
import fs from "node:fs";
import path from "node:path";

const dist = "dist";
const pages = {};

const walk = (dir) =>
  fs.readdirSync(dir, { withFileTypes: true }).forEach((e) => {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) return walk(p);
    if (!e.name.endsWith(".html")) return;
    const rel = path.relative(dist, p).split(path.sep).join("/");
    const route = "/" + rel.replace(/index\.html$/, "").replace(/\/$/, "");
    pages[route] = fs.readFileSync(p, "utf8");
  });
walk(dist);

const ids = Object.fromEntries(
  Object.entries(pages).map(([r, h]) => [r, new Set([...h.matchAll(/\sid="([^"]+)"/g)].map((m) => m[1]))]),
);

const norm = (href, from) => (href ? href.replace(/\/$/, "") || "/" : from);
let broken = 0;

for (const [route, html] of Object.entries(pages)) {
  for (const [, href, frag] of html.matchAll(/href="(\/[^"#]*)?#([^"]+)"/g)) {
    const target = norm(href, route);
    if (!(target in pages)) {
      console.error(`missing page   ${route} -> ${target}`);
      broken++;
    } else if (!ids[target].has(decodeURIComponent(frag))) {
      console.error(`missing anchor ${route} -> ${target}#${decodeURIComponent(frag)}`);
      broken++;
    }
  }
  for (const [, href] of html.matchAll(/href="(\/[^"#]*)"/g)) {
    const target = norm(href, route);
    if (!(target in pages) && !/^\/(favicon|_astro)/.test(target)) {
      console.error(`missing link   ${route} -> ${target}`);
      broken++;
    }
  }
}

console.log(`${Object.keys(pages).length} routes, ${broken} broken references`);
process.exit(broken ? 1 : 0);
