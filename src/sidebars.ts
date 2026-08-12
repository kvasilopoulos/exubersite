import { getCollection } from "astro:content";
import type { SidebarSection } from "./layouts/DocsLayout.astro";
import reference from "./data/reference.json";

export const GUIDE_PAGES = [
  { href: "/guide/methodology", label: "Methodology" },
  { href: "/guide/settings", label: "Settings" },
  { href: "/guide/critical-values", label: "Critical values" },
];

export async function replicationSection(): Promise<SidebarSection> {
  const pages = (await getCollection("replication")).sort((a, b) => a.data.order - b.data.order);
  return {
    title: "Replication",
    items: [
      { href: "/replication", label: "Overview" },
      ...pages.map((p) => ({ href: `/replication/${p.id}`, label: p.data.title })),
    ],
  };
}

export function referenceSection(): SidebarSection {
  return {
    title: "Reference",
    items: [
      { href: "/reference", label: "Index" },
      ...reference.r.map((g) => ({
        href: `/reference#${slugify(g.title)}`,
        label: g.title,
      })),
      { href: "/reference#python", label: "Python — pyexuber" },
      { href: "/reference#cpp", label: "C++ — exubercore" },
    ],
  };
}

export function guideSection(): SidebarSection {
  return { title: "Guide", items: [{ href: "/guide", label: "Overview" }, ...GUIDE_PAGES] };
}

export function slugify(s: string): string {
  return s
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-|-$/g, "");
}
