import { getCollection } from "astro:content";
import type { SidebarSection } from "./layouts/DocsLayout.astro";
import reference from "./data/reference.json";

export const GUIDE_PAGES = [
  { href: "/guide", label: "Getting started" },
  { href: "/guide/methodology", label: "Methodology" },
  { href: "/guide/critical-values", label: "Settings and critical values" },
  { href: "/guide/pipeline", label: "Results, tidying and plotting" },
  
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
      { href: "/reference#python", label: "Python (pyexuber)" },
      { href: "/reference#cpp", label: "C++ (exubercore)" },
    ],
  };
}

export function guideSection(): SidebarSection {
  return { title: "Introduction", items: GUIDE_PAGES };
}

// The method families beyond PSY: user guides built from the exuber vignettes.
export async function beyondPsySection(): Promise<SidebarSection> {
  const pages = (await getCollection("guide"))
    .filter((p) => p.data.group === "Beyond PSY")
    .sort((a, b) => a.data.order - b.data.order);
  return {
    title: "Beyond PSY",
    items: pages.map((p) => ({ href: `/guide/${p.id}`, label: p.data.title })),
  };
}

export function slugify(s: string): string {
  return s
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-|-$/g, "");
}
