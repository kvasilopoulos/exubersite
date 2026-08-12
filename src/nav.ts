// Single source of truth for the header. Active state is derived from the
// pathname, so adding a page here is the only edit needed.
export const links = [
  { href: "/methodology", label: "Methodology" },
  { href: "/settings", label: "Settings" },
  { href: "/critical-values", label: "Critical values" },
  { href: "/replication", label: "Replication" },
  { href: "/reference", label: "Reference" },
  { href: "/suite", label: "Suite" },
];

export const repo = "https://github.com/kvasilopoulos/exuber";

export const isActive = (pathname: string, href: string) => {
  const p = pathname.replace(/\/$/, "") || "/";
  return p === href || p.startsWith(href + "/");
};
