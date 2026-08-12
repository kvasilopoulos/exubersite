import { defineConfig } from "astro/config";
import tailwindcss from "@tailwindcss/vite";

export default defineConfig({
  site: "https://exuber.kvasilopoulos.com",
  markdown: {
    // Astro ships Shiki; a light theme to sit on the white ground. The panel
    // background is overridden in global.css so it matches the palette.
    shikiConfig: { theme: "github-light", wrap: false },
  },
  vite: {
    plugins: [tailwindcss()],
  },
});
