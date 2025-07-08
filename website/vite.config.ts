import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  // Use an empty string for the base to ensure relative paths.
  // This allows the site to be viewed correctly from the file system.
  base: '/GeneKnow/',
});
