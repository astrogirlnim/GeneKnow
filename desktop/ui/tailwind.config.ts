import type { Config } from 'tailwindcss'

const config: Config = {
  content: [
    './index.html',
    './src/**/*.{ts,tsx}',
    '../src-tauri/**/*.rs'
  ],
  theme: {
    extend: {},
  },
  plugins: [],
}

export default config 