# GeneKnow Desktop UI

This is the React + TypeScript frontend for the GeneKnow desktop application, built with Vite and integrated with Tauri for cross-platform desktop deployment.

## 🚀 Quick Start

### Prerequisites
- Node.js 20+
- pnpm (recommended package manager)

### Development Setup

1. **Install dependencies**
```bash
pnpm install
```

2. **Start development server**
```bash
pnpm dev
```
This starts the Vite development server at `http://localhost:5173`

3. **Run with Tauri (Desktop App)**
```bash
pnpm run tauri-dev
```
This starts both the frontend and the Tauri desktop application with hot reload.

## 📦 Scripts

- `pnpm dev` - Start Vite development server
- `pnpm build` - Build for production
- `pnpm lint` - Run ESLint
- `pnpm preview` - Preview production build
- `pnpm run tauri-dev` - Start Tauri development mode (frontend + desktop)
- `pnpm run tauri-build` - Build desktop application for distribution

## 🛠️ Tech Stack

- **Framework**: React 19 with TypeScript
- **Build Tool**: Vite 6.0
- **Styling**: Tailwind CSS 4.1
- **Desktop Integration**: Tauri 2.x
- **State Management**: React Router DOM 7.6
- **Real-time Communication**: Socket.IO Client

## 🏗️ Architecture

The UI communicates with the backend through multiple channels:

1. **Tauri IPC**: For desktop-specific operations (file system, OS integration)
2. **HTTP API**: For genomic pipeline operations (`http://localhost:5001`)
3. **WebSocket**: For real-time progress updates during analysis

## 📁 Project Structure

```
src/
├── api/                    # API integration layers
│   ├── geneknowPipeline.ts # GeneKnow pipeline API
│   ├── geneknowTauri.ts   # Tauri-specific API
│   └── genomicProcessing.ts # Genomic processing utilities
├── components/             # Reusable React components
│   ├── FirstRunSetup.tsx  # First-time setup wizard
│   ├── GenomicProcessing.tsx # Processing interface
│   └── Layout.tsx         # Main layout component
├── pages/                  # Page components
│   ├── DashboardPage.tsx  # Analysis results dashboard
│   ├── UploadPage.tsx     # File upload interface
│   └── WelcomePage.tsx    # Landing page
├── hooks/                  # Custom React hooks
│   └── useLogger.ts       # Logging utilities
├── App.tsx                # Main App component
└── main.tsx              # Application entry point
```

## 🔧 Configuration

- **Vite Config**: `vite.config.ts` - Build configuration
- **Tailwind Config**: `tailwind.config.ts` - Styling configuration
- **TypeScript Config**: `tsconfig.json` - TypeScript configuration
- **Tauri Config**: `../src-tauri/tauri.conf.json` - Desktop app configuration

## 🧪 Development

### Hot Reload
The development setup supports hot reload for both:
- React components (instant updates)
- Tauri backend (automatic restart)

### Environment Variables
Create a `.env` file for development:
```
VITE_API_URL=http://localhost:5001
VITE_WS_URL=http://localhost:5001
```

### API Integration
The UI integrates with the GeneKnow pipeline server running on port 5001. Make sure to start the pipeline server before testing upload functionality:

```bash
# In another terminal
cd ../../geneknow_pipeline
python enhanced_api_server.py
```

## 🚀 Production Build

Building for production creates optimized assets and a desktop application:

```bash
# Build frontend
pnpm build

# Build desktop app (includes frontend)
pnpm run tauri-build
```

The desktop application will be available in `../src-tauri/target/release/bundle/`

## 📱 Platform Support

- **macOS**: Native .dmg installer
- **Windows**: .msi installer
- **Linux**: AppImage and .deb packages

## 🔗 Related Documentation

- [Main Project README](../../README.md)
- [Tauri Integration Guide](../TESTING_GUIDE.md)
- [API Documentation](src/api/README.md)
- [Deployment Guide](../../docs/DEPLOYMENT_GUIDE.md)
