#!/bin/bash

# GenePredict Quick Start Script
# Simple wrapper around the main dev setup script with sensible defaults

set -e

# Colors for output
readonly GREEN='\033[0;32m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m' # No Color

echo -e "${BLUE}🧬 GenePredict Quick Start${NC}"
echo -e "${BLUE}=========================${NC}"
echo ""
echo -e "${GREEN}This script will:${NC}"
echo "  ✅ Check prerequisites (Node.js, Python, Rust)"
echo "  ✅ Set up development environment"
echo "  ✅ Install all dependencies"
echo "  ✅ Start Tauri desktop application"
echo ""
echo -e "${BLUE}For advanced options, use: ./scripts/dev-setup-and-run.sh --help${NC}"
echo ""

# Ask for confirmation unless --yes is passed
if [[ "$1" != "--yes" ]]; then
    read -p "Continue? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
fi

# Run the main setup script with Tauri mode (default)
exec "$(dirname "$0")/dev-setup-and-run.sh" -v 