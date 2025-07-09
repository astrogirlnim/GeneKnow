#!/bin/bash

# 🔧 Git Hooks Installation Script for GenePredict
# This script installs git hooks to prevent versioning conflicts

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
HOOKS_DIR="$PROJECT_ROOT/.git/hooks"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_info "Installing GenePredict git hooks..."

# Check if we're in a git repository
if [[ ! -d "$PROJECT_ROOT/.git" ]]; then
    print_error "This doesn't appear to be a git repository"
    exit 1
fi

# Create hooks directory if it doesn't exist
mkdir -p "$HOOKS_DIR"

# Install pre-push hook
print_info "Installing pre-push hook..."
cat > "$HOOKS_DIR/pre-push" << 'EOF'
#!/bin/bash

# 🔧 GenePredict Pre-Push Hook
# Prevents manual version tag conflicts

# Colors
RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m'

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

# Check if we're pushing version tags
while read local_ref local_sha remote_ref remote_sha; do
    if [[ "$remote_ref" =~ refs/tags/v[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
        TAG_NAME=$(basename "$remote_ref")
        
        print_warning "Detected version tag push: $TAG_NAME"
        print_warning "Manual version tags can cause conflicts with the release pipeline!"
        print_warning ""
        print_warning "Recommended workflow:"
        print_warning "1. Use the automated release pipeline (push to main)"
        print_warning "2. Use ./scripts/bump-version.sh for development"
        print_warning "3. See docs/VERSIONING_GUIDELINES.md for details"
        print_warning ""
        
        # Ask for confirmation
        echo -n "Are you sure you want to push this version tag? (y/N): "
        read -r response
        
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            print_error "Tag push cancelled. Use the automated release pipeline instead."
            exit 1
        fi
        
        print_success "Proceeding with manual tag push..."
    fi
done

# Check for version file consistency before push
TAURI_VERSION=""
PACKAGE_VERSION=""

if [[ -f "desktop/src-tauri/tauri.conf.json" ]]; then
    TAURI_VERSION=$(cat desktop/src-tauri/tauri.conf.json | jq -r '.version' 2>/dev/null || echo "")
fi

if [[ -f "desktop/ui/package.json" ]]; then
    PACKAGE_VERSION=$(cat desktop/ui/package.json | jq -r '.version' 2>/dev/null || echo "")
fi

if [[ -n "$TAURI_VERSION" && -n "$PACKAGE_VERSION" && "$TAURI_VERSION" != "$PACKAGE_VERSION" ]]; then
    print_error "Version mismatch detected!"
    print_error "tauri.conf.json: $TAURI_VERSION"
    print_error "package.json: $PACKAGE_VERSION"
    print_error ""
    print_error "Fix with: ./scripts/bump-version.sh patch"
    exit 1
fi

exit 0
EOF

# Make the hook executable
chmod +x "$HOOKS_DIR/pre-push"
print_success "Pre-push hook installed"

# Install prepare-commit-msg hook
print_info "Installing prepare-commit-msg hook..."
cat > "$HOOKS_DIR/prepare-commit-msg" << 'EOF'
#!/bin/bash

# 🔧 GenePredict Prepare-Commit-Msg Hook
# Adds helpful context to commit messages

COMMIT_MSG_FILE=$1
COMMIT_SOURCE=$2

# Only modify commit message if it's a regular commit (not merge, squash, etc.)
if [ -z "$COMMIT_SOURCE" ]; then
    # Check if this is a version bump commit
    if grep -q "package.json\|tauri.conf.json" <<< "$(git diff --cached --name-only)"; then
        # Get the current version from tauri.conf.json
        if [[ -f "desktop/src-tauri/tauri.conf.json" ]]; then
            VERSION=$(cat desktop/src-tauri/tauri.conf.json | jq -r '.version' 2>/dev/null || echo "")
            if [[ -n "$VERSION" && ! $(head -1 "$COMMIT_MSG_FILE") =~ "🔖 Bump version" ]]; then
                # Add version context if it's not already there
                echo "" >> "$COMMIT_MSG_FILE"
                echo "# Version: $VERSION" >> "$COMMIT_MSG_FILE"
                echo "# ℹ️  This appears to be a version-related commit" >> "$COMMIT_MSG_FILE"
                echo "# ℹ️  Use ./scripts/bump-version.sh for consistent versioning" >> "$COMMIT_MSG_FILE"
                echo "# ℹ️  See docs/VERSIONING_GUIDELINES.md for details" >> "$COMMIT_MSG_FILE"
            fi
        fi
    fi
fi
EOF

# Make the hook executable
chmod +x "$HOOKS_DIR/prepare-commit-msg"
print_success "Prepare-commit-msg hook installed"

print_success "Git hooks installation completed!"
print_info ""
print_info "Installed hooks:"
print_info "  • pre-push: Prevents manual version tag conflicts"
print_info "  • prepare-commit-msg: Adds version context to commits"
print_info ""
print_info "To uninstall hooks:"
print_info "  rm $HOOKS_DIR/pre-push"
print_info "  rm $HOOKS_DIR/prepare-commit-msg"
print_info ""
print_info "For more information, see docs/VERSIONING_GUIDELINES.md" 