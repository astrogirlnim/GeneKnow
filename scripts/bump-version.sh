#!/bin/bash

# 🔖 GeneKnow Version Bumping Script
# This script allows manual version bumping during development

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

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

# Function to show usage
show_usage() {
    echo "🔖 GeneKnow Version Bumping Script"
    echo ""
    echo "Usage: $0 [patch|minor|major]"
    echo ""
    echo "Options:"
    echo "  patch   - Bump patch version (0.1.0 → 0.1.1)"
    echo "  minor   - Bump minor version (0.1.0 → 0.2.0)"
    echo "  major   - Bump major version (0.1.0 → 1.0.0)"
    echo ""
    echo "If no option is provided, defaults to 'patch'"
    echo ""
    echo "Examples:"
    echo "  $0 patch   # 0.1.0 → 0.1.1"
    echo "  $0 minor   # 0.1.0 → 0.2.0"
    echo "  $0 major   # 0.1.0 → 1.0.0"
    echo ""
}

# Check if we're in the right directory
if [[ ! -f "$PROJECT_ROOT/desktop/src-tauri/tauri.conf.json" ]]; then
    print_error "This script must be run from the project root or scripts directory"
    exit 1
fi

# Get version bump type
VERSION_TYPE="${1:-patch}"

# Validate version type
if [[ ! "$VERSION_TYPE" =~ ^(patch|minor|major)$ ]]; then
    print_error "Invalid version type: $VERSION_TYPE"
    show_usage
    exit 1
fi

print_info "Starting version bump: $VERSION_TYPE"

# Get current version from tauri.conf.json
TAURI_CONFIG="$PROJECT_ROOT/desktop/src-tauri/tauri.conf.json"
PACKAGE_JSON="$PROJECT_ROOT/desktop/ui/package.json"

if [[ ! -f "$TAURI_CONFIG" ]]; then
    print_error "tauri.conf.json not found at $TAURI_CONFIG"
    exit 1
fi

if [[ ! -f "$PACKAGE_JSON" ]]; then
    print_error "package.json not found at $PACKAGE_JSON"
    exit 1
fi

# Check if jq is available
if ! command -v jq &> /dev/null; then
    print_error "jq is required but not installed. Install with: brew install jq"
    exit 1
fi

# Get current versions
CURRENT_VERSION=$(cat "$TAURI_CONFIG" | jq -r '.version')
PACKAGE_VERSION=$(cat "$PACKAGE_JSON" | jq -r '.version')

print_info "Current tauri.conf.json version: $CURRENT_VERSION"
print_info "Current package.json version: $PACKAGE_VERSION"

# Verify consistency
if [[ "$CURRENT_VERSION" != "$PACKAGE_VERSION" ]]; then
    print_warning "Version mismatch detected!"
    print_warning "tauri.conf.json: $CURRENT_VERSION"
    print_warning "package.json: $PACKAGE_VERSION"
    print_warning "Using tauri.conf.json version: $CURRENT_VERSION"
fi

# Parse current version (format: x.y.z)
IFS='.' read -r MAJOR MINOR PATCH <<< "$CURRENT_VERSION"

print_info "Current version components: MAJOR=$MAJOR, MINOR=$MINOR, PATCH=$PATCH"

# Calculate new version based on type
case "$VERSION_TYPE" in
    "major")
        NEW_MAJOR=$((MAJOR + 1))
        NEW_MINOR=0
        NEW_PATCH=0
        print_info "Major version bump: $MAJOR.$MINOR.$PATCH → $NEW_MAJOR.$NEW_MINOR.$NEW_PATCH"
        ;;
    "minor")
        NEW_MAJOR=$MAJOR
        NEW_MINOR=$((MINOR + 1))
        NEW_PATCH=0
        print_info "Minor version bump: $MAJOR.$MINOR.$PATCH → $NEW_MAJOR.$NEW_MINOR.$NEW_PATCH"
        ;;
    "patch"|*)
        NEW_MAJOR=$MAJOR
        NEW_MINOR=$MINOR
        NEW_PATCH=$((PATCH + 1))
        print_info "Patch version bump: $MAJOR.$MINOR.$PATCH → $NEW_MAJOR.$NEW_MINOR.$NEW_PATCH"
        ;;
esac

NEW_VERSION="$NEW_MAJOR.$NEW_MINOR.$NEW_PATCH"
NEW_TAG="v$NEW_VERSION"

print_info "New version: $NEW_VERSION"
print_info "New tag: $NEW_TAG"

# Check if tag already exists
if git tag -l | grep -q "^$NEW_TAG$"; then
    print_error "Tag $NEW_TAG already exists!"
    print_error "Please choose a different version or delete the existing tag first"
    exit 1
fi

# Update tauri.conf.json
print_info "Updating tauri.conf.json..."
jq --arg version "$NEW_VERSION" '.version = $version' "$TAURI_CONFIG" > "$TAURI_CONFIG.tmp"
mv "$TAURI_CONFIG.tmp" "$TAURI_CONFIG"
print_success "Updated tauri.conf.json to version $NEW_VERSION"

# Update package.json
print_info "Updating package.json..."
jq --arg version "$NEW_VERSION" '.version = $version' "$PACKAGE_JSON" > "$PACKAGE_JSON.tmp"
mv "$PACKAGE_JSON.tmp" "$PACKAGE_JSON"
print_success "Updated package.json to version $NEW_VERSION"

# Commit changes
print_info "Committing version changes..."
git add "$TAURI_CONFIG" "$PACKAGE_JSON"

if git diff --cached --quiet; then
    print_warning "No changes to commit"
else
    git commit -m "🔖 Bump version to $NEW_VERSION"
    print_success "Version files committed"
fi

# Create git tag
print_info "Creating git tag..."
git tag "$NEW_TAG"
print_success "Git tag created: $NEW_TAG"

print_success "Version bump completed successfully!"
print_info "Summary:"
print_info "  Previous version: $CURRENT_VERSION"
print_info "  New version: $NEW_VERSION"
print_info "  Tag: $NEW_TAG"
print_info "  Type: $VERSION_TYPE"
print_info ""
print_info "To push changes and tag:"
print_info "  git push origin main"
print_info "  git push origin $NEW_TAG"
print_info ""
print_info "To trigger a release, push to main or use GitHub Actions workflow_dispatch" 