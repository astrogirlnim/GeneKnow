# 🔖 GeneKnow Versioning Guidelines

## Overview

This document outlines the versioning strategy for GeneKnow to prevent tag conflicts and ensure smooth releases.

## 🎯 Core Principles

1. **Single Source of Truth**: `desktop/src-tauri/tauri.conf.json` is the primary version source
2. **Automated Versioning**: Let the release pipeline handle version bumping
3. **No Manual Tags**: Avoid creating version tags manually
4. **Consistent Formatting**: Always use semantic versioning (e.g., `v1.2.3`)

## 🔄 Versioning Workflow

### For Regular Development

1. **Make your changes** in a feature branch
2. **Create a pull request** to `main`
3. **Merge to main** - the release pipeline will automatically:
   - Increment the patch version
   - Update `tauri.conf.json` and `package.json`
   - Create the appropriate git tag
   - Build and release

### For Manual Version Bumping (Development Only)

If you need to manually bump versions during development:

```bash
# For patch version (0.1.0 → 0.1.1)
./scripts/bump-version.sh patch

# For minor version (0.1.0 → 0.2.0)
./scripts/bump-version.sh minor

# For major version (0.1.0 → 1.0.0)
./scripts/bump-version.sh major
```

**⚠️ Important**: Only use manual bumping for development. Production releases should use the automated pipeline.

## 🚫 What NOT to Do

### ❌ Never Create Tags Manually
```bash
# DON'T DO THIS - causes conflicts
git tag v0.1.5
git push origin v0.1.5
```

### ❌ Never Edit Version Files Directly
```bash
# DON'T DO THIS - causes inconsistencies
vim desktop/src-tauri/tauri.conf.json
```

### ❌ Never Force Push Tags
```bash
# DON'T DO THIS - breaks remote state
git push --force origin v0.1.5
```

## 🛠️ Troubleshooting Common Issues

### Issue: Tag Conflicts During Pull
```bash
# Error: ! [rejected] v0.1.2 -> v0.1.2 (would clobber existing tag)
```

**Solution**:
```bash
# Delete conflicting local tags
git tag -d v0.1.2 v0.1.3

# Pull with tags
git pull --tags origin main
```

### Issue: Version Mismatch Between Files
```bash
# Check current versions
cat desktop/src-tauri/tauri.conf.json | jq '.version'
cat desktop/ui/package.json | jq '.version'

# Fix with bump script
./scripts/bump-version.sh patch
```

### Issue: Release Pipeline Fails
1. Check if tags already exist
2. Verify version consistency
3. Check GitHub Actions logs
4. Contact team lead if needed

## 📋 Release Pipeline Details

### Automatic Triggers
- Push to `main` branch
- Changes to `desktop/**`, `*.py`, `requirements.txt`, or docs

### Manual Triggers
- GitHub Actions → "🚀 Release Pipeline" → "Run workflow"
- Choose version type: `patch`, `minor`, or `major`

### Pipeline Steps
1. **Version Calculation** - Determines next version
2. **Tag Management** - Creates/updates git tags
3. **Validation** - Runs tests and linting
4. **Build** - Creates cross-platform binaries
5. **Release** - Publishes to GitHub Releases

## 🔧 Developer Guidelines

### Before Starting Work
```bash
# Always start with latest main
git checkout main
git pull origin main

# Create feature branch
git checkout -b feature/your-feature
```

### Before Creating PR
```bash
# Check if versions are consistent
./scripts/bump-version.sh --check  # (if script supports it)

# Or manually check
cat desktop/src-tauri/tauri.conf.json | jq '.version'
cat desktop/ui/package.json | jq '.version'
```

### After Merge to Main
- Monitor GitHub Actions for successful release
- Check that new version appears in releases
- Verify download links work

## 📞 Getting Help

### Quick Fixes
1. **Tag conflicts**: Delete local tags and pull
2. **Version mismatch**: Run `./scripts/bump-version.sh patch`
3. **Build failures**: Check GitHub Actions logs

### When to Contact Team
- Release pipeline repeatedly fails
- Version numbers are completely wrong
- Need emergency hotfix release
- Major version changes

## 📊 Version Strategy

### Current Versioning Scheme
- **Major** (1.0.0): Breaking changes, major new features
- **Minor** (0.1.0): New features, backwards compatible
- **Patch** (0.0.1): Bug fixes, small improvements

### Release Cadence
- **Patch**: As needed (bug fixes)
- **Minor**: Weekly/bi-weekly (feature releases)
- **Major**: Monthly/quarterly (major milestones)

---

## 🔄 Quick Reference

### Most Common Commands
```bash
# Check current version
cat desktop/src-tauri/tauri.conf.json | jq '.version'

# Fix tag conflicts
git tag -d v0.1.2 v0.1.3
git pull --tags origin main

# Manual version bump (development only)
./scripts/bump-version.sh patch

# Check git tags
git tag -l

# Check remote tags
git ls-remote --tags origin
```

### Key Files
- `desktop/src-tauri/tauri.conf.json` - Primary version source
- `desktop/ui/package.json` - Frontend version (must match)
- `.github/workflows/release.yml` - Release pipeline
- `scripts/bump-version.sh` - Manual version bumping

---

**Remember**: When in doubt, let the automated pipeline handle versioning. Manual intervention should be rare and carefully coordinated with the team. 