# 🚀 Release Pipeline Trigger Guide

## Why Your Release Pipeline Might Not Be Triggering

### Current Situation
You're on the `deployment-refactor` branch, but the release pipeline only triggers on pushes to the `main` branch.

### Release Pipeline Triggers

The release pipeline (`release.yml`) triggers when:

1. **Push to `main` branch** with changes in these paths:
   - `desktop/**`
   - `geneknow_pipeline/**`
   - `langgraph/**`
   - `.github/workflows/release.yml`
   - `package.json`
   - `README.md`
   - `*.py`
   - `requirements.txt`
   - `requirements-lite.txt`
   - `docs/**`

2. **Manual workflow dispatch** via GitHub Actions UI

## How to Trigger a Release

### Option 1: Merge to Main (Recommended)

```bash
# 1. Create a PR from your branch
git push origin deployment-refactor

# 2. Go to GitHub and create a Pull Request
# 3. Review and merge the PR
# 4. The release pipeline will automatically trigger on merge
```

### Option 2: Manual Trigger

1. Go to your repository on GitHub
2. Click on "Actions" tab
3. Find "🚀 Release Pipeline" in the left sidebar
4. Click "Run workflow"
5. Select:
   - Branch: `main`
   - Version type: `patch`, `minor`, or `major`
6. Click "Run workflow"

### Option 3: Direct Push to Main (Not Recommended)

```bash
# Only if you have direct push permissions
git checkout main
git merge deployment-refactor
git push origin main
```

## Verifying the Pipeline Will Trigger

Before merging, ensure you have changes in the monitored paths:

```bash
# Check what files have changed
git diff --name-only main..deployment-refactor

# Look for files in:
# - desktop/
# - geneknow_pipeline/
# - .github/workflows/
# - docs/
```

## Common Issues and Solutions

### Issue 1: No Changes in Monitored Paths
If your changes are outside the monitored paths, the pipeline won't trigger.

**Solution**: Either:
- Add a small change to a monitored file (e.g., update `docs/`)
- Use manual workflow dispatch

### Issue 2: Branch Protection Rules
If `main` has protection rules, you must go through a PR.

**Solution**: Create and merge a PR as described in Option 1.

### Issue 3: Pipeline Failures
If the pipeline starts but fails:

1. **Validation errors**: Check ESLint and TypeScript
   ```bash
   cd desktop/ui
   pnpm lint
   pnpm exec tsc --noEmit
   ```

2. **Build errors**: Test locally first
   ```bash
   pnpm build
   pnpm run tauri-build
   ```

3. **Version conflicts**: The pipeline handles these automatically

## Testing Before Release

Always test locally before triggering a release:

```bash
# 1. Lint and type check
cd desktop/ui
pnpm lint
pnpm exec tsc --noEmit

# 2. Build frontend
pnpm build

# 3. Bundle Python
cd ../scripts
./bundle-python-optimized.sh

# 4. Check Rust
cd ../src-tauri
cargo check

# 5. Full build test (optional)
cd ../ui
pnpm run tauri-build
```

## Release Pipeline Features

When triggered, the pipeline will:

1. **Auto-increment version** (patch by default)
2. **Build for all platforms**:
   - Linux x64 (`.deb`, `.AppImage`)
   - Windows x64 (`.msi`, `.exe`)
   - macOS x64 (`.dmg`, `.app`)
   - macOS ARM64 (`.dmg`, `.app`)
3. **Create GitHub Release** with installers
4. **Clean up old artifacts** to save storage

## Best Practices

1. **Always create a PR** for visibility and review
2. **Test locally first** to avoid failed builds
3. **Use semantic commit messages** for better changelogs
4. **Monitor the Actions tab** during the release
5. **Download and test** the built artifacts

## Quick Checklist

- [ ] All tests pass locally
- [ ] No linting errors
- [ ] Frontend builds successfully
- [ ] Python bundling works
- [ ] Rust compiles without errors
- [ ] Changes are in monitored paths
- [ ] PR is ready to merge to `main` 