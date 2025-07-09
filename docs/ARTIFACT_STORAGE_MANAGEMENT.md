# 📦 GitHub Actions Artifact Storage Management Guide

## 🚨 Problem: Artifact Storage Quota Exceeded

When you see `Failed to CreateArtifact: Artifact storage quota has been hit`, it means your repository has exceeded GitHub's artifact storage limits.

## 📊 GitHub Storage Limits

| Plan | Storage Limit | Retention |
|------|---------------|-----------|
| **Free** | 500MB | 90 days |
| **Pro** | 1GB | 90 days |
| **Team** | 2GB | 90 days |
| **Enterprise** | 50GB | 90 days |

## 🔍 Check Current Storage Usage

### Via GitHub Web UI
1. Go to your repo → **Actions** → **Management** → **Storage**
2. View current usage and breakdown by workflow

### Via GitHub CLI
```bash
# Install GitHub CLI if not installed
brew install gh  # macOS
# or: choco install gh  # Windows
# or: sudo apt install gh  # Linux

# Check current storage usage
gh api repos/:owner/:repo/actions/cache/usage

# List all artifacts
gh api repos/:owner/:repo/actions/artifacts --jq '.artifacts[] | {name: .name, size_in_bytes: .size_in_bytes, created_at: .created_at}'
```

## 🧹 Immediate Solutions

### 1. **Manual Cleanup (Quickest)**
```bash
# Delete ALL artifacts older than 7 days
gh api repos/:owner/:repo/actions/artifacts --jq '.artifacts[] | select(.created_at < (now - 7*24*3600 | todate)) | .id' | xargs -I {} gh api -X DELETE repos/:owner/:repo/actions/artifacts/{}

# Delete specific artifact types
gh api repos/:owner/:repo/actions/artifacts --jq '.artifacts[] | select(.name | contains("Build")) | .id' | xargs -I {} gh api -X DELETE repos/:owner/:repo/actions/artifacts/{}
```

### 2. **Automated Cleanup (Implemented)**
We've added an automated cleanup job to your release pipeline:
- ✅ **Runs before every build** to free up space
- ✅ **Deletes artifacts older than 7 days**
- ✅ **Prevents future quota issues**

### 3. **Optimized Artifact Settings**
Updated your release workflow with:
- **Retention**: 7 days (was 90 days)
- **Compression**: Level 9 (maximum)
- **Cleanup**: Automatic before builds

## 🛠️ Workflow Optimizations Applied

### Release Pipeline Changes:
```yaml
# Added automatic cleanup job
cleanup:
  name: 🧹 Cleanup Old Artifacts
  runs-on: ubuntu-latest
  # Deletes artifacts older than 7 days before building

# Optimized artifact uploads
upload-artifact@v4:
  retention-days: 7        # Reduced from 90 days
  compression-level: 9     # Maximum compression
```

### PR Pipeline (Already Optimized):
- ✅ **No artifact uploads** - only builds for validation
- ✅ **Saves storage space** - builds are tested but not stored
- ✅ **Faster execution** - no upload time

## 📈 Storage Monitoring

### Set up notifications for storage usage:
```bash
# Create a monitoring script
cat > check_storage.sh << 'EOF'
#!/bin/bash
USAGE=$(gh api repos/:owner/:repo/actions/cache/usage --jq '.total_active_cache_size_in_bytes')
LIMIT=$((500 * 1024 * 1024))  # 500MB for free plan

if [ $USAGE -gt $((LIMIT * 80 / 100)) ]; then
    echo "⚠️  Storage usage is at $(($USAGE * 100 / $LIMIT))% of limit"
    echo "🧹 Consider running cleanup..."
fi
EOF

chmod +x check_storage.sh
```

## 🔧 Best Practices

### 1. **Minimize Artifact Sizes**
- Only upload essential files (not entire target directories)
- Use maximum compression
- Filter out unnecessary files

### 2. **Short Retention Periods**
- **PR builds**: No artifacts (validation only)
- **Release builds**: 7 days (enough for debugging)
- **Nightly builds**: 3 days max

### 3. **Selective Uploads**
```yaml
# Only upload specific files
path: |
  desktop/src-tauri/target/*/release/bundle/*.dmg
  desktop/src-tauri/target/*/release/bundle/*.msi
  desktop/src-tauri/target/*/release/bundle/*.deb
```

### 4. **Conditional Uploads**
```yaml
# Only upload on main branch
if: github.ref == 'refs/heads/main'
```

## 🚀 Next Steps

1. **Push the updated workflows** to trigger cleanup
2. **Monitor storage usage** regularly
3. **Consider upgrading plan** if needed for larger projects
4. **Review retention policies** periodically

## 📊 Storage Calculator

Estimate your storage needs:
```bash
# Calculate build sizes
find desktop/src-tauri/target -name "*.dmg" -o -name "*.msi" -o -name "*.deb" | xargs du -h

# Estimate monthly usage
# (Build size × Platforms × Releases per month × Retention days) / 30
```

## 🆘 Emergency Cleanup

If you hit quota and need immediate relief:
```bash
# Nuclear option: Delete ALL artifacts
gh api repos/:owner/:repo/actions/artifacts --jq '.artifacts[].id' | xargs -I {} gh api -X DELETE repos/:owner/:repo/actions/artifacts/{}

# Verify cleanup
gh api repos/:owner/:repo/actions/artifacts --jq '.total_count'
```

---

💡 **Pro Tip**: Set up a scheduled workflow to run cleanup weekly to prevent future issues. 