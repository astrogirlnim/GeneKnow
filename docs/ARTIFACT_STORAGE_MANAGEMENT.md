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

## 🚀 **GitHub Plan Upgrade Options**

### **Current Storage Limits by Plan:**

| Plan | Artifact Storage | Monthly Cost | Best For |
|------|------------------|--------------|----------|
| **Free** | 500MB | $0 | Small projects, personal use |
| **Pro** | 1GB | $4/month | Individual developers |
| **Team** | 2GB | $4/user/month | Small teams |
| **Enterprise** | 50GB | $21/user/month | Large organizations |

### **Upgrade Recommendations:**

**For GenePredict:**
- **GitHub Pro** ($4/month) - Doubles storage to 1GB, sufficient for most development
- **GitHub Team** ($4/user/month) - If you have collaborators
- **Enterprise** - Only if you plan massive scale with many releases

### **How to Upgrade:**
1. Go to your GitHub account → **Settings** → **Billing**
2. Choose **Change plan**
3. Select **GitHub Pro** for individual use
4. Billing is monthly, can cancel anytime

## 🔧 **Optimizations Applied to Your Workflows**

### **1. Selective Artifact Cleanup**
```yaml
# Smart cleanup strategy:
- Failed builds: Delete after 3 days
- PR builds: Delete after 1 day (testing only)
- Successful main builds: Keep for 14 days (until release)
```

### **2. Optimized Artifact Uploads**
**Before:** Uploaded entire bundle directory (~200MB per platform)
**After:** Only upload essential installer files (~50MB per platform)

```yaml
# Only upload installer files
path: |
  desktop/src-tauri/target/*/release/bundle/**/*.dmg
  desktop/src-tauri/target/*/release/bundle/**/*.msi
  desktop/src-tauri/target/*/release/bundle/**/*.deb
  desktop/src-tauri/target/*/release/bundle/**/*.AppImage
  desktop/src-tauri/target/*/release/bundle/**/*.exe
```

### **3. Shorter Retention Periods**
- **Artifacts**: 7 days (was 90 days)
- **Compression**: Maximum level 9
- **PR builds**: No artifacts uploaded (validation only)

### **4. Release vs Artifact Protection**
**Important:** Our cleanup **ONLY** affects temporary build artifacts, **NOT** permanent GitHub releases.

| Type | What It Is | Retention | Protected |
|------|------------|-----------|-----------|
| **GitHub Actions Artifacts** | Temporary build files | 1-14 days | ❌ (cleaned up) |
| **GitHub Releases** | Permanent download links | Forever | ✅ (protected) |

## 📊 **Storage Usage Estimation**

### **Before Optimization:**
```
Per release: 3 platforms × 200MB = 600MB
With 90-day retention: 600MB × 10 releases = 6GB ❌
```

### **After Optimization:**
```
Per release: 3 platforms × 50MB = 150MB
With 7-day retention: 150MB × 2 releases = 300MB ✅
```

### **With GitHub Pro (1GB):**
```
Headroom: 1GB - 300MB = 700MB available
Can handle: ~13 releases simultaneously
```

## 🔒 **Release Protection Strategy**

### **What Gets Deleted:**
- ❌ Failed PR builds (after 3 days)
- ❌ Successful PR builds (after 1 day)
- ❌ Old successful main builds (after 14 days)

### **What Stays Protected:**
- ✅ Recent successful main builds (14 days)
- ✅ GitHub releases (forever)
- ✅ Release assets (forever)

### **Release Creation Process:**
1. **Build** → Creates temporary artifact
2. **Test** → Validates build quality
3. **Release** → Moves artifacts to permanent release
4. **Cleanup** → Deletes temporary artifacts (not release!)

## 💰 **Cost-Benefit Analysis**

### **Free Plan + Optimizations:**
- **Cost**: $0/month
- **Storage**: 500MB
- **Releases**: ~3 simultaneous releases
- **Best for**: Solo development, infrequent releases

### **Pro Plan + Optimizations:**
- **Cost**: $4/month ($48/year)
- **Storage**: 1GB
- **Releases**: ~6 simultaneous releases
- **Best for**: Active development, regular releases

### **Recommendation:**
Start with **optimized Free plan** → upgrade to **Pro** when you hit limits consistently. 