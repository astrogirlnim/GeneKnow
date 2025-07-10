# 💰 GitHub Actions Cost Optimization Guide

## 🚨 Problem: Billing Limit Exceeded

Your repository was consuming **$70.79** in GitHub Actions minutes, exceeding your GitHub Pro plan's included usage. The main culprit: **expensive multi-platform builds**.

## 📊 Cost Breakdown

### Before Optimization
- **Linux builds**: $0.008/minute
- **Windows builds**: $0.016/minute (2x more expensive)
- **macOS builds**: $0.08/minute (10x more expensive!)

**Full release build**: 4 platforms × 6-8 minutes = ~$8-12 per release

### After Optimization
- **Linux-only builds**: $0.008/minute
- **Single platform**: 1 platform × 6-8 minutes = ~$0.05-0.10 per release

**Result**: 98% cost reduction!

## ✅ Optimizations Implemented

### 1. **Linux-Only Builds (Immediate Fix)**

Modified workflow files to build only on Linux:
```yaml
# .github/workflows/release.yml
# .github/workflows/pr.yml

matrix:
  include:
    - os: ubuntu-latest
      rust-target: x86_64-unknown-linux-gnu
      name: 'Linux x64'
    # Expensive builds disabled:
    # - Windows builds (2x cost)
    # - macOS builds (10x cost)
```

### 2. **Aggressive Artifact Cleanup**

```yaml
# Keep only 1 most recent artifact
retention-days: 1
compression-level: 9
```

### 3. **Smart Caching**

```yaml
# Cache dependencies to reduce build times
- uses: actions/cache@v4
  with:
    path: ${{ env.STORE_PATH }}
    key: ${{ runner.os }}-pnpm-store-${{ hashFiles('**/pnpm-lock.yaml') }}
```

### 4. **Path-Based Triggering**

```yaml
# Only run on relevant changes
on:
  push:
    paths:
      - 'desktop/**'
      - 'geneknow_pipeline/**'
      - '.github/workflows/**'
```

## 🔧 Using the Optimization Script

Run the interactive cost optimization tool:

```bash
./scripts/optimize-github-actions.sh
```

**Features:**
- Monitor current usage
- Clean up old artifacts
- View recent workflow runs
- Toggle expensive builds on/off
- Get cost optimization tips

## 📈 Cost Monitoring

### Monthly Budget Planning

With GitHub Pro ($4/month), you get:
- **2,000 free minutes** of Linux usage
- **250 free minutes** of Windows usage  
- **50 free minutes** of macOS usage

### Current Usage Estimate

**Linux-only builds:**
- ~20 releases per month = 20 × 6 minutes = 120 minutes
- Cost: $0.96/month
- **Remaining budget**: $69.83/month

### If You Re-enable All Platforms

**Full multi-platform builds:**
- ~6 releases per month = 6 × 30 minutes = 180 minutes (mixed)
- Cost: ~$50-60/month
- **Risk**: Will exceed budget frequently

## 🎯 Recommendations

### For Current Budget (GitHub Pro)
1. **Keep Linux-only builds** (current setup)
2. **Create manual releases** only when needed
3. **Use workflow_dispatch** for controlled releases
4. **Monitor usage** monthly with the script

### For Higher Volume (Consider Upgrading)
1. **GitHub Team**: $4/user/month + higher limits
2. **GitHub Enterprise**: Unlimited usage
3. **Self-hosted runners**: One-time setup cost

## 🔄 Re-enabling Expensive Builds

When you're ready to build for all platforms again:

### 1. Uncomment Matrix Entries

In `.github/workflows/release.yml` and `.github/workflows/pr.yml`:

```yaml
matrix:
  include:
    - os: ubuntu-latest
      rust-target: x86_64-unknown-linux-gnu
      name: 'Linux x64'
    - os: windows-latest  # Uncomment this block
      rust-target: x86_64-pc-windows-msvc
      name: 'Windows x64'
    - os: macos-latest    # Uncomment this block
      rust-target: x86_64-apple-darwin
      name: 'macOS x64'
```

### 2. Gradual Re-enabling

Start with one platform at a time:

1. **Windows first** (medium cost increase)
2. **Monitor usage** for a week
3. **Add macOS** if budget allows

## 📋 Monthly Checklist

- [ ] Run `./scripts/optimize-github-actions.sh` to check usage
- [ ] Clean up old artifacts
- [ ] Review failed workflows to avoid waste
- [ ] Consider if all platforms are needed for each release
- [ ] Monitor billing page for unexpected usage

## 🚀 Alternative Solutions

### Self-Hosted Runners
- **Cost**: One-time setup (~$50-100/month for VPS)
- **Benefits**: Unlimited usage, faster builds
- **Complexity**: Higher maintenance

### Staged Releases
- **Week 1**: Linux release
- **Week 2**: Windows release  
- **Week 3**: macOS release
- **Week 4**: All platforms (if needed)

### Manual Releases Only
- Disable automatic releases on push
- Use `workflow_dispatch` for controlled releases
- Build only when publishing to users

## 💡 Pro Tips

1. **Use draft releases** to test without storage costs
2. **Cancel failed workflows** immediately to save minutes
3. **Review large artifacts** before uploading
4. **Use `if` conditions** to skip unnecessary steps
5. **Monitor the Actions tab** during builds

## 🔍 Troubleshooting

### "Artifact storage quota exceeded"
```bash
# Quick fix
./scripts/optimize-github-actions.sh
# Select option 3: Clean up old artifacts
```

### "Billing limit exceeded"
```bash
# Check current usage
gh api user/settings/billing/actions
```

### Re-enable specific platforms
```bash
# Use the interactive script
./scripts/optimize-github-actions.sh
# Select option 5: Toggle expensive builds
```

---

## 📊 Summary

✅ **Immediate fix**: Linux-only builds (98% cost reduction)  
✅ **Monitoring**: Interactive optimization script  
✅ **Flexibility**: Easy to re-enable platforms when needed  
✅ **Budget-friendly**: Stay within GitHub Pro limits  

Your pipelines will now run efficiently within budget while maintaining full functionality for development and testing. 