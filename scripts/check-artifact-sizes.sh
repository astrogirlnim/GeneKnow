#!/bin/bash

# 📊 GitHub Actions Artifact Size Analysis
# Analyzes current artifact usage and provides storage recommendations

set -e

echo "📊 GitHub Actions Artifact Size Analysis"
echo "========================================"
echo ""

# Check if gh CLI is installed
if ! command -v gh &> /dev/null; then
    echo "❌ GitHub CLI (gh) not found. Please install it first:"
    echo "   brew install gh  # macOS"
    echo "   choco install gh # Windows"
    echo "   sudo apt install gh # Linux"
    exit 1
fi

# Get repository information
REPO_INFO=$(gh repo view --json nameWithOwner -q .nameWithOwner 2>/dev/null || echo "")
if [ -z "$REPO_INFO" ]; then
    echo "❌ Could not determine repository. Please run from a git repository."
    exit 1
fi

echo "📁 Repository: $REPO_INFO"
echo ""

# Function to format bytes
format_bytes() {
    local bytes=$1
    if [ $bytes -lt 1024 ]; then
        echo "${bytes}B"
    elif [ $bytes -lt 1048576 ]; then
        echo "$((bytes / 1024))KB"
    elif [ $bytes -lt 1073741824 ]; then
        echo "$((bytes / 1048576))MB"
    else
        echo "$((bytes / 1073741824))GB"
    fi
}

# Get current storage usage
echo "📊 Current Storage Usage:"
echo "========================"
USAGE_JSON=$(gh api repos/$REPO_INFO/actions/cache/usage 2>/dev/null || echo "{}")
if [ "$USAGE_JSON" != "{}" ]; then
    TOTAL_BYTES=$(echo "$USAGE_JSON" | jq -r '.total_active_caches_size_in_bytes // 0')
    TOTAL_COUNT=$(echo "$USAGE_JSON" | jq -r '.total_active_caches_count // 0')
    echo "💾 Cache Storage: $(format_bytes $TOTAL_BYTES) across $TOTAL_COUNT caches"
fi

# Get artifact information
echo ""
echo "📦 Artifact Breakdown:"
echo "====================="
ARTIFACTS_JSON=$(gh api repos/$REPO_INFO/actions/artifacts --paginate 2>/dev/null || echo '{"artifacts": []}')
TOTAL_ARTIFACTS=$(echo "$ARTIFACTS_JSON" | jq -r '.artifacts | length')
TOTAL_SIZE=$(echo "$ARTIFACTS_JSON" | jq -r '[.artifacts[].size_in_bytes] | add // 0')

echo "Total Artifacts: $TOTAL_ARTIFACTS"
echo "Total Size: $(format_bytes $TOTAL_SIZE)"
echo ""

# Show largest artifacts
echo "🏆 Largest Artifacts (Top 10):"
echo "$ARTIFACTS_JSON" | jq -r '.artifacts | sort_by(.size_in_bytes) | reverse | .[0:10] | .[] | "\(.name): \(.size_in_bytes / 1048576 | floor)MB (\(.created_at | split("T")[0]))"'
echo ""

# Show breakdown by platform
echo "📊 Platform Breakdown:"
echo "$ARTIFACTS_JSON" | jq -r '.artifacts | group_by(.name | split(" ")[0]) | .[] | {platform: .[0].name | split(" ")[0], count: length, size: ([.[].size_in_bytes] | add), latest: (.[0].created_at | split("T")[0])} | "\(.platform): \(.count) artifacts, \(.size / 1048576 | floor)MB (latest: \(.latest))"' 2>/dev/null || echo "No platform data available"
echo ""

# Calculate GitHub plan recommendations
echo "💰 GitHub Plan Analysis:"
echo "========================"
CURRENT_LIMIT=$((500 * 1024 * 1024))  # 500MB free plan
PRO_LIMIT=$((1024 * 1024 * 1024))     # 1GB pro plan

USAGE_PERCENT=$((TOTAL_SIZE * 100 / CURRENT_LIMIT))
echo "Current usage: $(format_bytes $TOTAL_SIZE) / 500MB (${USAGE_PERCENT}%)"

if [ $TOTAL_SIZE -gt $CURRENT_LIMIT ]; then
    OVERAGE=$((TOTAL_SIZE - CURRENT_LIMIT))
    echo "❌ OVER QUOTA by $(format_bytes $OVERAGE)"
    echo ""
    echo "🚀 RECOMMENDATION: Upgrade to GitHub Pro"
    echo "   - Cost: \$4/month (\$48/year)"
    echo "   - Storage: 1GB (doubles current limit)"
    echo "   - Your usage: $(format_bytes $TOTAL_SIZE) / 1GB ($(($TOTAL_SIZE * 100 / PRO_LIMIT))%)"
    echo ""
elif [ $TOTAL_SIZE -gt $((CURRENT_LIMIT * 80 / 100)) ]; then
    echo "⚠️  High usage - consider upgrading soon"
    echo ""
    echo "💡 RECOMMENDATION: Consider GitHub Pro"
    echo "   - Cost: \$4/month (\$48/year)"
    echo "   - Storage: 1GB (doubles current limit)"
    echo "   - Headroom: $(format_bytes $((PRO_LIMIT - TOTAL_SIZE))) additional space"
    echo ""
else
    echo "✅ Usage within limits"
    echo ""
fi

# Provide optimization suggestions
echo "🔧 Optimization Suggestions:"
echo "============================"
echo ""
echo "1. **Immediate Relief:**"
echo "   - Run: ./scripts/emergency-cleanup.sh"
echo "   - Select option 1 (delete all) or 2 (keep 1 recent)"
echo ""
echo "2. **Short-term (Free plan):**"
echo "   - Ultra-aggressive cleanup (already implemented)"
echo "   - 1-day retention (already implemented)"
echo "   - Maximum compression (already implemented)"
echo ""
echo "3. **Medium-term (Recommended):**"
echo "   - **Upgrade to GitHub Pro (\$4/month)**"
echo "   - Doubles storage to 1GB"
echo "   - No workflow changes needed"
echo ""
echo "4. **Long-term optimizations:**"
echo "   - Conditional uploads (only on releases)"
echo "   - Split platform builds"
echo "   - External storage solutions"
echo ""

# Show next steps
echo "🎯 Next Steps:"
echo "=============="
echo ""
if [ $TOTAL_SIZE -gt $CURRENT_LIMIT ]; then
    echo "🚨 URGENT: Over quota - immediate action needed"
    echo "1. Run emergency cleanup: ./scripts/emergency-cleanup.sh"
    echo "2. Upgrade to GitHub Pro: GitHub Settings → Billing → Change Plan"
    echo "3. Continue development without storage worries"
else
    echo "✅ You're within limits, but consider:"
    echo "1. Monitor usage regularly"
    echo "2. Consider GitHub Pro for peace of mind"
    echo "3. Set up automated cleanup"
fi
echo ""

echo "📋 Summary:"
echo "==========="
echo "Repository: $REPO_INFO"
echo "Artifacts: $TOTAL_ARTIFACTS"
echo "Storage: $(format_bytes $TOTAL_SIZE) / 500MB (${USAGE_PERCENT}%)"
echo "Status: $([ $TOTAL_SIZE -gt $CURRENT_LIMIT ] && echo "❌ OVER QUOTA" || echo "✅ Within limits")"
echo "Recommendation: $([ $TOTAL_SIZE -gt $((CURRENT_LIMIT * 80 / 100)) ] && echo "Upgrade to GitHub Pro" || echo "Continue monitoring")" 