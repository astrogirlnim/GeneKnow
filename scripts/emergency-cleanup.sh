#!/bin/bash

# 🚨 EMERGENCY GitHub Actions Artifact Cleanup
# Use this when you've hit the storage quota and need immediate relief

set -e

echo "🚨 EMERGENCY CLEANUP - GitHub Actions Storage Quota"
echo "=================================================="
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
echo "📊 Checking current repository..."
REPO_INFO=$(gh repo view --json nameWithOwner -q .nameWithOwner 2>/dev/null || echo "")

if [ -z "$REPO_INFO" ]; then
    echo "❌ Could not determine repository. Please run from a git repository."
    echo "   Or set the repository manually with: gh repo set-default owner/repo"
    exit 1
fi

echo "📁 Repository: $REPO_INFO"
echo ""

# Function to format bytes to human readable
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

# Check current storage usage
echo "📊 Checking storage usage..."
USAGE_JSON=$(gh api repos/$REPO_INFO/actions/cache/usage 2>/dev/null || echo "{}")
if [ "$USAGE_JSON" != "{}" ]; then
    TOTAL_BYTES=$(echo "$USAGE_JSON" | jq -r '.total_active_caches_size_in_bytes // 0')
    TOTAL_COUNT=$(echo "$USAGE_JSON" | jq -r '.total_active_caches_count // 0')
    echo "💾 Cache Storage: $(format_bytes $TOTAL_BYTES) across $TOTAL_COUNT caches"
fi

# Get all artifacts with detailed information
echo ""
echo "🔍 Fetching all artifacts..."
ARTIFACTS_JSON=$(gh api repos/$REPO_INFO/actions/artifacts --paginate 2>/dev/null || echo '{"artifacts": []}')
TOTAL_ARTIFACTS=$(echo "$ARTIFACTS_JSON" | jq -r '.artifacts | length')
TOTAL_SIZE=$(echo "$ARTIFACTS_JSON" | jq -r '[.artifacts[].size_in_bytes] | add // 0')

echo "📦 Found $TOTAL_ARTIFACTS artifacts using $(format_bytes $TOTAL_SIZE)"
echo ""

# Show artifact breakdown by workflow
echo "📊 Artifact breakdown by workflow:"
echo "$ARTIFACTS_JSON" | jq -r '.artifacts | group_by(.workflow_run.head_branch) | .[] | {branch: .[0].workflow_run.head_branch, count: length, size: ([.[].size_in_bytes] | add)} | "  \(.branch): \(.count) artifacts, \(.size / 1048576 | floor)MB"' 2>/dev/null || true
echo ""

# Emergency cleanup options
echo "🚨 EMERGENCY CLEANUP OPTIONS:"
echo "=============================="
echo ""
echo "1) Delete ALL artifacts (Nuclear option)"
echo "2) Keep only 1 most recent artifact"
echo "3) Delete artifacts older than 1 hour"
echo "4) Delete specific workflow artifacts"
echo "5) Custom cleanup (you specify criteria)"
echo "6) Exit without changes"
echo ""

read -p "Select option (1-6): " OPTION

case $OPTION in
    1)
        echo ""
        echo "⚠️  WARNING: This will delete ALL artifacts!"
        read -p "Are you SURE? Type 'DELETE ALL' to confirm: " CONFIRM
        if [ "$CONFIRM" = "DELETE ALL" ]; then
            echo "🗑️  Deleting ALL artifacts..."
            ARTIFACT_IDS=$(echo "$ARTIFACTS_JSON" | jq -r '.artifacts[].id')
            COUNT=0
            for id in $ARTIFACT_IDS; do
                if gh api -X DELETE repos/$REPO_INFO/actions/artifacts/$id 2>/dev/null; then
                    COUNT=$((COUNT + 1))
                    echo -ne "\r🗑️  Deleted $COUNT/$TOTAL_ARTIFACTS artifacts..."
                fi
            done
            echo ""
            echo "✅ Deleted $COUNT artifacts"
        else
            echo "❌ Cancelled"
        fi
        ;;
    
    2)
        echo ""
        echo "🔒 Keeping only the most recent artifact..."
        ARTIFACT_IDS=$(echo "$ARTIFACTS_JSON" | jq -r '.artifacts | sort_by(.created_at) | reverse | .[1:] | .[].id')
        COUNT=0
        TOTAL_TO_DELETE=$(echo "$ARTIFACT_IDS" | wc -l | xargs)
        for id in $ARTIFACT_IDS; do
            if gh api -X DELETE repos/$REPO_INFO/actions/artifacts/$id 2>/dev/null; then
                COUNT=$((COUNT + 1))
                echo -ne "\r🗑️  Deleted $COUNT/$TOTAL_TO_DELETE artifacts..."
            fi
        done
        echo ""
        echo "✅ Deleted $COUNT artifacts, kept 1 most recent"
        ;;
    
    3)
        echo ""
        echo "🕐 Deleting artifacts older than 1 hour..."
        ONE_HOUR_AGO=$(date -u -v-1H '+%Y-%m-%dT%H:%M:%SZ' 2>/dev/null || date -u -d '1 hour ago' '+%Y-%m-%dT%H:%M:%SZ')
        ARTIFACT_IDS=$(echo "$ARTIFACTS_JSON" | jq -r --arg date "$ONE_HOUR_AGO" '.artifacts[] | select(.created_at < $date) | .id')
        COUNT=0
        TOTAL_TO_DELETE=$(echo "$ARTIFACT_IDS" | wc -l | xargs)
        for id in $ARTIFACT_IDS; do
            if [ -n "$id" ] && gh api -X DELETE repos/$REPO_INFO/actions/artifacts/$id 2>/dev/null; then
                COUNT=$((COUNT + 1))
                echo -ne "\r🗑️  Deleted $COUNT artifacts older than 1 hour..."
            fi
        done
        echo ""
        echo "✅ Deleted $COUNT artifacts"
        ;;
    
    4)
        echo ""
        echo "📋 Available workflows:"
        echo "$ARTIFACTS_JSON" | jq -r '.artifacts | map(.workflow_run.name) | unique | .[]' | nl
        read -p "Enter workflow number to delete its artifacts: " WORKFLOW_NUM
        WORKFLOW_NAME=$(echo "$ARTIFACTS_JSON" | jq -r '.artifacts | map(.workflow_run.name) | unique | .[]' | sed -n "${WORKFLOW_NUM}p")
        if [ -n "$WORKFLOW_NAME" ]; then
            echo "🗑️  Deleting artifacts from workflow: $WORKFLOW_NAME"
            ARTIFACT_IDS=$(echo "$ARTIFACTS_JSON" | jq -r --arg wf "$WORKFLOW_NAME" '.artifacts[] | select(.workflow_run.name == $wf) | .id')
            COUNT=0
            for id in $ARTIFACT_IDS; do
                if gh api -X DELETE repos/$REPO_INFO/actions/artifacts/$id 2>/dev/null; then
                    COUNT=$((COUNT + 1))
                    echo -ne "\r🗑️  Deleted $COUNT artifacts..."
                fi
            done
            echo ""
            echo "✅ Deleted $COUNT artifacts from $WORKFLOW_NAME"
        else
            echo "❌ Invalid selection"
        fi
        ;;
    
    5)
        echo ""
        echo "📝 Custom cleanup - Enter your criteria"
        echo "Examples:"
        echo "  - Delete by name pattern: enter a keyword"
        echo "  - Delete by size: artifacts larger than X MB"
        read -p "Enter cleanup criteria: " CRITERIA
        # Add custom cleanup logic here based on criteria
        echo "⚠️  Custom cleanup not fully implemented in this version"
        ;;
    
    6)
        echo "👍 Exiting without changes"
        exit 0
        ;;
    
    *)
        echo "❌ Invalid option"
        exit 1
        ;;
esac

# Final storage check
echo ""
echo "📊 Checking final storage usage..."
FINAL_ARTIFACTS=$(gh api repos/$REPO_INFO/actions/artifacts --jq '.total_count' 2>/dev/null || echo "unknown")
echo "✅ Remaining artifacts: $FINAL_ARTIFACTS"
echo ""
echo "🚀 Storage cleanup completed!"
echo ""
echo "💡 Tips to prevent future issues:"
echo "  1. Use shorter retention periods (1-3 days)"
echo "  2. Only upload essential files"
echo "  3. Use maximum compression"
echo "  4. Consider GitHub Pro for more storage"
echo "  5. Run cleanup before each release build" 