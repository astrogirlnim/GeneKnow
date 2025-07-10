#!/bin/bash

# 🧹 Aggressive GitHub Actions Artifact Cleanup Script
# This script will free up storage space immediately by removing old artifacts

set -e

echo "🧹 Starting aggressive artifact cleanup..."

# Get repository info
REPO=$(gh repo view --json owner,name -q '.owner.login + "/" + .name')
echo "📍 Repository: $REPO"

# Get current storage usage
echo "📊 Current storage usage:"
gh api repos/$REPO/actions/cache/usage --jq '.total_active_cache_size_in_bytes / 1024 / 1024 | floor | tostring + " MB"'

# Get all artifacts
echo "📦 Fetching all artifacts..."
gh api repos/$REPO/actions/artifacts --paginate --jq '.artifacts[] | {id: .id, name: .name, size: (.size_in_bytes / 1024 / 1024 | floor), created_at: .created_at}' > /tmp/artifacts.json

# Count total artifacts
TOTAL_COUNT=$(cat /tmp/artifacts.json | wc -l)
echo "📊 Found $TOTAL_COUNT artifacts"

# Show top 10 largest artifacts
echo "📈 Top 10 largest artifacts:"
cat /tmp/artifacts.json | jq -r '. | [.name, .size, .created_at] | @tsv' | sort -k2 -nr | head -10 | while read name size created; do
    echo "  📦 $name: ${size}MB (created: $created)"
done

# Delete strategy: Keep only the most recent 3 artifacts
echo "🗑️  Deleting all but the most recent 3 artifacts..."

# Get artifact IDs to delete (all except the newest 3)
ARTIFACTS_TO_DELETE=$(cat /tmp/artifacts.json | jq -r '.id' | head -n -3)

DELETE_COUNT=0
for artifact_id in $ARTIFACTS_TO_DELETE; do
    artifact_name=$(cat /tmp/artifacts.json | jq -r "select(.id == $artifact_id) | .name")
    artifact_size=$(cat /tmp/artifacts.json | jq -r "select(.id == $artifact_id) | .size")
    
    echo "🗑️  Deleting: $artifact_name (${artifact_size}MB)"
    
    if gh api -X DELETE repos/$REPO/actions/artifacts/$artifact_id 2>/dev/null; then
        DELETE_COUNT=$((DELETE_COUNT + 1))
        echo "  ✅ Deleted successfully"
    else
        echo "  ⚠️  Failed to delete (may already be deleted)"
    fi
done

echo "📊 Cleanup summary:"
echo "  🗑️  Deleted $DELETE_COUNT artifacts"
echo "  🔒 Kept the 3 most recent artifacts"

# Check storage usage after cleanup
echo "📊 Storage usage after cleanup:"
gh api repos/$REPO/actions/cache/usage --jq '.total_active_cache_size_in_bytes / 1024 / 1024 | floor | tostring + " MB"'

echo "✅ Aggressive cleanup completed!"
echo "🚀 You can now re-run your build pipeline"

# Cleanup temp files
rm -f /tmp/artifacts.json 