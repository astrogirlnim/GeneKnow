#!/bin/bash

# 🧹 Quick GitHub Actions Artifact Cleanup
# Run this script to immediately free up storage space

set -e

echo "🧹 Quick artifact cleanup for GitHub Actions storage quota..."

# Get repository information
if ! command -v gh &> /dev/null; then
    echo "❌ GitHub CLI (gh) not found. Please install it first:"
    echo "   brew install gh"
    exit 1
fi

# Get all artifact IDs and delete them (except the most recent 2)
echo "📦 Fetching artifact list..."
ARTIFACT_IDS=$(gh api repos/:owner/:repo/actions/artifacts --jq '.artifacts | sort_by(.created_at) | reverse | .[2:] | .[].id')

if [ -z "$ARTIFACT_IDS" ]; then
    echo "✅ No old artifacts to delete"
    exit 0
fi

echo "🗑️  Deleting old artifacts..."
COUNT=0
for id in $ARTIFACT_IDS; do
    if gh api -X DELETE repos/:owner/:repo/actions/artifacts/$id 2>/dev/null; then
        COUNT=$((COUNT + 1))
        echo "  ✅ Deleted artifact $id"
    else
        echo "  ⚠️  Failed to delete artifact $id"
    fi
done

echo "📊 Cleanup complete: Deleted $COUNT artifacts"
echo "🚀 Storage space freed up - you can now retry your build!" 