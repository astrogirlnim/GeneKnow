#!/bin/bash

# 🔧 GitHub Actions Cost Optimization Script
# Helps monitor and optimize GitHub Actions usage to stay within budget

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}🔧 GitHub Actions Cost Optimization Tool${NC}"
echo "================================================"

# Check if GitHub CLI is installed
if ! command -v gh &> /dev/null; then
    echo -e "${RED}❌ GitHub CLI (gh) is not installed${NC}"
    echo "Install it from: https://cli.github.com/"
    exit 1
fi

# Check if authenticated
if ! gh auth status &> /dev/null; then
    echo -e "${RED}❌ Not authenticated with GitHub CLI${NC}"
    echo "Run: gh auth login"
    exit 1
fi

echo -e "${GREEN}✅ GitHub CLI is ready${NC}"
echo ""

# Get repository info
REPO_INFO=$(gh repo view --json name,owner)
REPO_NAME=$(echo $REPO_INFO | jq -r '.name')
REPO_OWNER=$(echo $REPO_INFO | jq -r '.owner.login')

echo -e "${BLUE}📊 Repository: ${REPO_OWNER}/${REPO_NAME}${NC}"
echo ""

# Function to show usage
show_usage() {
    echo -e "${YELLOW}📊 Current GitHub Actions Usage:${NC}"
    echo "=================================="
    
    # Get billing info (requires admin access)
    if gh api user/settings/billing/actions 2>/dev/null | jq -r '.total_minutes_used' &> /dev/null; then
        BILLING_INFO=$(gh api user/settings/billing/actions)
        TOTAL_MINUTES=$(echo $BILLING_INFO | jq -r '.total_minutes_used')
        INCLUDED_MINUTES=$(echo $BILLING_INFO | jq -r '.included_minutes')
        TOTAL_PAID=$(echo $BILLING_INFO | jq -r '.total_paid_for_minutes')
        
        echo -e "Total minutes used: ${RED}${TOTAL_MINUTES}${NC}"
        echo -e "Included minutes: ${GREEN}${INCLUDED_MINUTES}${NC}"
        echo -e "Paid minutes: ${YELLOW}${TOTAL_PAID}${NC}"
        
        # Calculate cost estimate
        if [ "$TOTAL_MINUTES" -gt "$INCLUDED_MINUTES" ]; then
            OVERAGE=$((TOTAL_MINUTES - INCLUDED_MINUTES))
            echo -e "⚠️  Over limit by: ${RED}${OVERAGE} minutes${NC}"
        else
            REMAINING=$((INCLUDED_MINUTES - TOTAL_MINUTES))
            echo -e "✅ Remaining: ${GREEN}${REMAINING} minutes${NC}"
        fi
    else
        echo -e "${YELLOW}⚠️  Cannot access billing info (requires admin access)${NC}"
    fi
    echo ""
}

# Function to show recent workflow runs
show_recent_runs() {
    echo -e "${YELLOW}📋 Recent Workflow Runs (Last 10):${NC}"
    echo "=================================="
    
    gh run list --limit 10 --json status,conclusion,name,createdAt,html_url | \
    jq -r '.[] | "\(.status) | \(.conclusion // "running") | \(.name) | \(.createdAt) | \(.html_url)"' | \
    while IFS='|' read -r status conclusion name created_at url; do
        case $conclusion in
            "success") icon="✅" ;;
            "failure") icon="❌" ;;
            "cancelled") icon="⏹️" ;;
            *) icon="🔄" ;;
        esac
        echo -e "${icon} ${name} (${created_at})"
    done
    echo ""
}

# Function to clean up old artifacts
cleanup_artifacts() {
    echo -e "${YELLOW}🧹 Cleaning up old artifacts...${NC}"
    echo "=================================="
    
    # List current artifacts
    ARTIFACTS=$(gh api repos/:owner/:repo/actions/artifacts --jq '.artifacts[] | select(.expired == false) | .id')
    ARTIFACT_COUNT=$(echo "$ARTIFACTS" | wc -l)
    
    if [ -z "$ARTIFACTS" ]; then
        echo -e "${GREEN}✅ No artifacts to clean up${NC}"
        return
    fi
    
    echo -e "Found ${YELLOW}${ARTIFACT_COUNT}${NC} artifacts"
    
    # Keep only the 3 most recent
    ARTIFACTS_TO_DELETE=$(echo "$ARTIFACTS" | tail -n +4)
    DELETE_COUNT=$(echo "$ARTIFACTS_TO_DELETE" | wc -l)
    
    if [ -z "$ARTIFACTS_TO_DELETE" ]; then
        echo -e "${GREEN}✅ Only keeping 3 most recent artifacts${NC}"
        return
    fi
    
    echo -e "Deleting ${RED}${DELETE_COUNT}${NC} old artifacts..."
    
    echo "$ARTIFACTS_TO_DELETE" | while read -r artifact_id; do
        if [ -n "$artifact_id" ]; then
            gh api -X DELETE "repos/:owner/:repo/actions/artifacts/$artifact_id"
            echo -e "  ${RED}❌${NC} Deleted artifact $artifact_id"
        fi
    done
    
    echo -e "${GREEN}✅ Artifact cleanup complete${NC}"
    echo ""
}

# Function to show cost optimization tips
show_optimization_tips() {
    echo -e "${BLUE}💡 Cost Optimization Tips:${NC}"
    echo "========================="
    echo ""
    echo -e "${GREEN}✅ Already Implemented:${NC}"
    echo "• Linux-only builds (10x cheaper than macOS)"
    echo "• Aggressive artifact cleanup (1-day retention)"
    echo "• Maximum compression for artifacts"
    echo "• Caching for dependencies"
    echo ""
    echo -e "${YELLOW}🔧 Additional Optimizations:${NC}"
    echo "• Use 'workflow_dispatch' for manual releases only"
    echo "• Skip CI on documentation-only changes"
    echo "• Use 'pull_request' paths to limit when workflows run"
    echo "• Consider self-hosted runners for heavy builds"
    echo ""
    echo -e "${RED}💰 Cost Breakdown (per minute):${NC}"
    echo "• Linux: \$0.008/minute"
    echo "• Windows: \$0.016/minute (2x more expensive)"
    echo "• macOS: \$0.08/minute (10x more expensive!)"
    echo ""
}

# Function to enable/disable expensive builds
toggle_expensive_builds() {
    echo -e "${BLUE}🔧 Toggle Expensive Builds${NC}"
    echo "========================="
    echo ""
    echo "Current status: Linux-only builds (cost-optimized)"
    echo ""
    echo "Options:"
    echo "1. Keep Linux-only (recommended for budget)"
    echo "2. Enable all platforms (high cost)"
    echo "3. Enable Windows only (medium cost)"
    echo "4. Enable macOS only (highest cost)"
    echo ""
    read -p "Choose option (1-4): " choice
    
    case $choice in
        1)
            echo -e "${GREEN}✅ Keeping Linux-only builds${NC}"
            ;;
        2)
            echo -e "${YELLOW}⚠️  This will significantly increase costs!${NC}"
            read -p "Are you sure? (y/N): " confirm
            if [[ $confirm =~ ^[Yy]$ ]]; then
                echo "To enable all platforms, uncomment the matrix entries in:"
                echo "• .github/workflows/release.yml"
                echo "• .github/workflows/pr.yml"
                echo -e "${RED}💸 Expected cost: ~\$8-12 per release${NC}"
            fi
            ;;
        3)
            echo -e "${YELLOW}⚠️  This will double your costs compared to Linux-only${NC}"
            echo "Uncomment Windows entries in workflow files"
            echo -e "${YELLOW}💸 Expected cost: ~\$2-4 per release${NC}"
            ;;
        4)
            echo -e "${RED}⚠️  This is the most expensive option!${NC}"
            echo "Uncomment macOS entries in workflow files"
            echo -e "${RED}💸 Expected cost: ~\$6-10 per release${NC}"
            ;;
        *)
            echo -e "${RED}❌ Invalid option${NC}"
            ;;
    esac
    echo ""
}

# Main menu
while true; do
    echo -e "${BLUE}🎯 What would you like to do?${NC}"
    echo "1. 📊 Show current usage"
    echo "2. 📋 Show recent workflow runs"
    echo "3. 🧹 Clean up old artifacts"
    echo "4. 💡 Show optimization tips"
    echo "5. 🔧 Toggle expensive builds"
    echo "6. 🚪 Exit"
    echo ""
    read -p "Choose option (1-6): " choice
    
    case $choice in
        1) show_usage ;;
        2) show_recent_runs ;;
        3) cleanup_artifacts ;;
        4) show_optimization_tips ;;
        5) toggle_expensive_builds ;;
        6) echo -e "${GREEN}👋 Goodbye!${NC}"; exit 0 ;;
        *) echo -e "${RED}❌ Invalid option${NC}" ;;
    esac
    
    echo ""
    read -p "Press Enter to continue..."
    echo ""
done 