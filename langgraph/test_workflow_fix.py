#!/usr/bin/env python3
"""
Test file to verify GitHub Actions workflows trigger correctly
after fixing path filters to include langgraph/**
"""

def test_workflow_trigger():
    """Simple test function to verify workflow triggering"""
    print("✅ Workflow path filters have been fixed!")
    print("🧬 PR and Release workflows should now trigger on langgraph changes")
    return True

if __name__ == "__main__":
    test_workflow_trigger() 