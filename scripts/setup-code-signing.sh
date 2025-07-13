#!/bin/bash

# 🍎 macOS Code Signing Setup Script for GeneKnow
# This script helps developers set up code signing for macOS releases

set -e

echo "🔐 GeneKnow macOS Code Signing Setup"
echo "===================================="
echo

# Check if running on macOS
if [[ "$OSTYPE" != "darwin"* ]]; then
    echo "❌ This script is for macOS only"
    exit 1
fi

# Check for required tools
if ! command -v security &> /dev/null; then
    echo "❌ security command not found"
    exit 1
fi

# Function to check if developer account is set up
check_developer_account() {
    echo "📋 Checking Apple Developer Account setup..."
    
    # Check for certificates
    CERT_COUNT=$(security find-identity -v -p codesigning | grep -c "Developer ID Application" || echo "0")
    
    if [[ $CERT_COUNT -eq 0 ]]; then
        echo "❌ No Developer ID Application certificates found"
        echo "📝 You need to:"
        echo "   1. Join Apple Developer Program (\$99/year)"
        echo "   2. Create a Developer ID Application certificate"
        echo "   3. Download and install it in Keychain"
        echo
        echo "🔗 More info: https://developer.apple.com/developer-id/"
        return 1
    else
        echo "✅ Found $CERT_COUNT Developer ID Application certificate(s)"
        security find-identity -v -p codesigning | grep "Developer ID Application"
        return 0
    fi
}

# Function to update Tauri config with signing identity
update_tauri_config() {
    local cert_name="$1"
    local config_file="desktop/src-tauri/tauri.conf.json"
    
    echo "📝 Updating Tauri configuration..."
    
    # Create backup
    cp "$config_file" "$config_file.backup"
    
    # Update signingIdentity using jq if available, otherwise sed
    if command -v jq &> /dev/null; then
        tmp=$(mktemp)
        jq ".bundle.macOS.signingIdentity = \"$cert_name\"" "$config_file" > "$tmp"
        mv "$tmp" "$config_file"
        echo "✅ Updated signingIdentity to: $cert_name"
    else
        echo "⚠️  jq not found, please manually update signingIdentity in $config_file"
        echo "   Set: \"signingIdentity\": \"$cert_name\""
    fi
}

# Function to set up environment variables
setup_environment() {
    echo "🔧 Setting up environment variables..."
    
    # Apple ID for notarization
    read -p "📧 Enter your Apple ID email: " APPLE_ID
    
    # App-specific password
    echo "🔑 You'll need an app-specific password for notarization"
    echo "   Create one at: https://appleid.apple.com/account/manage"
    read -s -p "🔐 Enter your app-specific password: " APP_PASSWORD
    echo
    
    # Team ID (optional)
    read -p "👥 Enter your Team ID (optional): " TEAM_ID
    
    # Create .env file for local development
    cat > .env.local << EOF
# Apple Developer Configuration
APPLE_ID=${APPLE_ID}
APPLE_PASSWORD=${APP_PASSWORD}
APPLE_TEAM_ID=${TEAM_ID}

# Tauri Environment Variables
TAURI_SIGNING_IDENTITY=Developer ID Application
TAURI_SIGNING_IDENTITY_PASSWORD=${APP_PASSWORD}
EOF
    
    echo "✅ Environment variables saved to .env.local"
    echo "⚠️  Keep this file secure and don't commit it to git!"
}

# Function to test code signing
test_code_signing() {
    echo "🧪 Testing code signing setup..."
    
    # Build the app
    echo "🔨 Building GeneKnow..."
    cd desktop/ui
    
    if ! pnpm run tauri-build; then
        echo "❌ Build failed"
        return 1
    fi
    
    # Check if the app is signed
    local app_path="desktop/src-tauri/target/release/bundle/macos/GeneKnow.app"
    
    if [[ -d "$app_path" ]]; then
        echo "🔍 Checking code signature..."
        codesign -v -v "$app_path"
        
        if [[ $? -eq 0 ]]; then
            echo "✅ Code signing successful!"
            
            # Check notarization (if configured)
            echo "🧾 Checking notarization status..."
            spctl -a -v "$app_path"
            
            if [[ $? -eq 0 ]]; then
                echo "✅ App passes Gatekeeper checks!"
            else
                echo "⚠️  App needs notarization for distribution"
                echo "   Run: xcrun altool --notarize-app --file GeneKnow.dmg"
            fi
        else
            echo "❌ Code signing failed"
            return 1
        fi
    else
        echo "❌ Built app not found at $app_path"
        return 1
    fi
}

# Function to show next steps
show_next_steps() {
    echo
    echo "🎉 Setup complete! Next steps:"
    echo "==============================="
    echo
    echo "1. 🔨 Build signed releases:"
    echo "   cd desktop/ui && pnpm run tauri-build"
    echo
    echo "2. 📦 Create notarized DMG:"
    echo "   xcrun altool --notarize-app --file GeneKnow.dmg"
    echo
    echo "3. 🚀 Distribute:"
    echo "   - Upload to GitHub releases"
    echo "   - Update download website"
    echo "   - Test on clean macOS system"
    echo
    echo "4. 📝 Update documentation:"
    echo "   - Remove bypass instructions"
    echo "   - Update installation guide"
    echo
    echo "🔗 Useful links:"
    echo "   - Apple Developer: https://developer.apple.com"
    echo "   - Tauri Signing: https://tauri.app/v1/guides/distribution/sign"
    echo "   - Notarization: https://developer.apple.com/documentation/notarization"
}

# Main execution
main() {
    echo "Starting code signing setup..."
    echo
    
    # Check if developer account is set up
    if ! check_developer_account; then
        echo "❌ Please set up your Apple Developer account first"
        exit 1
    fi
    
    # Get certificate name
    echo
    echo "📋 Available certificates:"
    security find-identity -v -p codesigning | grep "Developer ID Application"
    echo
    
    read -p "📝 Enter the certificate name to use: " CERT_NAME
    
    if [[ -z "$CERT_NAME" ]]; then
        echo "❌ Certificate name required"
        exit 1
    fi
    
    # Update configuration
    update_tauri_config "$CERT_NAME"
    
    # Set up environment
    setup_environment
    
    # Test the setup
    echo
    read -p "🧪 Test the setup by building now? (y/n): " -n 1 -r
    echo
    
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        test_code_signing
    fi
    
    # Show next steps
    show_next_steps
}

# Run main function
main "$@" 