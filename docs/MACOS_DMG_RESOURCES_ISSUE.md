# macOS DMG Bundled Resources Issue

## Problem Summary

When installing GeneKnow from an unsigned DMG on macOS, the bundled Python runtime and scripts are stripped out by macOS Gatekeeper security, causing the app to fail at startup.

### What's Happening

1. **Build Time**: The release pipeline successfully bundles 453MB of resources including:
   - `python_runtime/` - Complete Python 3.13 installation
   - `geneknow_pipeline/` - All pipeline code and ML models
   - `start_api_server.sh` - Startup script
   - ML models and databases

2. **DMG Creation**: Resources are correctly included in the DMG (verified in CI/CD)

3. **Installation**: When users install the DMG, macOS Gatekeeper:
   - Identifies the app as unsigned
   - Strips out "suspicious" executable content
   - Only preserves the `desktop/` subdirectory
   - Result: App crashes because Python runtime is missing

## Root Cause

The app is not code-signed with an Apple Developer ID certificate. macOS treats unsigned apps as potentially dangerous and removes executable content during installation.

Current `tauri.conf.json` settings:
```json
"macOS": {
  "signingIdentity": null,      // No signing certificate
  "hardenedRuntime": false,     // Not hardened
  "entitlements": null          // No entitlements
}
```

## Solutions

### Solution 1: Proper Code Signing (Recommended for Production)

1. **Get Apple Developer Account**
   - Cost: $99/year
   - Sign up at: https://developer.apple.com/programs/

2. **Create Developer ID Certificate**
   - In Xcode: Preferences → Accounts → Manage Certificates
   - Create "Developer ID Application" certificate

3. **Update tauri.conf.json**
   ```json
   "macOS": {
     "signingIdentity": "Developer ID Application: Your Name (TEAMID)",
     "hardenedRuntime": true,
     "entitlements": "./entitlements.plist"
   }
   ```

4. **Create entitlements.plist**
   ```xml
   <?xml version="1.0" encoding="UTF-8"?>
   <!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
   <plist version="1.0">
   <dict>
     <key>com.apple.security.cs.allow-unsigned-executable-memory</key>
     <true/>
     <key>com.apple.security.cs.allow-dyld-environment-variables</key>
     <true/>
   </dict>
   </plist>
   ```

5. **Notarize the App**
   - After building, submit to Apple for notarization
   - This allows the app to run without security warnings

### Solution 2: Manual Installation Workaround (For Testing)

For users testing unsigned builds:

1. **Install the DMG normally**

2. **Remove quarantine attributes**
   ```bash
   xattr -cr /Applications/GeneKnow.app
   ```

3. **Manually restore bundled resources**
   ```bash
   # Check if resources are missing
   ls -la /Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/
   
   # If missing, extract from DMG manually
   hdiutil attach GeneKnow_*.dmg
   cp -R /Volumes/GeneKnow/GeneKnow.app/Contents/Resources/_up_/bundled_resources/* \
        /Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/
   hdiutil detach /Volumes/GeneKnow
   ```

4. **First run**: Right-click the app and select "Open" (bypasses Gatekeeper)

### Solution 3: Alternative Distribution Methods

1. **Direct .app Bundle**
   - Distribute as a zipped .app instead of DMG
   - Users still need to remove quarantine: `xattr -cr GeneKnow.app`
   - Resources less likely to be stripped

2. **Homebrew Cask**
   - Create a Homebrew formula
   - Handles quarantine removal automatically
   - Better for technical users

3. **pkg Installer**
   - Use productbuild to create a .pkg installer
   - Can include post-install scripts to fix permissions
   - Still requires signing for smooth experience

## Verification Steps

To verify if resources are properly included:

```bash
# Check bundled resources in installed app
ls -la /Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/

# Should show:
# - start_api_server.sh
# - python_runtime/
# - geneknow_pipeline/
# - desktop/

# Check Python executable
/Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/python_runtime/bin/python3 --version

# Test startup script
/Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/start_api_server.sh
```

## Temporary Fix for Current Users

Until the app is properly signed, users can:

1. Download the release DMG
2. Mount it but don't drag to Applications yet
3. Open Terminal and run:
   ```bash
   # Copy app with resources intact
   sudo cp -R /Volumes/GeneKnow/GeneKnow.app /Applications/
   
   # Remove quarantine
   sudo xattr -cr /Applications/GeneKnow.app
   
   # Fix permissions
   sudo chmod -R 755 /Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/
   ```

## Long-term Solution

The sustainable solution is to implement proper code signing:

1. Set up Apple Developer account
2. Configure automatic signing in CI/CD
3. Implement notarization workflow
4. Test with macOS 10.15+ (Catalina and newer)

This will ensure smooth installation without manual intervention and maintain the security model that macOS users expect.

## Related Files

- `.github/workflows/release.yml` - Build pipeline
- `desktop/src-tauri/tauri.conf.json` - Tauri configuration
- `desktop/scripts/bundle-python-optimized.sh` - Resource bundling
- `desktop/scripts/setup-code-signing.sh` - Signing setup helper 