# Installing GeneKnow from ZIP on macOS

## Overview

Starting with version 0.1.4, GeneKnow for macOS is distributed as a ZIP archive instead of a DMG file. This change helps preserve the bundled Python runtime and resources that are required for the app to function properly.

## Installation Steps

### 1. Download the ZIP File

1. Visit the [GeneKnow releases page](https://github.com/astrogirlnim/GeneKnow/releases/latest)
2. Download the file named `GeneKnow-[version]-macos-[arch].zip`
   - For Intel Macs: `x86_64` 
   - For Apple Silicon (M1/M2/M3): `aarch64`

### 2. Extract the Application

1. **Locate the downloaded ZIP file** (usually in your Downloads folder)
2. **Double-click the ZIP file** to extract it
   - macOS will automatically extract the contents
   - You should see `GeneKnow.app` appear

### 3. Remove Quarantine (Required for Unsigned Apps)

Since GeneKnow is not yet code-signed with an Apple Developer ID, you need to remove the quarantine flag:

1. **Open Terminal** (found in Applications → Utilities)
2. **Run this command**:
   ```bash
   xattr -cr ~/Downloads/GeneKnow.app
   ```
   Replace `~/Downloads/GeneKnow.app` with the actual path if you extracted it elsewhere.

### 4. Move to Applications

1. **Drag `GeneKnow.app`** to your Applications folder
2. The app is now installed

### 5. First Launch

For the first launch, you need to bypass Gatekeeper:

1. **Right-click** on GeneKnow.app in Applications
2. Select **"Open"** from the context menu
3. In the security dialog, click **"Open"**

After the first launch, you can open GeneKnow normally by double-clicking.

## Verification

To verify the app was installed correctly with all resources:

1. Open Terminal
2. Run:
   ```bash
   ls -la /Applications/GeneKnow.app/Contents/Resources/_up_/bundled_resources/
   ```
3. You should see:
   - `start_api_server.sh`
   - `python_runtime/`
   - `geneknow_pipeline/`
   - `desktop/`

## Troubleshooting

### "GeneKnow is damaged and can't be opened"

This message appears because the app is not code-signed. To fix:

1. Make sure you removed quarantine: `xattr -cr /Applications/GeneKnow.app`
2. Try right-clicking and selecting "Open" instead of double-clicking

### Python Runtime Missing

If the app launches but shows errors about missing Python:

1. The ZIP extraction may have failed to preserve the bundled resources
2. Try re-downloading and extracting the ZIP file
3. Make sure you're using macOS's built-in extraction (not a third-party tool)

### App Doesn't Launch

Check Console.app for error messages:
1. Open Console (Applications → Utilities)
2. Filter for "GeneKnow"
3. Look for error messages when trying to launch

## Why ZIP Instead of DMG?

We switched from DMG to ZIP distribution because:

1. **Preserves Resources**: macOS Gatekeeper strips executable content from unsigned DMGs during installation
2. **Simpler Process**: No mounting/unmounting required
3. **Better Compatibility**: Works consistently across different macOS versions
4. **Maintains Structure**: The app bundle structure remains intact with all Python runtime files

## Security Note

GeneKnow is currently distributed without an Apple Developer ID certificate. This means:

- You'll see security warnings on first launch
- You need to explicitly allow the app to run
- The app hasn't been notarized by Apple

We're working on obtaining proper code signing certificates for future releases to provide a smoother installation experience.

## Future Improvements

In upcoming releases, we plan to:

1. Obtain an Apple Developer ID for proper code signing
2. Implement app notarization for enhanced security
3. Potentially offer Homebrew installation as an alternative
4. Return to DMG distribution once signing is implemented

## Need Help?

If you encounter issues:

1. Check our [GitHub Issues](https://github.com/astrogirlnim/GeneKnow/issues)
2. Review the [main installation guide](../README.md)
3. Contact support with:
   - Your macOS version
   - The exact error message
   - Console.app logs related to GeneKnow 