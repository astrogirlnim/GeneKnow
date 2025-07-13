# 🍎 macOS Installation Guide for GeneKnow

This guide helps macOS users install and run GeneKnow while handling security prompts.

## 🔒 Understanding the Security Warning

When you first try to open GeneKnow on macOS, you may see this error:

> **"GeneKnow" Not Opened**  
> Apple could not verify "GeneKnow" is free of malware that may harm your Mac or compromise your privacy.

**This is normal and expected!** Here's why:

- **GeneKnow is currently unsigned** - We don't have an Apple Developer Certificate yet
- **macOS Gatekeeper** protects against unsigned applications by default
- **It's completely safe** - GeneKnow is open-source and processes data locally only

## ✅ How to Install GeneKnow Safely

### Method 1: Right-Click to Open (Recommended)

1. **Don't double-click** the GeneKnow app initially
2. **Right-click** (or Control-click) on the GeneKnow application
3. Select **"Open"** from the context menu
4. You'll see a different dialog with an **"Open"** button
5. Click **"Open"** to confirm
6. **GeneKnow will now run and be trusted permanently**

### Method 2: System Preferences Override

1. Try to open GeneKnow normally (it will be blocked)
2. Go to **System Preferences** → **Security & Privacy** → **General**
3. You'll see: *"GeneKnow was blocked from use because it is not from an identified developer"*
4. Click **"Open Anyway"** next to this message
5. Confirm by clicking **"Open"** in the dialog

### Method 3: Terminal Command (Advanced Users)

```bash
# Remove the quarantine attribute that triggers Gatekeeper
sudo xattr -rd com.apple.quarantine /path/to/GeneKnow.app

# Example if it's in Applications folder:
sudo xattr -rd com.apple.quarantine /Applications/GeneKnow.app
```

## 🛡️ Is This Safe?

**Yes, it's completely safe** to bypass Gatekeeper for GeneKnow because:

✅ **Open Source**: Full source code available on GitHub  
✅ **Privacy-First**: No data ever leaves your device  
✅ **Local Processing**: All genomic analysis runs on your Mac  
✅ **No Network Calls**: Works completely offline  
✅ **Transparent**: You can verify exactly what the app does  

## 🔐 Why Isn't GeneKnow Code Signed?

Code signing requires:
- **Apple Developer Account** ($99/year)
- **Certificate Management** (complex setup)
- **Notarization Process** (additional steps)

For open-source privacy-focused apps like GeneKnow, this creates barriers without significant security benefits since:
- The source code is publicly auditable
- The app runs completely locally
- No external data transmission occurs

## 🚀 First Launch Experience

After bypassing Gatekeeper, your first launch will:

1. **Initialize the local database** (1-3 minutes)
2. **Set up the Python runtime** (bundled with the app)
3. **Start the analysis pipeline** (runs locally on your Mac)
4. **Show the welcome screen** when ready

## 🛠️ Troubleshooting

### "App is damaged and can't be opened"
```bash
# Remove quarantine and try again
sudo xattr -rd com.apple.quarantine /Applications/GeneKnow.app
```

### "No permission to open the application"
```bash
# Fix permissions
sudo chmod +x /Applications/GeneKnow.app/Contents/MacOS/GeneKnow
```

### App won't start after bypass
1. Check Console.app for error messages
2. Restart your Mac
3. Try the terminal command method

## 📋 System Requirements

- **macOS 10.13** or later
- **8GB RAM** minimum (16GB recommended)
- **2GB free disk space**
- **Intel or Apple Silicon** Mac

## 🔄 Future Updates

Once we implement proper code signing:
- No more security warnings
- Automatic updates through the app
- App Store distribution (planned)

## 🆘 Still Having Issues?

If you continue having problems:

1. **Check our GitHub Issues**: [github.com/astrogirlnim/GeneKnow/issues](https://github.com/astrogirlnim/GeneKnow/issues)
2. **Create a new issue** with:
   - Your macOS version
   - Error messages
   - Steps you've tried
3. **Email support**: [Insert support email when available]

## 📚 Additional Resources

- [GeneKnow User Guide](./USER_GUIDE.md)
- [Privacy Policy](./PRIVACY_POLICY.md)
- [Source Code](https://github.com/astrogirlnim/GeneKnow)

---

**Remember**: GeneKnow is designed for complete privacy - your genetic data never leaves your Mac! 🔒 