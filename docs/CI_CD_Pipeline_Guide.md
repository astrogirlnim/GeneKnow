# 🧬 GeneKnow CI/CD Pipeline Guide

## 📋 Overview

The GeneKnow CI/CD pipeline is designed to ensure high-quality, secure, and automated deployment of our privacy-first genomic risk assessment platform. Our pipeline consists of two main workflows:

1. **🔍 PR Validation Pipeline** - Validates code quality, security, and cross-platform compatibility for pull requests
2. **🚀 Release Pipeline** - Automates versioning, building, and deploying production releases

## 🏗️ Pipeline Architecture

### 🔍 PR Validation Pipeline (`.github/workflows/pr.yml`)

**Trigger:** Pull requests to `main` or `develop` branches  
**Purpose:** Ensure code quality and prevent broken builds before merging

```mermaid
graph TD
    A[PR Created] --> B[Code Quality & Linting]
    B --> C[Cross-Platform Build Testing]
    B --> D[Security Audit]
    C --> E[Linux x64 Build]
    C --> F[Windows x64 Build]
    C --> G[macOS x64 Build]
    D --> H[Node.js Security Audit]
    D --> I[Rust Security Audit]
    E --> J[PR Summary]
    F --> J
    G --> J
    H --> J
    I --> J
```

**Jobs:**
- **🧪 Code Quality & Linting**: ESLint, TypeScript compilation, bundle size analysis
- **🔨 Build Testing**: Multi-platform debug builds (Linux, Windows, macOS)
- **🛡️ Security Audit**: Dependency vulnerability scanning for Node.js and Rust
- **📋 PR Summary**: Generates comprehensive validation results

### 🚀 Release Pipeline (`.github/workflows/release.yml`)

**Trigger:** Pushes to `main` branch or manual workflow dispatch  
**Purpose:** Automated versioning and production releases

### 🌐 Download Website Pipeline (`.github/workflows/deploy-website.yml`)

**Trigger:** Changes to `download_website/` files, manual dispatch, or completion of Release Pipeline  
**Purpose:** Automatically updates the download website with latest releases

```mermaid
graph TD
    A[Main Branch Push] --> B[Determine Version]
    B --> C[Update Version Files]
    C --> D[Commit & Tag]
    D --> E[Validation & Testing]
    E --> F[Multi-Platform Build]
    F --> G[Linux Production Build]
    F --> H[Windows Production Build]
    F --> I[macOS Production Build]
    G --> J[Create GitHub Release]
    H --> J
    I --> J
    J --> K[Upload Artifacts]
    J --> L[Trigger Website Update]
    L --> M[Deploy Latest Website]
```

**Jobs:**
- **🏷️ Version Management**: Semantic versioning based on conventional commits
- **🧪 Validation & Testing**: Final validation before production builds
- **🔨 Build Release**: Production builds for Linux, Windows, and macOS
- **🚀 Create GitHub Release**: GitHub release with installers and release notes

## 🚀 Usage Instructions

### For Developers

#### Creating a Pull Request
1. Create a feature branch from `develop`
2. Make your changes following conventional commit patterns
3. Push to GitHub and create a PR to `develop` or `main`
4. The PR pipeline will automatically run and validate your changes
5. Review the PR summary for validation results

#### Conventional Commit Messages
Our pipeline uses conventional commits for automatic versioning:

```bash
# Patch version (0.1.0 -> 0.1.1)
git commit -m "fix: resolve parsing issue with VCF files"
git commit -m "docs: update installation instructions"

# Minor version (0.1.0 -> 0.2.0)
git commit -m "feat: add support for BAM file format"
git commit -m "feat: implement risk heatmap visualization"

# Major version (0.1.0 -> 1.0.0)
git commit -m "feat!: change API structure for ML models"
git commit -m "feat: add new feature

BREAKING CHANGE: This changes the public API"
```

### For Maintainers

#### Automated Releases
1. Merge approved PRs to `main` branch
2. The release pipeline automatically:
   - Calculates the next version based on commit messages
   - Updates version files (`package.json`, `tauri.conf.json`, `Cargo.toml`)
   - Creates a git tag
   - Builds production releases for all platforms
   - Creates a GitHub release with installers

#### Manual Releases
1. Go to **Actions** → **Release Pipeline**
2. Click **Run workflow**
3. Select the version bump type (patch, minor, major)
4. Click **Run workflow**

## ⚙️ Configuration & Setup

### Required GitHub Secrets

For the pipeline to work correctly, configure these secrets in your repository:

#### **Optional: Code Signing (Recommended for Production)**
```bash
# For Tauri code signing (macOS/Windows)
TAURI_SIGNING_PRIVATE_KEY="-----BEGIN PRIVATE KEY-----\n...\n-----END PRIVATE KEY-----"
TAURI_SIGNING_PRIVATE_KEY_PASSWORD="your-key-password"
```

### Environment Variables

The pipeline uses these environment variables (configured in workflow files):

```yaml
NODE_VERSION: '20'          # Node.js version
RUST_VERSION: '1.77.2'      # Rust version
PNPM_VERSION: '9.15.1'      # pnpm version
```

### Platform-Specific Dependencies

#### Linux (Ubuntu)
```bash
sudo apt-get install -y libgtk-3-dev libwebkit2gtk-4.0-dev librsvg2-dev
```

#### macOS
- Xcode Command Line Tools (automatically installed)
- Universal binary support for Intel and Apple Silicon

#### Windows
- Microsoft Visual C++ Build Tools (automatically installed)
- Windows SDK (automatically installed)

## 🛠️ Troubleshooting

### Common Issues and Solutions

#### **Build Failures**

**Issue:** TypeScript compilation errors
```bash
error TS2304: Cannot find name 'SomeType'
```
**Solution:** Check imports and type definitions in `desktop/ui/src/`

**Issue:** Rust compilation errors
```bash
error[E0433]: failed to resolve: use of undeclared crate or module
```
**Solution:** Check dependencies in `desktop/src-tauri/Cargo.toml`

**Issue:** Frontend build failures
```bash
✘ [ERROR] Could not resolve "some-package"
```
**Solution:** Run `pnpm install` and ensure all dependencies are in `package.json`

#### **Security Audit Failures**

**Issue:** Node.js vulnerabilities
```bash
found 3 high severity vulnerabilities
```
**Solution:** Run `pnpm audit fix` or update vulnerable packages

**Issue:** Rust vulnerabilities
```bash
error: package 'some-crate' contains a security vulnerability
```
**Solution:** Update `Cargo.toml` dependencies or use `cargo audit fix`

#### **Cross-Platform Issues**

**Issue:** Linux builds failing
```bash
Package 'gtk+-3.0' not found
```
**Solution:** System dependencies are automatically installed in CI

**Issue:** macOS builds failing
```bash
error: linking with `cc` failed: exit status: 1
```
**Solution:** Ensure Xcode Command Line Tools are installed

**Issue:** Windows builds failing
```bash
LINK : fatal error LNK1181: cannot open input file
```
**Solution:** Check that all dependencies are compatible with Windows

### Debug Mode

To debug pipeline issues:

1. **Enable Debug Logging:**
   ```bash
   RUST_LOG=debug cargo tauri dev
   ```

2. **Check Workflow Logs:**
   - Go to **Actions** tab in GitHub
   - Click on the failed workflow
   - Expand the failed step to see detailed logs

3. **Local Testing:**
   ```bash
   # Test ESLint
   cd desktop/ui && pnpm lint

   # Test TypeScript compilation
   cd desktop/ui && pnpm exec tsc --noEmit

   # Test build
   cd desktop/ui && pnpm build

   # Test Tauri build
   cd desktop/ui && pnpm tauri-build
   ```

### Performance Optimization

#### **Caching Strategy**
Our pipeline uses multiple caching layers:
- **pnpm Cache**: Stores Node.js dependencies
- **Cargo Cache**: Stores Rust dependencies and build artifacts
- **GitHub Actions Cache**: Stores intermediate build files

#### **Build Optimization**
- **Debug builds** for PRs (faster feedback)
- **Release builds** only for main branch (slower but optimized)
- **Parallel builds** across platforms
- **Incremental compilation** for Rust

## 🔧 Maintenance

### Regular Tasks

#### **Monthly:**
- Review and update Node.js/Rust versions
- Update GitHub Actions to latest versions
- Review security audit results
- Clean up old release artifacts

#### **Quarterly:**
- Update dependencies to latest stable versions
- Review and optimize build times
- Update documentation
- Review CI/CD metrics and usage

#### **Annually:**
- Evaluate new CI/CD tools and practices
- Review security practices
- Update code signing certificates
- Audit pipeline permissions

### Monitoring and Metrics

#### **Key Metrics to Monitor:**
- Build success rate (target: >95%)
- Build time (target: <10 minutes)
- Test coverage (target: >80%)
- Security vulnerability count (target: 0 high/critical)

#### **Alerting:**
- Failed builds on main branch
- Security vulnerabilities
- Long build times (>15 minutes)
- Quota usage warnings

## 📚 Best Practices

### Code Quality
- Always run `pnpm lint` before committing
- Use TypeScript strict mode
- Write meaningful commit messages
- Keep pull requests focused and small

### Security
- Regularly update dependencies
- Use `pnpm audit` and `cargo audit`
- Keep secrets secure and rotate regularly
- Review third-party dependencies

### Performance
- Optimize bundle sizes
- Use appropriate caching strategies
- Monitor build times
- Profile application performance

### Documentation
- Update documentation with code changes
- Keep README files current
- Document breaking changes
- Maintain changelog

## 🚀 Advanced Features

### Custom Workflows

You can create custom workflows for specific needs:

```yaml
# .github/workflows/custom.yml
name: Custom Workflow
on:
  workflow_dispatch:
    inputs:
      environment:
        description: 'Environment to deploy'
        required: true
        default: 'staging'
        type: choice
        options:
          - staging
          - production
```

### Release Channels

Configure different release channels:

```yaml
# Alpha releases
- name: Alpha Release
  if: contains(github.ref, 'alpha')
  run: |
    echo "Building alpha release..."
```

### Integration Testing

Add integration tests to the pipeline:

```yaml
- name: Integration Tests
  run: |
    echo "Running integration tests..."
    # Add your integration test commands here
```

## 📞 Support

If you encounter issues with the CI/CD pipeline:

1. **Check this documentation** for common solutions
2. **Review workflow logs** in GitHub Actions
3. **Search existing issues** in the repository
4. **Create a new issue** with detailed information:
   - Workflow name and run ID
   - Error messages and logs
   - Steps to reproduce
   - Expected vs actual behavior

## 🏷️ Version History

- **v1.0.0**: Initial pipeline implementation
- **v1.1.0**: Added cross-platform build matrix
- **v1.2.0**: Enhanced security auditing
- **v1.3.0**: Added automated versioning
- **v1.4.0**: Improved caching and performance

---

🧬 **GeneKnow** - Privacy-First Genomic Risk Assessment Platform

*This documentation is maintained by the GeneKnow development team. Last updated: January 2025* 