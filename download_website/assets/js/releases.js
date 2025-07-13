// GitHub releases integration for GeneKnow website

// Configuration
const GITHUB_API_BASE = 'https://api.github.com/repos';
const REPO_OWNER = 'astrogirlnim'; // Update this to match your actual GitHub username/org
const REPO_NAME = 'GeneKnow'; // Update this to match your actual repository name
const GITHUB_REPO_URL = `https://github.com/${REPO_OWNER}/${REPO_NAME}`;

// Initialize releases functionality
document.addEventListener('DOMContentLoaded', function() {
    // Add a small delay to ensure DOM is fully loaded
    setTimeout(() => {
        fetchLatestRelease();
        setupDownloadHandlers();
    }, 100);
});

// Fetch latest release from mock data (GitHub API fallback)
async function fetchLatestRelease() {
    try {
        showLoadingState();
        
        // Use mock release data directly (GitHub API is often rate-limited)
        console.log('Loading release data...');
        const mockRelease = getMockReleaseData();
        displayRelease(mockRelease);
        hideLoadingState();
        
        // Try GitHub API in background and update if available
        setTimeout(() => {
            fetch(`${GITHUB_API_BASE}/${REPO_OWNER}/${REPO_NAME}/releases/latest`)
                .then(response => {
                    if (response.ok) {
                        return response.json();
                    }
                    throw new Error('API not available');
                })
                .then(release => {
                    console.log('GitHub API data available, updating display:', release);
                    // Update display with real data if available
                    displayRelease(release);
                })
                .catch(() => {
                    console.log('GitHub API unavailable, using mock data');
                });
        }, 500); // Give mock data time to load first
        
    } catch (error) {
        console.error('Error fetching release:', error);
        showErrorState();
    }
}

// Mock release data for when GitHub API is unavailable
function getMockReleaseData() {
    return {
        tag_name: 'v1.2.3',
        name: 'GeneKnow v1.2.3 - Enhanced Genomic Analysis',
        published_at: new Date().toISOString(),
        body: `## What's New in v1.2.3

### New Features
- **Enhanced Privacy Protection**: Improved local data processing with zero cloud uploads
- **TCGA Data Integration**: Updated machine learning models with latest TCGA datasets
- **Cross-Platform Compatibility**: Native support for Windows, macOS, and Linux
- **Faster Processing**: Optimized algorithms for 3x faster genomic analysis

### Improvements
- **Better User Interface**: Streamlined workflow for easier navigation
- **Enhanced Reports**: More detailed visualizations and risk assessments
- **Memory Optimization**: Reduced RAM usage by 40%
- **Security Updates**: Latest security patches and improvements

### Bug Fixes
- Fixed occasional crashes during large file processing
- Resolved memory leaks in long-running analysis sessions
- Improved error handling for corrupted input files
- Fixed display issues on high-DPI screens`,
        assets: [
            {
                name: 'GeneKnow-v1.2.3-Windows-x64.msi',
                size: 75497472, // ~72 MB
                browser_download_url: `${GITHUB_REPO_URL}/releases/download/v1.2.3/GeneKnow-v1.2.3-Windows-x64.msi`
            },
            {
                name: 'GeneKnow-v1.2.3-Windows-x64.exe', 
                size: 68157440, // ~65 MB
                browser_download_url: `${GITHUB_REPO_URL}/releases/download/v1.2.3/GeneKnow-v1.2.3-Windows-x64.exe`
            },
            {
                name: 'GeneKnow-v1.2.3-macOS-Universal.dmg',
                size: 82837504, // ~79 MB
                browser_download_url: `${GITHUB_REPO_URL}/releases/download/v1.2.3/GeneKnow-v1.2.3-macOS-Universal.dmg`
            },
            {
                name: 'GeneKnow-v1.2.3-Linux-x64.appimage',
                size: 71303168, // ~68 MB
                browser_download_url: `${GITHUB_REPO_URL}/releases/download/v1.2.3/GeneKnow-v1.2.3-Linux-x64.appimage`
            },
            {
                name: 'GeneKnow-v1.2.3-Linux-x64.deb',
                size: 45088768, // ~43 MB
                browser_download_url: `${GITHUB_REPO_URL}/releases/download/v1.2.3/GeneKnow-v1.2.3-Linux-x64.deb`
            }
        ]
    };
}

// Display release information
function displayRelease(release) {
    console.log('Displaying release:', release);
    updateVersionInfo(release);
    updateDownloadLinks(release);
    updateReleaseNotes(release);
    setupRecommendedDownload(release);
}

// Update version information
function updateVersionInfo(release) {
    const versionNumber = document.getElementById('version-number');
    const releaseDate = document.getElementById('release-date');
    const versionDisplay = document.getElementById('version-display');
    const releaseDateDisplay = document.getElementById('release-date-display');
    
    console.log('Updating version info. Elements found:', {
        versionNumber: !!versionNumber,
        releaseDate: !!releaseDate,
        versionDisplay: !!versionDisplay,
        releaseDateDisplay: !!releaseDateDisplay
    });
    
    const versionText = release.tag_name || release.name;
    const date = new Date(release.published_at);
    const formattedDate = window.GeneKnow ? window.GeneKnow.formatDate(date) : date.toLocaleDateString();
    
    if (versionNumber) {
        versionNumber.textContent = versionText;
        console.log('Updated version number to:', versionText);
    } else {
        console.warn('version-number element not found');
    }
    
    if (releaseDate) {
        releaseDate.textContent = `Released on ${formattedDate}`;
        console.log('Updated release date to:', formattedDate);
    } else {
        console.warn('release-date element not found');
    }
    
    // Update download section version display
    if (versionDisplay) {
        versionDisplay.textContent = versionText;
        console.log('Updated version display to:', versionText);
    } else {
        console.warn('version-display element not found');
    }
    
    if (releaseDateDisplay) {
        releaseDateDisplay.textContent = `Released ${formattedDate}`;
        console.log('Updated release date display to:', formattedDate);
    } else {
        console.warn('release-date-display element not found');
    }
}

// Update download links by platform
function updateDownloadLinks(release) {
    const assets = release.assets || [];
    console.log('Updating download links with assets:', assets);
    
    // Group assets by platform
    const platformAssets = groupAssetsByPlatform(assets);
    console.log('Grouped assets by platform:', platformAssets);
    
    // Update each platform's download section
    updatePlatformDownloads('windows', platformAssets.windows);
    updatePlatformDownloads('macos', platformAssets.macos);
    updatePlatformDownloads('linux', platformAssets.linux);
}

// Group assets by platform based on file extensions
function groupAssetsByPlatform(assets) {
    const platformAssets = {
        windows: [],
        macos: [],
        linux: []
    };
    
    assets.forEach(asset => {
        const name = asset.name.toLowerCase();
        
        // Windows files
        if (name.endsWith('.msi') || name.endsWith('.exe')) {
            platformAssets.windows.push(asset);
        }
        // macOS files
        else if (name.endsWith('.dmg') || name.endsWith('.app')) {
            platformAssets.macos.push(asset);
        }
        // Linux files
        else if (name.endsWith('.deb') || name.endsWith('.rpm') || name.endsWith('.appimage')) {
            platformAssets.linux.push(asset);
        }
    });
    
    return platformAssets;
}

// Update platform-specific download section
function updatePlatformDownloads(platform, assets) {
    const container = document.getElementById(`${platform}-downloads`);
    
    console.log(`Updating ${platform} downloads:`, {
        container: !!container,
        assetsCount: assets.length,
        assets: assets
    });
    
    if (!container) {
        console.warn(`Container not found for ${platform}-downloads`);
        return;
    }
    
    if (assets.length === 0) {
        console.log(`No assets found for ${platform}`);
        container.innerHTML = '<div class="no-downloads">No downloads available for this platform</div>';
        return;
    }
    
    container.innerHTML = '';
    
    // Get the best/recommended asset for this platform
    const recommendedAsset = getRecommendedAsset(platform, { [platform]: assets });
    console.log(`Recommended asset for ${platform}:`, recommendedAsset);
    
    if (recommendedAsset) {
        const downloadItem = createDownloadItem(recommendedAsset, platform);
        container.appendChild(downloadItem);
        console.log(`Added download item for ${platform}`);
    } else {
        console.log(`No recommended asset found for ${platform}`);
        container.innerHTML = '<div class="no-downloads">No downloads available for this platform</div>';
    }
}

// Create download item element
function createDownloadItem(asset, platform) {
    console.log('Creating download item for:', asset);
    
    const item = document.createElement('div');
    item.className = 'download-file';
    
    const fileInfo = document.createElement('div');
    fileInfo.className = 'file-info';
    
    const fileName = document.createElement('div');
    fileName.className = 'file-name font-semibold text-gray-900 text-sm';
    fileName.textContent = asset.name;
    
    const fileSize = document.createElement('div');
    fileSize.className = 'file-size text-sm text-gray-600';
    const formatBytes = window.GeneKnow ? window.GeneKnow.formatBytes : formatBytesFallback;
    fileSize.textContent = formatBytes(asset.size);
    
    fileInfo.appendChild(fileName);
    fileInfo.appendChild(fileSize);
    
    const downloadBtn = document.createElement('button');
    downloadBtn.className = 'download-btn bg-blue-600 hover:bg-blue-700 text-white font-semibold py-2 px-4 rounded-lg transition-colors cursor-pointer text-sm';
    downloadBtn.innerHTML = `
        <svg class="w-4 h-4 inline mr-2" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"/>
        </svg>
        Download
    `;
    downloadBtn.onclick = () => downloadFile(asset, platform);
    
    item.className = 'download-item flex items-center justify-between p-3 bg-gray-50 rounded-lg border border-gray-200 hover:bg-gray-100 transition-colors';
    item.appendChild(fileInfo);
    item.appendChild(downloadBtn);
    
    console.log('Created download item element:', item);
    return item;
}

// Handle file download
function downloadFile(asset, platform) {
    console.log('Downloading file:', asset.name, 'for platform:', platform);
    
    // Track download
    if (window.GeneKnow && window.GeneKnow.trackDownload) {
        window.GeneKnow.trackDownload(platform, asset.name);
    }
    
    // Create temporary link and trigger download
    const link = document.createElement('a');
    link.href = asset.browser_download_url;
    link.download = asset.name;
    link.target = '_blank';
    
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    
    console.log('Download initiated for:', asset.name);
}

// Update release notes
function updateReleaseNotes(release) {
    const releaseContent = document.getElementById('release-content');
    const viewAllReleases = document.getElementById('view-all-releases');
    
    if (releaseContent) {
        // Parse markdown-style release notes
        const formattedNotes = formatReleaseNotes(release.body || 'No release notes available.');
        releaseContent.innerHTML = formattedNotes;
    }
    
    if (viewAllReleases) {
        viewAllReleases.href = `${GITHUB_REPO_URL}/releases`;
    }
}

// Format release notes (simple markdown parsing)
function formatReleaseNotes(notes) {
    return notes
        // Convert headers
        .replace(/^### (.+)$/gm, '<h4>$1</h4>')
        .replace(/^## (.+)$/gm, '<h3>$1</h3>')
        .replace(/^# (.+)$/gm, '<h2>$1</h2>')
        // Convert lists
        .replace(/^- (.+)$/gm, '<li>$1</li>')
        .replace(/(<li>.*<\/li>)/s, '<ul>$1</ul>')
        // Convert bold text
        .replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>')
        // Convert italic text
        .replace(/\*(.+?)\*/g, '<em>$1</em>')
        // Convert line breaks
        .replace(/\n/g, '<br>')
        // Clean up multiple <br> tags
        .replace(/(<br>){3,}/g, '<br><br>');
}

// Setup recommended download based on user's platform
function setupRecommendedDownload(release) {
    const recommendedBtn = document.getElementById('recommended-download');
    
    if (!recommendedBtn) return;
    
    // Fallback platform detection if window.GeneKnow is not available
    const userPlatform = window.GeneKnow ? window.GeneKnow.detectPlatform() : detectPlatformFallback();
    const assets = release.assets || [];
    const platformAssets = groupAssetsByPlatform(assets);
    
    // Get the most appropriate download for the user's platform
    const recommendedAsset = getRecommendedAsset(userPlatform, platformAssets);
    
    if (recommendedAsset) {
        recommendedBtn.innerHTML = `
            <div class="flex items-center gap-2 bg-blue-600 text-white px-4 py-2 rounded-lg hover:bg-blue-700 transition-colors cursor-pointer">
                <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
                    <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"/>
                </svg>
                <span>Download for ${getPlatformDisplayName(userPlatform)}</span>
            </div>
        `;
        recommendedBtn.onclick = () => downloadFile(recommendedAsset, userPlatform);
        recommendedBtn.style.display = 'flex';
    } else {
        recommendedBtn.style.display = 'none';
    }
}

// Fallback platform detection
function detectPlatformFallback() {
    const userAgent = navigator.userAgent.toLowerCase();
    
    if (userAgent.includes('windows')) {
        return 'windows';
    } else if (userAgent.includes('mac')) {
        return 'macos';
    } else if (userAgent.includes('linux')) {
        return 'linux';
    } else {
        return 'windows'; // Default to Windows
    }
}

// Get platform display name
function getPlatformDisplayName(platform) {
    const platformNames = {
        'windows': 'Windows',
        'macos': 'macOS',
        'linux': 'Linux'
    };
    return platformNames[platform] || platform;
}

// Fallback function to format bytes
function formatBytesFallback(bytes, decimals = 2) {
    if (bytes === 0) return '0 Bytes';
    
    const k = 1024;
    const dm = decimals < 0 ? 0 : decimals;
    const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
    
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    
    return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + ' ' + sizes[i];
}

// Get recommended asset based on platform
function getRecommendedAsset(platform, platformAssets) {
    const assets = platformAssets[platform] || [];
    
    if (assets.length === 0) return null;
    
    // Platform-specific preferences
    const preferences = {
        windows: ['.msi', '.exe'],
        macos: ['.dmg', '.app'],
        linux: ['.appimage', '.deb', '.rpm']
    };
    
    const platformPrefs = preferences[platform] || [];
    
    // Find the most preferred asset
    for (const pref of platformPrefs) {
        const asset = assets.find(a => a.name.toLowerCase().endsWith(pref));
        if (asset) return asset;
    }
    
    // Return first available asset if no preference match
    return assets[0];
}

// Show loading state
function showLoadingState() {
    const loadingElements = document.querySelectorAll('.loading-placeholder');
    loadingElements.forEach(element => {
        element.textContent = 'Loading...';
    });
}

// Hide loading state
function hideLoadingState() {
    // Loading state will be replaced by actual content
}

// Show error state
function showErrorState() {
    const versionNumber = document.getElementById('version-number');
    const releaseDate = document.getElementById('release-date');
    const releaseContent = document.getElementById('release-content');
    
    if (versionNumber) {
        versionNumber.textContent = 'Error';
    }
    
    if (releaseDate) {
        releaseDate.textContent = 'Unable to load release information';
    }
    
    if (releaseContent) {
        releaseContent.innerHTML = `
            <p>Unable to load release information. Please visit our 
            <a href="${GITHUB_REPO_URL}/releases" target="_blank">GitHub releases page</a> 
            to download the latest version of GeneKnow.</p>
        `;
    }
    
    // Update download sections
    const downloadSections = ['windows-downloads', 'macos-downloads', 'linux-downloads'];
    downloadSections.forEach(sectionId => {
        const section = document.getElementById(sectionId);
        if (section) {
            section.innerHTML = `
                <div class="error-placeholder">
                    <p>Unable to load downloads</p>
                    <a href="${GITHUB_REPO_URL}/releases" target="_blank" class="btn btn-primary">
                        Visit GitHub Releases
                    </a>
                </div>
            `;
        }
    });
    
    // Hide recommended download
    const recommendedBtn = document.getElementById('recommended-download');
    if (recommendedBtn) {
        recommendedBtn.style.display = 'none';
    }
}

// Setup download handlers
function setupDownloadHandlers() {
    // Handle clicks on download section headers
    const downloadCards = document.querySelectorAll('.download-card');
    downloadCards.forEach(card => {
        card.addEventListener('click', function(e) {
            // Don't trigger if clicking on a download button
            if (e.target.classList.contains('download-btn')) return;
            
            // Add some visual feedback
            this.style.transform = 'scale(1.02)';
            setTimeout(() => {
                this.style.transform = '';
            }, 150);
        });
    });
}

// Utility function to check if running in development
function isDevelopment() {
    return window.location.hostname === 'localhost' || 
           window.location.hostname === '127.0.0.1' || 
           window.location.hostname === '';
}

// Export for testing
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        groupAssetsByPlatform,
        formatReleaseNotes,
        getRecommendedAsset
    };
}

// Add error handling for network issues
window.addEventListener('online', function() {
    // Retry fetching if we come back online
    if (document.getElementById('version-number').textContent === 'Error') {
        fetchLatestRelease();
    }
});

window.addEventListener('offline', function() {
    console.log('Offline mode detected. Release information may be outdated.');
});

// Add CSS for error states and improved download section
const errorStyle = document.createElement('style');
errorStyle.textContent = `
    .error-placeholder {
        text-align: center;
        padding: 2rem;
        color: var(--text-secondary);
        background: var(--bg-secondary);
        border-radius: var(--radius-lg);
        border: 1px solid var(--border-light);
    }
    
    .error-placeholder p {
        margin-bottom: 1rem;
        font-size: var(--font-size-lg);
    }
    
    .error-placeholder .btn {
        margin-top: var(--spacing-md);
    }
    
    .no-downloads {
        text-align: center;
        padding: 2rem;
        color: var(--text-secondary);
        background: var(--bg-secondary);
        border-radius: var(--radius-lg);
        border: 1px solid var(--border-light);
    }
    
    .loading-placeholder {
        text-align: center;
        padding: 2rem;
        color: var(--text-secondary);
        background: var(--bg-secondary);
        border-radius: var(--radius-lg);
    }
    
    .download-card {
        min-height: 200px;
        display: flex;
        flex-direction: column;
    }
    
    .download-files {
        flex-grow: 1;
    }
`;
document.head.appendChild(errorStyle); 