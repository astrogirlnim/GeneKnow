// GitHub releases integration for GenePredict website

// Configuration
const GITHUB_API_BASE = 'https://api.github.com/repos';
const REPO_OWNER = 'astrogirlnim'; // Update this to match your actual GitHub username/org
const REPO_NAME = 'GeneKnow'; // Update this to match your actual repository name
const GITHUB_REPO_URL = `https://github.com/${REPO_OWNER}/${REPO_NAME}`;

// Initialize releases functionality
document.addEventListener('DOMContentLoaded', function() {
    fetchLatestRelease();
    setupDownloadHandlers();
});

// Fetch latest release from GitHub API
async function fetchLatestRelease() {
    try {
        showLoadingState();
        
        const response = await fetch(`${GITHUB_API_BASE}/${REPO_OWNER}/${REPO_NAME}/releases/latest`);
        
        if (!response.ok) {
            throw new Error(`HTTP error! status: ${response.status}`);
        }
        
        const release = await response.json();
        displayRelease(release);
        hideLoadingState();
        
    } catch (error) {
        console.error('Error fetching release:', error);
        showErrorState();
    }
}

// Display release information
function displayRelease(release) {
    updateVersionInfo(release);
    updateDownloadLinks(release);
    updateReleaseNotes(release);
    setupRecommendedDownload(release);
}

// Update version information
function updateVersionInfo(release) {
    const versionNumber = document.getElementById('version-number');
    const releaseDate = document.getElementById('release-date');
    
    if (versionNumber) {
        versionNumber.textContent = release.tag_name;
    }
    
    if (releaseDate) {
        const date = new Date(release.published_at);
        releaseDate.textContent = `Released on ${window.GenePredict.formatDate(date)}`;
    }
}

// Update download links by platform
function updateDownloadLinks(release) {
    const assets = release.assets || [];
    
    // Group assets by platform
    const platformAssets = groupAssetsByPlatform(assets);
    
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
    
    if (!container) return;
    
    if (assets.length === 0) {
        container.innerHTML = '<div class="no-downloads">No downloads available for this platform</div>';
        return;
    }
    
    container.innerHTML = '';
    
    assets.forEach(asset => {
        const downloadItem = createDownloadItem(asset, platform);
        container.appendChild(downloadItem);
    });
}

// Create download item element
function createDownloadItem(asset, platform) {
    const item = document.createElement('div');
    item.className = 'download-file';
    
    const fileInfo = document.createElement('div');
    fileInfo.className = 'file-info';
    
    const fileName = document.createElement('div');
    fileName.className = 'file-name';
    fileName.textContent = asset.name;
    
    const fileSize = document.createElement('div');
    fileSize.className = 'file-size';
    fileSize.textContent = window.GenePredict.formatBytes(asset.size);
    
    fileInfo.appendChild(fileName);
    fileInfo.appendChild(fileSize);
    
    const downloadBtn = document.createElement('button');
    downloadBtn.className = 'download-btn';
    downloadBtn.textContent = 'Download';
    downloadBtn.onclick = () => downloadFile(asset, platform);
    
    item.appendChild(fileInfo);
    item.appendChild(downloadBtn);
    
    return item;
}

// Handle file download
function downloadFile(asset, platform) {
    // Track download
    window.GenePredict.trackDownload(platform, asset.name);
    
    // Create temporary link and trigger download
    const link = document.createElement('a');
    link.href = asset.browser_download_url;
    link.download = asset.name;
    link.target = '_blank';
    
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
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
    
    const userPlatform = window.GenePredict.detectPlatform();
    const assets = release.assets || [];
    const platformAssets = groupAssetsByPlatform(assets);
    
    // Get the most appropriate download for the user's platform
    const recommendedAsset = getRecommendedAsset(userPlatform, platformAssets);
    
    if (recommendedAsset) {
        recommendedBtn.onclick = () => downloadFile(recommendedAsset, userPlatform);
        recommendedBtn.style.display = 'flex';
    } else {
        recommendedBtn.style.display = 'none';
    }
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
            to download the latest version.</p>
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

// Add CSS for error states
const errorStyle = document.createElement('style');
errorStyle.textContent = `
    .error-placeholder {
        text-align: center;
        padding: 2rem;
        color: var(--text-secondary);
    }
    
    .error-placeholder p {
        margin-bottom: 1rem;
    }
    
    .no-downloads {
        text-align: center;
        padding: 2rem;
        color: var(--text-secondary);
        font-style: italic;
    }
    
    .download-file {
        transition: transform 0.2s ease;
    }
    
    .download-file:hover {
        transform: translateX(4px);
    }
`;
document.head.appendChild(errorStyle); 