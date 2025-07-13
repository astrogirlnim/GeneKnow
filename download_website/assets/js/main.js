// Main JavaScript functionality for GeneKnow website

document.addEventListener('DOMContentLoaded', function() {
    // Initialize all functionality
    initMobileMenu();
    initPlatformDetection();
    initSmoothScrolling();
    initScrollEffects();
    initAnimations();
});

// Mobile menu functionality
function initMobileMenu() {
    const mobileToggle = document.querySelector('.mobile-menu-toggle');
    const nav = document.querySelector('.nav');

    if (mobileToggle && nav) {
        mobileToggle.addEventListener('click', function() {
            nav.classList.toggle('mobile-open');
            this.classList.toggle('active');
        });

        // Close mobile menu when clicking on nav links
        const navLinks = nav.querySelectorAll('.nav-link');
        navLinks.forEach(link => {
            link.addEventListener('click', function() {
                nav.classList.remove('mobile-open');
                mobileToggle.classList.remove('active');
            });
        });

        // Close mobile menu when clicking outside
        document.addEventListener('click', function(event) {
            if (!nav.contains(event.target) && !mobileToggle.contains(event.target)) {
                nav.classList.remove('mobile-open');
                mobileToggle.classList.remove('active');
            }
        });
    }
}

// Platform detection
function initPlatformDetection() {
    const platform = detectPlatform();
    const detectedOS = document.getElementById('detected-os');
    const recommendedOS = document.getElementById('recommended-os');

    if (detectedOS && recommendedOS) {
        const platformInfo = getPlatformInfo(platform);
        detectedOS.textContent = platformInfo.name;
        recommendedOS.textContent = platformInfo.name;
    }
}

function detectPlatform() {
    const userAgent = navigator.userAgent.toLowerCase();

    if (userAgent.includes('windows')) {
        return 'windows';
    } else if (userAgent.includes('mac')) {
        return 'macos';
    } else if (userAgent.includes('linux')) {
        return 'linux';
    } else {
        return 'unknown';
    }
}

function getPlatformInfo(platform) {
    const platformMap = {
        'windows': { name: 'Windows' },
        'macos': { name: 'macOS' },
        'linux': { name: 'Linux' },
        'unknown': { name: 'Your OS' }
    };

    return platformMap[platform] || platformMap['unknown'];
}

// Smooth scrolling for anchor links
function initSmoothScrolling() {
    const links = document.querySelectorAll('a[href^="#"]');

    links.forEach(link => {
        link.addEventListener('click', function(e) {
            e.preventDefault();

            const targetId = this.getAttribute('href');
            const targetElement = document.querySelector(targetId);

            if (targetElement) {
                const offsetTop = targetElement.offsetTop - 80; // Account for fixed header

                window.scrollTo({
                    top: offsetTop,
                    behavior: 'smooth'
                });
            }
        });
    });
}

// Scroll effects
function initScrollEffects() {
    const header = document.querySelector('.header');

    window.addEventListener('scroll', function() {
        const scrollTop = window.pageYOffset;

        // Keep header consistently white
        if (header) {
            header.style.background = '#FFFFFF';
            header.style.backdropFilter = 'none';
            if (scrollTop > 50) {
                header.style.boxShadow = '0 2px 4px -1px rgba(0, 0, 0, 0.06)';
            } else {
                header.style.boxShadow = '0 1px 2px 0 rgba(0, 0, 0, 0.03)';
            }
        }
    });
}

// Animations and effects
function initAnimations() {
    // Intersection Observer for fade-in animations
    const observerOptions = {
        threshold: 0.1,
        rootMargin: '0px 0px -50px 0px'
    };

    const observer = new IntersectionObserver(function(entries) {
        entries.forEach(entry => {
            if (entry.isIntersecting) {
                entry.target.classList.add('fade-in');
            }
        });
    }, observerOptions);

    // Observe elements for animation
    const animatedElements = document.querySelectorAll(
        '.feature-card, .download-card, .stat, .about-feature'
    );

    animatedElements.forEach(element => {
        observer.observe(element);
    });
}

// Utility functions
function formatBytes(bytes, decimals = 2) {
    if (bytes === 0) return '0 Bytes';

    const k = 1024;
    const dm = decimals < 0 ? 0 : decimals;
    const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB'];

    const i = Math.floor(Math.log(bytes) / Math.log(k));

    return parseFloat((bytes / Math.pow(k, i)).toFixed(dm)) + ' ' + sizes[i];
}

function formatDate(dateString) {
    const date = new Date(dateString);
    return date.toLocaleDateString('en-US', {
        year: 'numeric',
        month: 'long',
        day: 'numeric'
    });
}

// Error handling
function showError(message) {
    console.error('GeneKnow Website Error:', message);

    // You can add more sophisticated error handling here
    // For example, showing a toast notification or modal
}

// Loading states
function showLoading(element) {
    if (element) {
        element.classList.add('loading');
    }
}

function hideLoading(element) {
    if (element) {
        element.classList.remove('loading');
    }
}

// Light mode only - dark theme detection disabled

// Analytics and tracking (placeholder)
function trackEvent(event, category, action, label) {
    // Add your analytics tracking code here
    console.log('Track Event:', { event, category, action, label });
}

// Download tracking
function trackDownload(platform, fileName) {
    trackEvent('download', 'file', 'download', `${platform}-${fileName}`);
}

// Export functions for use in other scripts
window.GeneKnow = {
    formatBytes,
    formatDate,
    detectPlatform,
    getPlatformInfo,
    showError,
    showLoading,
    hideLoading,
    trackDownload,
    trackEvent
};

// Add fade-in CSS class
const style = document.createElement('style');
style.textContent = `
    .fade-in {
        animation: fadeIn 0.6s ease-in-out;
    }
    
    @keyframes fadeIn {
        from {
            opacity: 0;
            transform: translateY(20px);
        }
        to {
            opacity: 1;
            transform: translateY(0);
        }
    }
`;
document.head.appendChild(style);
