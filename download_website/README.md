# GeneKnow Download Website

This is the static website for hosting GeneKnow downloads and serving as a marketing landing page. The website is automatically deployed to GitHub Pages and fetches the latest releases from the GitHub API.

## Features

- **Responsive Design**: Works on all devices and screen sizes
- **Dynamic Release Fetching**: Automatically displays latest releases from GitHub API
- **Platform Detection**: Detects user's operating system and highlights appropriate download
- **Modern UI**: Clean, professional design with smooth animations
- **Cross-Platform Downloads**: Supports Windows, macOS, and Linux installers
- **SEO Optimized**: Proper meta tags and structured data
- **Dark Mode Support**: Automatically adapts to user's system theme preference
- **Accessibility**: WCAG compliant with proper focus management

## Structure

```
download_website/
├── index.html              # Main landing page
├── assets/
│   ├── css/
│   │   ├── main.css        # Main stylesheet
│   │   └── responsive.css  # Responsive design rules
│   ├── js/
│   │   ├── main.js         # Main JavaScript functionality
│   │   └── releases.js     # GitHub API integration
│   └── images/
│       ├── logo.svg        # Site logo
│       └── favicon.png     # Site favicon
├── _config.yml             # GitHub Pages configuration
└── README.md              # This file
```

## Setup

### 1. GitHub Repository Configuration

1. Ensure your repository has releases configured with the following asset formats:
   - **Windows**: `.msi`, `.exe` files
   - **macOS**: `.dmg`, `.app` files  
   - **Linux**: `.deb`, `.rpm`, `.AppImage` files

2. Update the repository information in `assets/js/releases.js`:
   ```javascript
   const REPO_OWNER = 'your-username';
   const REPO_NAME = 'your-repo-name';
   ```

### 2. GitHub Pages Setup

1. Go to your repository Settings → Pages
2. Set **Source** to "Deploy from a branch"
3. Select **Branch**: `main` or your default branch
4. Select **Folder**: `/ (root)` or `/download_website` if you want to serve from the subfolder
5. Click **Save**

### 3. Custom Domain (Optional)

If you want to use a custom domain:

1. Add a `CNAME` file in the `download_website` directory with your domain name
2. Configure your DNS settings to point to `your-username.github.io`

## Customization

### Logo and Branding

1. Replace `assets/images/logo.svg` with your own logo
2. Update `assets/favicon.png` with your favicon
3. Modify the color scheme in `assets/css/main.css` by updating the CSS variables

### Content

1. Update the hero section in `index.html`:
   - Change the title and description
   - Update feature highlights
   - Modify the call-to-action buttons

2. Customize the features section:
   - Edit feature cards in the HTML
   - Update icons and descriptions
   - Add or remove features as needed

3. Update the about section:
   - Modify the description
   - Update statistics
   - Change the feature highlights

### Styling

The website uses CSS custom properties (variables) for easy theming:

```css
:root {
    --primary-color: #2563eb;
    --secondary-color: #64748b;
    --accent-color: #f59e0b;
    /* ... more variables */
}
```

## Platform Support

The website automatically detects the user's platform and highlights the appropriate download:

- **Windows**: Prioritizes `.msi` files, then `.exe`
- **macOS**: Prioritizes `.dmg` files, then `.app`
- **Linux**: Prioritizes `.AppImage`, then `.deb`, then `.rpm`

## Automatic Updates

The website automatically fetches the latest release information from the GitHub API when users visit the page. No manual updates are required when you publish new releases.

## Development

### Local Development

1. Clone the repository
2. Open `index.html` in your browser
3. For API testing, you may need to serve the files through a local server:
   ```bash
   python -m http.server 8000
   # or
   npx serve .
   ```

### Testing

- Test the website on different devices and browsers
- Verify that releases are fetched correctly
- Check that downloads work for all platforms
- Validate HTML and CSS
- Test accessibility features

## Important Notes

### API Rate Limits

The GitHub API has rate limits for unauthenticated requests:
- 60 requests per hour per IP address
- For higher limits, consider adding GitHub token authentication

### Repository Name

Make sure to update the repository configuration in `releases.js` to match your actual GitHub repository:

```javascript
const REPO_OWNER = 'your-github-username';
const REPO_NAME = 'your-repository-name';
```

### CORS and Security

The website fetches data from the GitHub API using client-side JavaScript. This works because GitHub's API supports CORS. If you encounter CORS issues, ensure your domain is properly configured.

## Performance

- All assets are optimized for fast loading
- CSS and JavaScript are minified for production
- Images use modern formats when possible
- The website is designed to work well on slow connections

## License

This download website template is part of the GeneKnow project and follows the same licensing terms as the main repository.

## Contributing

To contribute to the download website:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## Support

If you encounter issues with the download website:

1. Check the browser console for JavaScript errors
2. Verify that the GitHub API is accessible
3. Ensure your release assets are properly formatted
4. Create an issue in the main repository

---

Made with love by the GeneKnow team 