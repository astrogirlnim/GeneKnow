import { BrowserRouter as Router, Routes, Route } from 'react-router-dom'
import WelcomePage from './pages/WelcomePage'
import FeaturesPage from './pages/FeaturesPage'
import PrivacyPage from './pages/PrivacyPage'
import HowItWorksPage from './pages/HowItWorksPage'
import UploadPage from './pages/UploadPage'
import DashboardPage from './pages/DashboardPage'
import ClinicalViewPage from './pages/ClinicalViewPage'
import ConfidenceCheckTest from './components/ConfidenceCheckTest' // TEMPORARY FOR TESTING
import { apiConfig } from './api/apiConfig';
import { listen } from '@tauri-apps/api/event';
import { useEffect } from 'react';

// --- Main App Component ---
export default function App() {
  useEffect(() => {
    // Initialize API configuration
    apiConfig.initialize().catch(console.error);

    // Listen for server status events
    const setupListeners = async () => {
      const unlistenReady = await listen('api-server-ready', () => {
        console.log('API server is ready');
        // Re-initialize API config to get the port
        apiConfig.reset();
        apiConfig.initialize().catch(console.error);
      });

      const unlistenError = await listen<string>('api-server-error', (event) => {
        console.error('API server error:', event.payload);
      });

      // Cleanup listeners on unmount
      return () => {
        unlistenReady();
        unlistenError();
      };
    };

    setupListeners();
  }, []);

  return (
    <Router>
      <Routes>
        <Route path="/" element={<WelcomePage />} />
        <Route path="/features" element={<FeaturesPage />} />
        <Route path="/privacy" element={<PrivacyPage />} />
        <Route path="/how-it-works" element={<HowItWorksPage />} />
        <Route path="/upload" element={<UploadPage />} />
        <Route path="/dashboard" element={<DashboardPage />} />
        <Route path="/clinical" element={<ClinicalViewPage />} />
        <Route path="/test-confidence" element={<ConfidenceCheckTest />} /> {/* TEMPORARY FOR TESTING */}
      </Routes>
    </Router>
  )
}

