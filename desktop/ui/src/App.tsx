import { BrowserRouter as Router, Routes, Route } from 'react-router-dom'
import WelcomePage from './pages/WelcomePage'
import FeaturesPage from './pages/FeaturesPage'
import PrivacyPage from './pages/PrivacyPage'
import HowItWorksPage from './pages/HowItWorksPage'
import UploadPage from './pages/UploadPage'
import DashboardPage from './pages/DashboardPage'
import ClinicalViewPage from './pages/ClinicalViewPage'

// --- Main App Component ---
export default function App() {
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
      </Routes>
    </Router>
  )
}

