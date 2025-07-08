import { HashRouter as Router, Routes, Route } from 'react-router-dom';
import OnboardingPage from './pages/OnboardingPage';
// import WelcomePage from './pages/WelcomePage';
// import FeaturesPage from './pages/FeaturesPage';
// import HowItWorksPage from './pages/HowItWorksPage';
// import PrivacyPage from './pages/PrivacyPage';
import Layout from './components/Layout';
import './App.css';
import './index.css';

function App() {
  return (
    <Router>
      <Layout>
        <Routes>
          <Route path="/" element={<OnboardingPage />} />
          {/* You can add back other routes here if needed */}
          {/* 
          <Route path="/welcome" element={<WelcomePage />} />
          <Route path="/features" element={<FeaturesPage />} />
          <Route path="/how-it-works" element={<HowItWorksPage />} />
          <Route path="/privacy" element={<PrivacyPage />} />
          */}
        </Routes>
      </Layout>
    </Router>
  );
}

export default App;
