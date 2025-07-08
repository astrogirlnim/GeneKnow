import React, { useState } from 'react';
import UploadOverlay from './UploadOverlay';
import DeepDive from './DeepDive';

// Type definitions
interface Alert {
  type: 'critical' | 'warning' | 'info' | 'success';
  title: string;
  desc: string;
}

interface PatientData {
  name: string;
  age: number;
  sex: string;
  condition: string;
  riskScore: string;
  details: string;
  alerts: Alert[];
}

const UploadPage: React.FC = () => {
  const [uploadStatus, setUploadStatus] = useState('ready');
  const [isDashboardEnabled, setIsDashboardEnabled] = useState(false);
  const [currentPatient, setCurrentPatient] = useState<PatientData | null>(null);

  const handleUploadComplete = (data: PatientData) => {
    setCurrentPatient(data);
    setIsDashboardEnabled(true);
    setUploadStatus('complete');
  };

  return (
    <div>
      {isDashboardEnabled && currentPatient ? (
        <DeepDive 
          currentPatient={currentPatient}
          uploadStatus={uploadStatus}
        />
      ) : (
        <UploadOverlay 
          onUploadComplete={handleUploadComplete}
        />
      )}
    </div>
  );
};

export default UploadPage; 