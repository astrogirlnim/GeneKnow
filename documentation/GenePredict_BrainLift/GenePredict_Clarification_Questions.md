
# GenePredict – Clarification Questions

## 🧬 Genomic Model Scope
1. Should the prototype focus only on **breast cancer risk** or support multiple diseases from day one?
2. Should we:
   - Use **precomputed variants** from tools like DeepVariant (faster)?
   - Or support **raw VCF-to-prediction** (more impressive, but complex)?

## 🖥️ User Experience
3. Should the clinician UI include **EHR-style filters** (e.g., filter by risk score, condition, or patient age)?
4. Will this version be **clinician-only**, or should patients interact with the tool directly as well?

## 🌐 NLP Use
5. Should Llama 3.1 reports include **actionable clinical recommendations** (e.g., “Schedule annual screening”)?
6. Do you want reports to include **citations or references** for transparency?

## 🧪 Evaluation & Success
7. How will we define **model performance**? (e.g., AUROC, F1 score, ClinVar accuracy benchmark)
8. Would you like a **feedback mechanism** (e.g., clinicians rate usefulness or interpretability of reports)?
