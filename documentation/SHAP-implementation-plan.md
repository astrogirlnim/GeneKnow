# **SHAP Sanity-Check Node Implementation Plan**

## **Objective**

To create a robust, self-contained node that validates the output of the Risk Model. This node will analyze the *reasoning* behind each prediction using SHAP and flag results that are statistically plausible but biologically questionable.

## **Node Signature**

* **Inputs**:  
  1. model\_object: The trained TensorFlow/PyTorch model object.  
  2. patient\_feature\_vector: A 1D NumPy array or Pandas Series containing the feature values for the single patient being scored.  
  3. feature\_names: A list of strings corresponding to the columns of the feature vector. **Crucially, these names must contain metadata (e.g., variant\_BRCA1\_p.185delAG\_ClinVarP)**.  
  4. predicted\_risk\_score: The float output from the Risk Model node.  
* **Output**:  
  * A single JSON object containing the prediction\_score and the validation block, as previously designed.

## **Core Implementation Steps**

\# Pseudocode for the node's main function

def run\_sanity\_check(model, patient\_vector, feature\_names, risk\_score):  
    \# 1\. Initialize SHAP Explainer  
    \# Use a pre-computed background dataset for reference  
    background\_data \= load\_background\_dataset('path/to/background.npy')  
    explainer \= shap.DeepExplainer(model, background\_data)

    \# 2\. Calculate SHAP values for the single patient instance  
    shap\_values \= explainer.shap\_values(patient\_vector)\[0\] \# Index \[0\] for single-output models

    \# 3\. Get Top Contributing Features  
    \# Combine feature names, their actual values, and their SHAP contributions  
    contributors \= get\_top\_contributors(feature\_names, patient\_vector, shap\_values, top\_n=5)

    \# 4\. Apply Sanity Rules Engine  
    validation\_status, reason\_code, explanation \= apply\_sanity\_rules(  
        risk\_score, patient\_vector, feature\_names, contributors  
    )

    \# 5\. Construct Final JSON Output  
    output\_json \= format\_output(  
        risk\_score, validation\_status, reason\_code, explanation, contributors  
    )

    return output\_json

## **Key Decisions & Technical Details**

1. **SHAP Explainer Choice**:  
   * **Decision**: Use shap.DeepExplainer.  
   * **Reason**: It's specifically optimized for TensorFlow and PyTorch models, making it significantly faster than model-agnostic explainers like KernelExplainer.  
2. **Background Dataset**:  
   * **Decision**: Create a representative background dataset from the training data. This is **the most critical prerequisite**.  
   * **Action Item**: Before implementing the node, run a K-Means clustering algorithm (k=100) on the training set and use the cluster centroids as the background data. This provides a robust baseline for SHAP calculations. Save this as a .npy file to be loaded by the node.  
3. **Feature Naming Convention**:  
   * **Decision**: Enforce a strict feature naming convention in the upstream Feature Vector Builder node.  
   * **Format**: type\_name\_metadata.  
   * **Examples**:  
     * variant\_BRCA1\_p.185delAG\_ClinVarP (P for Pathogenic)  
     * variant\_rs12345\_ClinVarB (B for Benign)  
     * PRS\_BreastCancer  
4. **Sanity Rules Engine (apply\_sanity\_rules function)**:  
   * **High-Risk Rule Logic**:  
     * IF risk\_score \> 0.6:  
     * AND no feature in top\_contributors contains '\_ClinVarP':  
     * RETURN 'FLAG\_FOR\_REVIEW', 'HIGH\_RISK\_WEAK\_EVIDENCE', ...  
   * **Low-Risk Rule Logic**:  
     * IF risk\_score \< 0.1:  
     * AND any feature in patient\_feature\_vector contains '\_ClinVarP':  
     * RETURN 'FLAG\_FOR\_REVIEW', 'LOW\_RISK\_IGNORED\_EVIDENCE', ...  
   * **Default Case**:  
     * ELSE: RETURN 'PASS', ...

## **Dependencies**

The node's environment will require the following libraries:

* shap  
* numpy  
* pandas  
* tensorflow or torch (matching the Risk Model)

## **Action Plan**

1. **\[Prerequisite\]** Modify the Feature Vector Builder to produce the required metadata-rich feature names.  
2. **\[Prerequisite\]** Generate and save the K-Means summarized background dataset.  
3. **\[Step 1\]** Implement the core node function (run\_sanity\_check) with placeholders for helper functions.  
4. **\[Step 2\]** Write the get\_top\_contributors helper function to process and sort SHAP values.  
5. **\[Step 3\]** Implement the apply\_sanity\_rules engine with the defined logic.  
6. **\[Step 4\]** Write the format\_output function to build the final JSON object.  
7. **\[Step 5\]** Write unit tests for each sanity rule using mock data to ensure they trigger correctly.