#!/usr/bin/env python3
"""
Analysis of Static Model Outputs for Fusion Layer

This script generates synthetic data representing the 5 static model outputs
and creates visualizations to help determine whether linear or logistic 
regression would be better for the fusion layer.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, roc_auc_score, accuracy_score
from fusion_layer import create_synthetic_training_data, StaticModelInputs
import warnings
warnings.filterwarnings('ignore')

def analyze_static_model_outputs(n_samples=5000):
    """
    Generate and analyze static model outputs to determine best regression approach.
    """
    print("🧬 Analyzing Static Model Outputs for Fusion Layer")
    print("=" * 60)
    
    # Generate synthetic training data
    print(f"📊 Generating {n_samples} synthetic samples...")
    training_data = create_synthetic_training_data(n_samples=n_samples)
    
    # Convert to DataFrame for analysis
    data_rows = []
    for inputs, risk_score in training_data:
        row = inputs.to_dict()
        row['risk_score'] = risk_score
        data_rows.append(row)
    
    df = pd.DataFrame(data_rows)
    
    # Convert categorical to numerical for correlation analysis
    df_numeric = df.copy()
    clinvar_map = {'pathogenic': 3, 'benign': 0, 'uncertain': 1, 'not_found': 0.5}
    df_numeric['clinvar_numeric'] = df_numeric['clinvar_classification'].map(clinvar_map)
    
    print("✅ Data generation complete")
    print(f"   Risk score range: {df['risk_score'].min():.3f} - {df['risk_score'].max():.3f}")
    print(f"   Risk score mean: {df['risk_score'].mean():.3f}")
    
    # Create comprehensive plots
    fig = plt.figure(figsize=(20, 24))
    
    # 1. Distribution of each static model output
    print("\n📈 Creating distribution plots...")
    
    # PRS Score Distribution
    plt.subplot(5, 4, 1)
    plt.hist(df['prs_score'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
    plt.title('PRS Score Distribution')
    plt.xlabel('PRS Score')
    plt.ylabel('Frequency')
    
    # ClinVar Classification Distribution
    plt.subplot(5, 4, 2)
    clinvar_counts = df['clinvar_classification'].value_counts()
    plt.bar(clinvar_counts.index, clinvar_counts.values, color='lightcoral', alpha=0.7)
    plt.title('ClinVar Classification Distribution')
    plt.xlabel('Classification')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    
    # CADD Score Distribution
    plt.subplot(5, 4, 3)
    plt.hist(df['cadd_score'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
    plt.title('CADD Score Distribution')
    plt.xlabel('CADD Score')
    plt.ylabel('Frequency')
    
    # TCGA Enrichment Distribution
    plt.subplot(5, 4, 4)
    plt.hist(df['tcga_enrichment'], bins=30, alpha=0.7, color='orange', edgecolor='black')
    plt.title('TCGA Enrichment Distribution')
    plt.xlabel('TCGA Enrichment')
    plt.ylabel('Frequency')
    
    # Gene Burden Score Distribution
    plt.subplot(5, 4, 5)
    plt.hist(df['gene_burden_score'], bins=15, alpha=0.7, color='purple', edgecolor='black')
    plt.title('Gene Burden Score Distribution')
    plt.xlabel('Gene Burden Score')
    plt.ylabel('Frequency')
    
    # Risk Score Distribution
    plt.subplot(5, 4, 6)
    plt.hist(df['risk_score'], bins=30, alpha=0.7, color='red', edgecolor='black')
    plt.title('Risk Score Distribution (Target)')
    plt.xlabel('Risk Score')
    plt.ylabel('Frequency')
    
    # 2. Correlation analysis
    print("🔍 Analyzing correlations...")
    
    plt.subplot(5, 4, 7)
    numeric_cols = ['prs_score', 'cadd_score', 'tcga_enrichment', 'gene_burden_score', 'clinvar_numeric', 'risk_score']
    corr_matrix = df_numeric[numeric_cols].corr()
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, square=True)
    plt.title('Feature Correlation Matrix')
    
    # 3. Scatter plots against risk score
    print("📊 Creating scatter plots...")
    
    # PRS vs Risk Score
    plt.subplot(5, 4, 8)
    plt.scatter(df['prs_score'], df['risk_score'], alpha=0.5, s=1)
    z = np.polyfit(df['prs_score'], df['risk_score'], 1)
    p = np.poly1d(z)
    plt.plot(df['prs_score'], p(df['prs_score']), "r--", alpha=0.8)
    plt.xlabel('PRS Score')
    plt.ylabel('Risk Score')
    plt.title('PRS vs Risk Score')
    
    # CADD vs Risk Score
    plt.subplot(5, 4, 9)
    plt.scatter(df['cadd_score'], df['risk_score'], alpha=0.5, s=1)
    z = np.polyfit(df['cadd_score'], df['risk_score'], 1)
    p = np.poly1d(z)
    plt.plot(df['cadd_score'], p(df['cadd_score']), "r--", alpha=0.8)
    plt.xlabel('CADD Score')
    plt.ylabel('Risk Score')
    plt.title('CADD vs Risk Score')
    
    # TCGA vs Risk Score
    plt.subplot(5, 4, 10)
    plt.scatter(df['tcga_enrichment'], df['risk_score'], alpha=0.5, s=1)
    z = np.polyfit(df['tcga_enrichment'], df['risk_score'], 1)
    p = np.poly1d(z)
    plt.plot(df['tcga_enrichment'], p(df['tcga_enrichment']), "r--", alpha=0.8)
    plt.xlabel('TCGA Enrichment')
    plt.ylabel('Risk Score')
    plt.title('TCGA vs Risk Score')
    
    # Gene Burden vs Risk Score
    plt.subplot(5, 4, 11)
    plt.scatter(df['gene_burden_score'], df['risk_score'], alpha=0.5, s=1)
    z = np.polyfit(df['gene_burden_score'], df['risk_score'], 1)
    p = np.poly1d(z)
    plt.plot(df['gene_burden_score'], p(df['gene_burden_score']), "r--", alpha=0.8)
    plt.xlabel('Gene Burden Score')
    plt.ylabel('Risk Score')
    plt.title('Gene Burden vs Risk Score')
    
    # 4. ClinVar boxplot
    plt.subplot(5, 4, 12)
    df.boxplot(column='risk_score', by='clinvar_classification', ax=plt.gca())
    plt.title('Risk Score by ClinVar Classification')
    plt.suptitle('')
    plt.xticks(rotation=45)
    
    # 5. Model comparison - prepare data
    print("🔬 Comparing Linear vs Logistic Regression...")
    
    # One-hot encode ClinVar for modeling
    clinvar_dummies = pd.get_dummies(df['clinvar_classification'], prefix='clinvar')
    X = pd.concat([
        df[['prs_score', 'cadd_score', 'tcga_enrichment', 'gene_burden_score']],
        clinvar_dummies
    ], axis=1)
    
    # Prepare targets for both regression types
    y_continuous = df['risk_score']  # For linear regression
    y_binary = (df['risk_score'] > 0.5).astype(int)  # For logistic regression
    
    # Split data
    X_train, X_test, y_cont_train, y_cont_test = train_test_split(X, y_continuous, test_size=0.2, random_state=42)
    _, _, y_bin_train, y_bin_test = train_test_split(X, y_binary, test_size=0.2, random_state=42)
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Linear Regression
    linear_reg = LinearRegression()
    linear_reg.fit(X_train_scaled, y_cont_train)
    y_pred_linear = linear_reg.predict(X_test_scaled)
    
    linear_mse = mean_squared_error(y_cont_test, y_pred_linear)
    linear_r2 = r2_score(y_cont_test, y_pred_linear)
    
    # Logistic Regression
    logistic_reg = LogisticRegression(random_state=42, max_iter=1000)
    logistic_reg.fit(X_train_scaled, y_bin_train)
    y_pred_logistic = logistic_reg.predict(X_test_scaled)
    y_pred_proba = logistic_reg.predict_proba(X_test_scaled)[:, 1]
    
    logistic_acc = accuracy_score(y_bin_test, y_pred_logistic)
    logistic_auc = roc_auc_score(y_bin_test, y_pred_proba)
    
    # 6. Prediction plots
    plt.subplot(5, 4, 13)
    plt.scatter(y_cont_test, y_pred_linear, alpha=0.5)
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlabel('Actual Risk Score')
    plt.ylabel('Predicted Risk Score')
    plt.title(f'Linear Regression\nMSE: {linear_mse:.4f}, R²: {linear_r2:.4f}')
    
    plt.subplot(5, 4, 14)
    plt.scatter(y_bin_test, y_pred_proba, alpha=0.5)
    plt.plot([0, 1], [0, 1], 'r--')
    plt.xlabel('Actual Binary Risk (>0.5)')
    plt.ylabel('Predicted Probability')
    plt.title(f'Logistic Regression\nACC: {logistic_acc:.4f}, AUC: {logistic_auc:.4f}')
    
    # 7. Residual analysis for linear regression
    plt.subplot(5, 4, 15)
    residuals = y_cont_test - y_pred_linear
    plt.scatter(y_pred_linear, residuals, alpha=0.5)
    plt.axhline(y=0, color='r', linestyle='--')
    plt.xlabel('Predicted Values')
    plt.ylabel('Residuals')
    plt.title('Linear Regression Residuals')
    
    # 8. Feature importance comparison
    plt.subplot(5, 4, 16)
    feature_names = X.columns
    linear_importance = np.abs(linear_reg.coef_)
    logistic_importance = np.abs(logistic_reg.coef_[0])
    
    x_pos = np.arange(len(feature_names))
    width = 0.35
    
    plt.bar(x_pos - width/2, linear_importance, width, label='Linear', alpha=0.7)
    plt.bar(x_pos + width/2, logistic_importance, width, label='Logistic', alpha=0.7)
    plt.xlabel('Features')
    plt.ylabel('Absolute Coefficient')
    plt.title('Feature Importance Comparison')
    plt.xticks(x_pos, feature_names, rotation=45, ha='right')
    plt.legend()
    
    # 9. Risk score distribution analysis
    plt.subplot(5, 4, 17)
    plt.hist(df['risk_score'], bins=50, alpha=0.7, density=True, label='Actual')
    plt.axvline(df['risk_score'].mean(), color='red', linestyle='--', label=f'Mean: {df["risk_score"].mean():.3f}')
    plt.axvline(df['risk_score'].median(), color='green', linestyle='--', label=f'Median: {df["risk_score"].median():.3f}')
    plt.xlabel('Risk Score')
    plt.ylabel('Density')
    plt.title('Risk Score Distribution Analysis')
    plt.legend()
    
    # 10. Risk thresholds analysis
    plt.subplot(5, 4, 18)
    thresholds = [0.25, 0.5, 0.75]
    threshold_names = ['Low', 'Moderate', 'High', 'Very High']
    
    # Calculate counts for each category
    low_count = (df['risk_score'] <= 0.25).sum()
    moderate_count = ((df['risk_score'] > 0.25) & (df['risk_score'] <= 0.5)).sum()
    high_count = ((df['risk_score'] > 0.5) & (df['risk_score'] <= 0.75)).sum()
    very_high_count = (df['risk_score'] > 0.75).sum()
    
    counts = [low_count, moderate_count, high_count, very_high_count]
    
    plt.bar(threshold_names, counts, alpha=0.7, color=['green', 'yellow', 'orange', 'red'])
    plt.xlabel('Risk Category')
    plt.ylabel('Count')
    plt.title('Risk Category Distribution')
    
    # 11. ClinVar impact analysis
    plt.subplot(5, 4, 19)
    clinvar_means = df.groupby('clinvar_classification')['risk_score'].mean()
    plt.bar(clinvar_means.index, clinvar_means.values, alpha=0.7)
    plt.xlabel('ClinVar Classification')
    plt.ylabel('Mean Risk Score')
    plt.title('Mean Risk Score by ClinVar')
    plt.xticks(rotation=45)
    
    # 12. Non-linearity test
    plt.subplot(5, 4, 20)
    # Test for non-linearity by comparing polynomial fits
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.pipeline import Pipeline
    
    # Create polynomial features for PRS (most important continuous feature)
    degrees = [1, 2, 3]
    colors = ['blue', 'red', 'green']
    
    for degree, color in zip(degrees, colors):
        poly_features = PolynomialFeatures(degree=degree)
        X_poly = poly_features.fit_transform(df[['prs_score']])
        poly_reg = LinearRegression()
        poly_reg.fit(X_poly, df['risk_score'])
        
        # Create smooth curve for plotting
        prs_range = np.linspace(df['prs_score'].min(), df['prs_score'].max(), 100).reshape(-1, 1)
        prs_poly = poly_features.transform(prs_range)
        risk_pred = poly_reg.predict(prs_poly)
        
        plt.plot(prs_range, risk_pred, color=color, label=f'Degree {degree}', linewidth=2)
    
    plt.scatter(df['prs_score'], df['risk_score'], alpha=0.3, s=1, color='gray')
    plt.xlabel('PRS Score')
    plt.ylabel('Risk Score')
    plt.title('Polynomial Fit Analysis (PRS)')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('static_model_analysis.png', dpi=300, bbox_inches='tight')
    print(f"📊 Comprehensive analysis plot saved as 'static_model_analysis.png'")
    
    # Print detailed analysis results
    print("\n" + "="*60)
    print("📋 ANALYSIS SUMMARY")
    print("="*60)
    
    print("\n🔢 Data Statistics:")
    print(f"   Total samples: {len(df):,}")
    print(f"   Risk score mean: {df['risk_score'].mean():.3f}")
    print(f"   Risk score std: {df['risk_score'].std():.3f}")
    print(f"   Risk score range: {df['risk_score'].min():.3f} - {df['risk_score'].max():.3f}")
    
    print("\n📊 Feature Correlations with Risk Score:")
    correlations = df_numeric[['prs_score', 'cadd_score', 'tcga_enrichment', 'gene_burden_score', 'clinvar_numeric']].corrwith(df_numeric['risk_score'])
    for feature, corr in correlations.items():
        print(f"   {feature}: {corr:.3f}")
    
    print("\n🎯 Model Performance Comparison:")
    print(f"   Linear Regression:")
    print(f"     - MSE: {linear_mse:.4f}")
    print(f"     - R²: {linear_r2:.4f}")
    print(f"   Logistic Regression (binary classification):")
    print(f"     - Accuracy: {logistic_acc:.4f}")
    print(f"     - AUC: {logistic_auc:.4f}")
    
    print("\n📈 Risk Category Distribution:")
    risk_categories = pd.cut(df['risk_score'], bins=[0, 0.25, 0.5, 0.75, 1.0], labels=['Low', 'Moderate', 'High', 'Very High'])
    category_counts = risk_categories.value_counts()
    for category, count in category_counts.items():
        pct = count / len(df) * 100
        print(f"   {category}: {count} ({pct:.1f}%)")
    
    print("\n🔍 Recommendation:")
    
    # Decision logic
    if linear_r2 > 0.7 and linear_mse < 0.05:
        if df['risk_score'].std() > 0.2:  # Good spread in continuous values
            recommendation = "LINEAR REGRESSION"
            reason = f"Strong linear relationship (R²={linear_r2:.3f}) with good continuous prediction capability"
        else:
            recommendation = "LOGISTIC REGRESSION"
            reason = "Limited continuous variation suggests binary classification is more appropriate"
    elif logistic_auc > 0.85:
        recommendation = "LOGISTIC REGRESSION"
        reason = f"Excellent binary classification performance (AUC={logistic_auc:.3f})"
    else:
        recommendation = "GRADIENT BOOSTING (current choice)"
        reason = "Complex non-linear relationships suggest tree-based methods are optimal"
    
    print(f"   Recommended approach: {recommendation}")
    print(f"   Reasoning: {reason}")
    
    print("\n💡 Additional Insights:")
    print(f"   - ClinVar 'pathogenic' variants have {df[df['clinvar_classification']=='pathogenic']['risk_score'].mean():.3f} mean risk")
    print(f"   - CADD scores above 20 correlate with {df[df['cadd_score']>20]['risk_score'].mean():.3f} mean risk")
    print(f"   - Data shows {'linear' if max(correlations) > 0.6 else 'non-linear'} relationships")
    
    return df, {
        'linear_mse': linear_mse,
        'linear_r2': linear_r2,
        'logistic_acc': logistic_acc,
        'logistic_auc': logistic_auc,
        'recommendation': recommendation
    }

if __name__ == "__main__":
    # Run the analysis
    df, results = analyze_static_model_outputs(n_samples=5000)
    
    print(f"\n✅ Analysis complete! Check 'static_model_analysis.png' for detailed visualizations.") 