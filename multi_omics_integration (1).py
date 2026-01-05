"""
Multi-Omics Data Integration Pipeline
Integrates genomics, transcriptomics, and proteomics data
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns
import matplotlib.pyplot as plt

class MultiOmicsIntegrator:
    """
    Integrate multiple omics layers for systems biology analysis
    """
    
    def __init__(self):
        self.genomics_data = None
        self.transcriptomics_data = None
        self.proteomics_data = None
        self.integrated_data = None
        
    def load_genomics(self, filepath: str) -> pd.DataFrame:
        """Load and preprocess genomics data (variant calls)"""
        df = pd.read_csv(filepath, sep='\t')
        # Filter for high-quality variants
        df_filtered = df[df['QUAL'] > 30]
        self.genomics_data = df_filtered
        return df_filtered
    
    def load_transcriptomics(self, filepath: str) -> pd.DataFrame:
        """Load and normalize RNA-seq data"""
        df = pd.read_csv(filepath, index_col=0)
        # Log2 transformation and normalization
        df_log = np.log2(df + 1)
        self.transcriptomics_data = df_log
        return df_log
    
    def load_proteomics(self, filepath: str) -> pd.DataFrame:
        """Load and normalize proteomics data"""
        df = pd.read_csv(filepath, index_col=0)
        scaler = StandardScaler()
        df_scaled = pd.DataFrame(
            scaler.fit_transform(df),
            index=df.index,
            columns=df.columns
        )
        self.proteomics_data = df_scaled
        return df_scaled
    
    def integrate_layers(self) -> pd.DataFrame:
        """
        Integrate multiple omics layers
        Returns combined feature matrix
        """
        # Merge on common samples
        merged = pd.concat([
            self.transcriptomics_data,
            self.proteomics_data
        ], axis=1, join='inner')
        
        self.integrated_data = merged
        return merged
    
    def dimensionality_reduction(self, n_components: int = 2):
        """
        PCA for visualization and feature reduction
        """
        pca = PCA(n_components=n_components)
        transformed = pca.fit_transform(self.integrated_data)
        
        explained_var = pca.explained_variance_ratio_
        print(f"Explained variance: {explained_var}")
        
        return transformed, pca
    
    def identify_biomarkers(self, threshold: float = 2.0) -> pd.DataFrame:
        """
        Identify potential biomarkers based on variance
        """
        variances = self.integrated_data.var()
        high_var_features = variances[variances > threshold]
        
        biomarkers = self.integrated_data[high_var_features.index]
        return biomarkers.sort_values(by=biomarkers.columns[0], ascending=False)

def main():
    """Example workflow"""
    integrator = MultiOmicsIntegrator()
    
    # Load data layers
    integrator.load_genomics('data/variants.vcf')
    integrator.load_transcriptomics('data/rnaseq_counts.csv')
    integrator.load_proteomics('data/proteomics.csv')
    
    # Integrate
    integrated = integrator.integrate_layers()
    print(f"Integrated data shape: {integrated.shape}")
    
    # Dimensionality reduction
    transformed, pca = integrator.dimensionality_reduction()
    
    # Find biomarkers
    biomarkers = integrator.identify_biomarkers()
    biomarkers.to_csv('results/candidate_biomarkers.csv')
    
    print(f"Candidate biomarkers identified: {len(biomarkers)}")

if __name__ == "__main__":
    main()
