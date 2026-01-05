# Multi-Omics Integration Pipeline

**Scalable Bioinformatics Pipeline for Integrating Multi-Modal Omics Data**

Computational framework for integrating genomics, transcriptomics, proteomics, and epigenomics data to identify biomarkers and understand complex biological systems.

## Overview

This repository contains production-ready pipelines for multi-omics data integration, designed for systems biology research, biomarker discovery, and precision medicine applications. The framework supports multiple data types and provides comprehensive tools for quality control, normalization, integration, and visualization.

## Key Features

✅ **Multi-modal data integration** - Genomics, transcriptomics, proteomics, metabolomics, epigenomics  
✅ **Automated QC pipelines** - Quality control and filtering for each omics layer  
✅ **Advanced normalization** - TMM, DESeq2, quantile normalization, batch correction  
✅ **Dimensionality reduction** - PCA, t-SNE, UMAP for visualization  
✅ **Biomarker discovery** - Statistical testing and machine learning approaches  
✅ **Pathway enrichment** - GO, KEGG, Reactome pathway analysis  
✅ **Interactive visualizations** - Publication-ready plots and dashboards  
✅ **Scalable architecture** - Handles datasets from 10s to 1000s of samples  

## Use Cases

### 1. Cancer Research
- Integrate tumor genomics, transcriptomics, and proteomics
- Identify driver mutations and their functional effects
- Discover therapeutic targets

### 2. Neurodegenerative Disease
- Multi-omics profiling of Alzheimer's, Parkinson's patients
- Identify disease biomarkers across molecular layers
- Understand disease mechanisms

### 3. Drug Discovery
- Predict drug response from multi-omics profiles
- Identify resistance mechanisms
- Stratify patient populations

### 4. Systems Biology
- Reconstruct regulatory networks
- Identify key molecular players
- Understand emergent properties

## Repository Structure

```
Multi-Omics-Pipeline/
├── data/
│   ├── raw/
│   │   ├── genomics/
│   │   │   ├── vcf_files/              # Variant call format files
│   │   │   └── metadata.csv
│   │   ├── transcriptomics/
│   │   │   ├── rnaseq_counts.csv       # RNA-seq count matrix
│   │   │   └── sample_info.csv
│   │   ├── proteomics/
│   │   │   ├── protein_abundance.csv
│   │   │   └── peptide_data.csv
│   │   ├── epigenomics/
│   │   │   ├── atac_seq/               # ATAC-seq peaks
│   │   │   └── chip_seq/               # ChIP-seq data
│   │   └── metabolomics/
│   │       └── metabolite_levels.csv
│   ├── processed/
│   │   ├── normalized_rnaseq.csv
│   │   ├── filtered_variants.csv
│   │   └── standardized_proteins.csv
│   └── integrated/
│       ├── combined_matrix.csv
│       └── integration_metadata.json
├── pipelines/
│   ├── 01_genomics_processing.py
│   ├── 02_transcriptomics_processing.py
│   ├── 03_proteomics_processing.py
│   ├── 04_epigenomics_processing.py
│   ├── 05_data_integration.py
│   ├── 06_biomarker_discovery.py
│   └── 07_pathway_analysis.py
├── modules/
│   ├── qc/
│   │   ├── genomics_qc.py
│   │   ├── rnaseq_qc.py
│   │   └── proteomics_qc.py
│   ├── normalization/
│   │   ├── rna_normalization.py
│   │   ├── protein_normalization.py
│   │   └── batch_correction.py
│   ├── integration/
│   │   ├── mofa.py                     # Multi-Omics Factor Analysis
│   │   ├── snf.py                      # Similarity Network Fusion
│   │   └── concatenation.py
│   ├── analysis/
│   │   ├── differential_expression.py
│   │   ├── variant_annotation.py
│   │   └── correlation_analysis.py
│   └── visualization/
│       ├── heatmaps.py
│       ├── volcano_plots.py
│       ├── network_plots.py
│       └── interactive_dash.py
├── workflows/
│   ├── cancer_analysis.py
│   ├── neurodegeneration_pipeline.py
│   └── biomarker_discovery.py
├── results/
│   ├── figures/
│   ├── tables/
│   └── reports/
├── tests/
│   └── test_pipelines.py
├── config/
│   ├── pipeline_config.yaml
│   └── analysis_parameters.json
├── docs/
│   ├── installation.md
│   ├── tutorials/
│   └── api_reference.md
├── requirements.txt
├── environment.yml
└── README.md
```

## Installation

### Prerequisites

```bash
# Python 3.8+
python --version

# Conda (recommended for environment management)
conda --version
```

### Quick Start

```bash
# Clone repository
git clone https://github.com/DASPAPIYA28/Multi-Omics-Pipeline.git
cd Multi-Omics-Pipeline

# Create conda environment
conda env create -f environment.yml
conda activate multiomics

# Or use pip
pip install -r requirements.txt

# Install package in development mode
pip install -e .
```

### Dependencies

**Core:**
- Python 3.8+
- pandas >= 1.3.0
- numpy >= 1.21.0
- scipy >= 1.7.0

**Bioinformatics:**
- biopython >= 1.79
- pysam >= 0.16.0
- pyvcf3 >= 1.0.3
- pyranges >= 0.0.111

**Statistics & Machine Learning:**
- scikit-learn >= 0.24.0
- statsmodels >= 0.13.0
- lifelines >= 0.26.0

**Visualization:**
- matplotlib >= 3.4.0
- seaborn >= 0.11.0
- plotly >= 5.0.0
- dash >= 2.0.0

**Specialized:**
- DESeq2 (via rpy2)
- edgeR (via rpy2)
- limma (via rpy2)

## Usage

### 1. Quick Example: Basic Integration

```python
from modules.integration import MultiOmicsIntegrator
from modules.visualization import IntegrationVisualizer

# Initialize integrator
integrator = MultiOmicsIntegrator()

# Load data layers
genomics = integrator.load_genomics('data/raw/genomics/variants.vcf')
rnaseq = integrator.load_transcriptomics('data/raw/transcriptomics/counts.csv')
proteins = integrator.load_proteomics('data/raw/proteomics/abundance.csv')

# Quality control
genomics_qc = integrator.qc_genomics(genomics, min_quality=30)
rnaseq_qc = integrator.qc_rnaseq(rnaseq, min_counts=10)
proteins_qc = integrator.qc_proteomics(proteins, missing_threshold=0.3)

# Normalize
rnaseq_norm = integrator.normalize_rnaseq(rnaseq_qc, method='DESeq2')
proteins_norm = integrator.normalize_proteomics(proteins_qc, method='quantile')

# Integrate
integrated_data = integrator.integrate([
    genomics_qc,
    rnaseq_norm,
    proteins_norm
], method='concatenation')

# Visualize
viz = IntegrationVisualizer()
viz.plot_pca(integrated_data, color_by='sample_type')
viz.plot_heatmap(integrated_data, cluster=True)
```

### 2. RNA-seq Processing Pipeline

```python
from pipelines import RNAseqProcessor

# Initialize processor
processor = RNAseqProcessor(
    count_file='data/raw/transcriptomics/counts.csv',
    metadata='data/raw/transcriptomics/metadata.csv'
)

# Quality control
processor.filter_low_counts(min_count=10, min_samples=3)
processor.detect_outliers(method='PCA', threshold=3)

# Normalization
processor.normalize(method='TMM')  # or 'DESeq2', 'quantile'

# Differential expression
de_results = processor.differential_expression(
    design='~ condition',
    contrast=['treatment', 'control']
)

# Save results
processor.save_processed_data('data/processed/rnaseq_normalized.csv')
de_results.to_csv('results/tables/de_genes.csv')
```

### 3. Variant Analysis Pipeline

```python
from pipelines import GenomicsProcessor

# Initialize
processor = GenomicsProcessor()

# Load VCF
variants = processor.load_vcf('data/raw/genomics/variants.vcf')

# Filter variants
filtered = processor.filter_variants(
    min_quality=30,
    min_depth=10,
    max_missing=0.1
)

# Annotate
annotated = processor.annotate_variants(
    reference='GRCh38',
    databases=['dbSNP', 'ClinVar', 'gnomAD']
)

# Identify high-impact variants
high_impact = processor.filter_high_impact(
    impact_levels=['HIGH', 'MODERATE']
)

# Save
annotated.to_csv('data/processed/annotated_variants.csv')
```

### 4. Proteomics Processing

```python
from pipelines import ProteomicsProcessor

# Initialize
processor = ProteomicsProcessor(
    abundance_file='data/raw/proteomics/protein_abundance.csv'
)

# Impute missing values
processor.impute_missing(method='MinProb')  # or 'KNN', 'SVD'

# Normalize
processor.normalize(method='quantile')  # or 'median', 'TMM'

# Batch correction
processor.correct_batch_effects(batch_column='batch')

# Differential abundance
da_results = processor.differential_abundance(
    design='~ condition',
    contrast=['disease', 'healthy']
)

# Save
processor.save_processed_data('data/processed/proteomics_normalized.csv')
```

### 5. Complete Integration Workflow

```python
from workflows import MultiOmicsWorkflow

# Initialize workflow
workflow = MultiOmicsWorkflow(
    project_name='cancer_study',
    output_dir='results/cancer_analysis'
)

# Add data layers
workflow.add_layer('genomics', 'data/raw/genomics/variants.vcf')
workflow.add_layer('transcriptomics', 'data/raw/transcriptomics/counts.csv')
workflow.add_layer('proteomics', 'data/raw/proteomics/abundance.csv')
workflow.add_layer('epigenomics', 'data/raw/epigenomics/atac_peaks.bed')

# Configure pipeline
workflow.configure(
    qc_params={'min_quality': 30, 'min_counts': 10},
    normalization={'rnaseq': 'DESeq2', 'proteomics': 'quantile'},
    integration_method='MOFA'  # or 'SNF', 'concatenation'
)

# Run complete analysis
results = workflow.run(
    identify_biomarkers=True,
    pathway_analysis=True,
    create_report=True
)

# Generate report
workflow.generate_report(
    output_file='results/cancer_analysis/report.html',
    include_figures=True
)
```

### 6. Biomarker Discovery

```python
from modules.analysis import BiomarkerDiscovery

# Initialize
discovery = BiomarkerDiscovery(
    integrated_data='data/integrated/combined_matrix.csv',
    metadata='data/integrated/sample_metadata.csv'
)

# Feature selection
selected_features = discovery.select_features(
    method='random_forest',  # or 'lasso', 'elastic_net'
    n_features=100,
    cv_folds=5
)

# Build classifier
model = discovery.build_classifier(
    features=selected_features,
    algorithm='random_forest',
    optimize_hyperparameters=True
)

# Evaluate
performance = discovery.evaluate_model(
    model=model,
    test_set='holdout',
    metrics=['accuracy', 'AUC', 'sensitivity', 'specificity']
)

# Identify top biomarkers
top_biomarkers = discovery.rank_features(
    method='importance',
    n_top=20
)

print(f"Model AUC: {performance['AUC']:.3f}")
print(f"Top biomarkers: {top_biomarkers['feature'].tolist()}")
```

## Advanced Features

### Multi-Omics Factor Analysis (MOFA)

```python
from modules.integration import MOFAIntegration

# Initialize MOFA
mofa = MOFAIntegration(n_factors=10)

# Prepare data views
views = {
    'genomics': genomics_data,
    'transcriptomics': rnaseq_data,
    'proteomics': protein_data
}

# Run MOFA
model = mofa.fit(views)

# Extract factors
factors = mofa.get_factors()

# Analyze factor weights
weights = mofa.get_weights()

# Visualize
mofa.plot_variance_explained()
mofa.plot_factor_correlations()
mofa.plot_top_features(factor=1, n_features=20)
```

### Pathway Enrichment Analysis

```python
from modules.analysis import PathwayEnrichment

# Initialize
enrichment = PathwayEnrichment()

# Load differentially expressed genes
de_genes = pd.read_csv('results/tables/de_genes.csv')
significant_genes = de_genes[de_genes['padj'] < 0.05]['gene_id'].tolist()

# Run enrichment
go_results = enrichment.run_go_enrichment(
    gene_list=significant_genes,
    background='genome',
    ontology='BP'  # Biological Process
)

kegg_results = enrichment.run_kegg_enrichment(
    gene_list=significant_genes,
    organism='hsa'  # Homo sapiens
)

# Visualize
enrichment.plot_dotplot(go_results, top_n=20)
enrichment.plot_network(kegg_results)
```

### Network Analysis

```python
from modules.analysis import NetworkAnalysis

# Initialize
network = NetworkAnalysis()

# Build correlation network
correlation_net = network.build_correlation_network(
    data=integrated_data,
    method='spearman',
    threshold=0.7
)

# Identify modules
modules = network.detect_communities(
    network=correlation_net,
    algorithm='louvain'
)

# Find hub genes
hubs = network.identify_hubs(
    network=correlation_net,
    method='degree',
    top_n=50
)

# Visualize
network.plot_network(
    correlation_net,
    color_by='module',
    size_by='degree'
)
```

## Configuration

### Pipeline Configuration (YAML)

```yaml
# config/pipeline_config.yaml

project:
  name: "my_multiomics_study"
  output_dir: "results/"

data:
  genomics:
    file: "data/raw/genomics/variants.vcf"
    quality_threshold: 30
    depth_threshold: 10
  
  transcriptomics:
    file: "data/raw/transcriptomics/counts.csv"
    min_counts: 10
    normalization: "DESeq2"
  
  proteomics:
    file: "data/raw/proteomics/abundance.csv"
    missing_threshold: 0.3
    normalization: "quantile"

integration:
  method: "MOFA"  # or "SNF", "concatenation"
  n_factors: 10
  
analysis:
  differential_analysis: true
  pathway_enrichment: true
  biomarker_discovery:
    enabled: true
    n_features: 100
    classifier: "random_forest"
  
visualization:
  create_plots: true
  interactive_dashboard: true
```

## Example Workflows

### Cancer Subtype Classification

```python
# Complete workflow for cancer subtype analysis
from workflows import CancerAnalysis

analysis = CancerAnalysis(
    genomics='data/raw/genomics/tumor_variants.vcf',
    rnaseq='data/raw/transcriptomics/tumor_expression.csv',
    proteomics='data/raw/proteomics/tumor_proteins.csv',
    metadata='data/raw/metadata/clinical_data.csv'
)

# Run analysis
results = analysis.run(
    identify_subtypes=True,
    n_clusters=4,
    survival_analysis=True,
    biomarker_discovery=True
)

# Generate report
analysis.create_report('results/cancer_subtypes_report.html')
```

### Neurodegenerative Disease Biomarkers

```python
from workflows import NeurodegenerationPipeline

# Initialize
pipeline = NeurodegenerationPipeline(
    disease='Alzheimers',
    stages=['early', 'moderate', 'advanced']
)

# Load data
pipeline.load_multiomics_data(
    genomics='data/genomics/AD_variants.vcf',
    rnaseq='data/rnaseq/brain_expression.csv',
    proteomics='data/proteomics/csf_proteins.csv',
    metabolomics='data/metabolomics/plasma_metabolites.csv'
)

# Run pipeline
results = pipeline.run_analysis(
    stage_comparison=True,
    longitudinal_analysis=True,
    biomarker_validation=True
)

# Identify progression markers
progression_markers = pipeline.identify_progression_biomarkers(
    baseline='healthy',
    follow_up_stages=['early', 'moderate', 'advanced']
)
```

## Performance & Scalability

### Dataset Sizes Tested

| Omics Layer | Samples | Features | Processing Time |
|-------------|---------|----------|-----------------|
| Genomics | 1,000 | 5M variants | ~2 hours |
| RNA-seq | 500 | 20K genes | ~15 minutes |
| Proteomics | 200 | 5K proteins | ~5 minutes |
| Metabolomics | 300 | 500 metabolites | ~2 minutes |

### Optimization Tips

```python
# Enable parallel processing
from modules.utils import enable_parallel_processing
enable_parallel_processing(n_cores=8)

# Use chunking for large datasets
processor.process_in_chunks(
    chunk_size=1000,
    n_jobs=4
)

# Memory-efficient processing
processor.use_sparse_matrices(True)
processor.enable_lazy_loading(True)
```

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{das2024multiomics,
  author = {Das, Papiya},
  title = {Multi-Omics Integration Pipeline: Scalable Framework for Systems Biology},
  year = {2024},
  url = {https://github.com/DASPAPIYA28/Multi-Omics-Pipeline}
}
```

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file.

## Contact

**Papiya Das**  
MSc Biotechnology, University of Leeds  
Email: papiya2810das@gmail.com  
LinkedIn: [linkedin.com/in/p--d](https://linkedin.com/in/p--d)  
GitHub: [@DASPAPIYA28](https://github.com/DASPAPIYA28)

## Acknowledgments

- Bioinformatics community for open-source tools
- University of Leeds for computational resources
- Contributors and users of this pipeline

---

**Note:** This repository is actively maintained. For questions, issues, or feature requests, please open an issue on GitHub.
