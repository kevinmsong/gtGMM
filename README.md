# gtGMM v2.0: Topological Data Analysis for Gene Networks

**Unsupervised discovery of functional gene modules through network topology**

`gtGMM` is a Python library for unsupervised topological analysis of gene expression networks. It transforms gene networks into continuous topological landscapes, enabling quantitative analysis of network organization using geometric and topological invariants. By integrating protein-protein interaction networks with expression data, `gtGMM` reveals functional architecture without requiring prior pathway knowledge or gene annotations.

This major v2.0 release includes a complete refactoring of the original codebase, offering a cleaner API, enhanced performance, and new features for comprehensive, publication-quality analysis.

**Authors:** Kevin Song, John Zhang, Lei Ye MD PhD, Jianyi Zhang MD PhD  
**Institution:** Department of Biomedical Engineering, The University of Alabama at Birmingham  
**Contact:** kmsong@uab.edu

---

## Key Features in v2.0

*   **High-Resolution Terrains:** Default `500x500` resolution for detailed topological analysis.
*   **Adaptive GMM Optimization:** Automatic selection of the optimal number of GMM components using Bayesian Information Criterion (BIC).
*   **Integrated TDA:** Seamlessly compute topological data analysis (TDA) metrics on `GeneTerrain` objects.
*   **Pathway Enrichment:** A new `EnrichmentAnalyzer` for robust pathway enrichment analysis with FDR correction.
*   **Streamlined API:** A clean, high-level API with convenience functions for common workflows.
*   **Publication-Quality Visualizations:** Enhanced plotting functions for creating publication-ready figures.

---

## Installation

```bash
git clone https://github.com/kevinmsong/gtGMM.git
cd gtGMM
pip install -e .
```

**Requirements:** Python 3.8+, NumPy, Pandas, scikit-learn, SciPy, Matplotlib, NetworkX, requests, gseapy

---

## Quick Start

This example demonstrates a complete analysis workflow for the GSE130699 cardiac regeneration dataset.

```python
import gtgmm

# 1. Load pre-processed cardiac regeneration data
genes, p1_expr, p8_expr = gtgmm.data.load_gse130699()

# 2. Create and build the P1 GeneTerrain
# This automatically fetches interactions, computes layout, and optimizes sigma
gt_p1 = gtgmm.create_terrain(genes, p1_expr, species=10090)

# 3. Run GMM decomposition and TDA
gmm_p1 = gt_p1.run_gmm(max_components=8, min_components=3)
tda_p1 = gt_p1.run_tda()

# 4. Perform pathway enrichment analysis
pathways = gtgmm.data.get_cardiac_pathways()
enrichment_analyzer = gtgmm.EnrichmentAnalyzer(pathways, background_genes=genes)

components_p1 = gmm_p1.extract_gene_components(gt_p1)
enrichment_results = enrichment_analyzer.analyze_components(components_p1, condition="P1")

# 5. Visualize the results
vis = gtgmm.TerrainVisualizer()
vis.plot_3d_terrain(gt_p1, title="P1 GeneTerrain")
vis.plot_component_map(gmm_p1, gt_p1, title="P1 GMM Components")

print("P1 TDA Summary:")
print(tda_p1.get_tda_summary())

print("\nTop 5 Enriched Pathways in P1:")
print(enrichment_results.head())
```

---

## Documentation

Full documentation, including API references and tutorials, is available in the `docs/` directory.

---

## Citing gtGMM

If you use `gtGMM` in your research, please cite:

```bibtex
@software{gtgmm2024,
  author = {Song, Kevin and Zhang, John and Ye, Lei and Zhang, Jianyi},
  title = {gtGMM: Topological Data Analysis for Gene Expression Networks},
  year = {2024},
  institution = {Department of Biomedical Engineering, The University of Alabama at Birmingham},
  url = {https://github.com/kevinmsong/gtGMM}
}
```

---

## License

MIT License - see the `LICENSE` file for details.

