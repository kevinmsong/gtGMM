# gtGMM: Topological Data Analysis for Gene Networks

**Unsupervised discovery of functional gene modules through network topology**

`gtGMM` is a Python library for unsupervised topological analysis of gene expression networks. It transforms gene networks into continuous topological landscapes, enabling quantitative analysis of network organization using geometric and topological invariants. In short, genes are placed on a 2D map by their protein–protein interactions, expression becomes height, and the resulting "terrain" is split into modules with a Gaussian mixture model (GMM) and summarized with topological data analysis (TDA). By integrating protein-protein interaction networks with expression data, `gtGMM` reveals functional architecture without requiring prior pathway knowledge or gene annotations.

This major v2.0 release includes a complete refactoring of the original codebase, offering a cleaner API, enhanced performance, and new features for comprehensive, publication-quality analysis.

---

## Key Features in v2.0

*   **High-Resolution Terrains:** Default `500x500` resolution for detailed topological analysis.
*   **Adaptive GMM Optimization:** Automatic selection of the number of GMM components by silhouette score (BIC is also recorded for each candidate).
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

**Requirements:** Python 3.8+, NumPy, pandas, scikit-learn, SciPy, Matplotlib, NetworkX, requests, gseapy (installed automatically by `pip install -e .`).

---

## Quick Start

This example builds a terrain for a small cardiac gene set and runs GMM decomposition, TDA, and pathway enrichment. It uses only functions that exist in `gtgmm` v2.0.

```python
import gtgmm
from gtgmm.gse130699_data import load_gse130699_expression

# 1. Pick genes: the first 5 genes from each built-in cardiac pathway
pathways = gtgmm.data.get_cardiac_pathways()
genes = sorted({g for members in pathways.values() for g in members[:5]})

# 2. Expression values for P1 (regenerative) neonatal mouse hearts.
#    Pass data_dir=... to read GSE130699_P1_expression.csv; without it, the function
#    generates illustrative values based on published GSE130699 results (not real data).
p1_expr = load_gse130699_expression(genes, condition="P1")

# 3. Build the GeneTerrain (fetches STRING interactions, computes layout, optimizes sigma)
gt_p1 = gtgmm.create_terrain(genes, p1_expr, species=10090)  # 10090 = mouse, 9606 = human

# 4. Adaptive GMM: choose the number of components by silhouette score, then fit
gmm_p1 = gtgmm.GMMOptimizer(max_components=8)
terrain_data = gt_p1.get_terrain_data()
k = gmm_p1.optimize(terrain_data, min_components=3)
gmm_p1.fit(terrain_data, n_components=k)
components_p1 = gmm_p1.extract_gene_components(gt_p1)

# 5. Topological summary (curvature, Betti numbers, Euler characteristic, ...)
tda_p1 = gtgmm.TerrainTDA(gt_p1)
print(tda_p1.get_tda_summary())

# 6. Pathway enrichment per GMM component (FDR-corrected)
analyzer = gtgmm.EnrichmentAnalyzer(pathways, background_genes=genes)
print(analyzer.analyze_components(components_p1, condition="P1").head())

# 7. Plots
vis = gtgmm.TerrainVisualizer()
vis.plot_terrain_3d(gt_p1, title="P1 GeneTerrain")
vis.plot_gmm_components(gt_p1, gmm_p1, title="P1 GMM Components")
```

For a full P1 vs P8 comparison, see [`examples/cardiac_regeneration_example.py`](examples/cardiac_regeneration_example.py) or call `gtgmm.analyze_regeneration(genes_p1, expr_p1, genes_p8, expr_p8)`.

Building a terrain calls the STRING API, so the first run needs network access.

---

## License

MIT License - see the `LICENSE` file for details.

