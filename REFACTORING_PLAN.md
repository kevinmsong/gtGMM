# gtGMM Refactoring and Enhancement Plan

**Objective:** Refactor the `gtGMM` repository to create a clean, lightweight, and production-ready package. This new version will integrate all optimizations and enhancements developed during the cardiac regeneration analysis, including the seamless integration of `TerrainTDA`.

---

## 1. Guiding Principles

*   **Modularity:** Each core component (data, terrain, GMM, TDA, visualization) will be a self-contained module with a clear API.
*   **Usability:** The package should be easy to install and use, with a high-level API for common workflows.
*   **Performance:** Code will be optimized for speed and memory efficiency, especially for large datasets and high-resolution terrains.
*   **Extensibility:** The architecture will be designed to easily accommodate new features and algorithms.
*   **Production-Ready:** The package will include proper documentation, tests, and packaging infrastructure.

---

## 2. Proposed Package Structure

The new package structure will follow standard Python best practices:

```
/gtGMM/
├── .github/                # GitHub actions for CI/CD
│   └── workflows/
│       └── python-package.yml
├── docs/                   # User documentation (Sphinx)
│   ├── source/
│   │   ├── conf.py
│   │   ├── index.rst
│   │   └── ...
│   └── Makefile
├── examples/               # Jupyter notebooks and scripts
│   ├── 1_basic_usage.ipynb
│   └── 2_cardiac_regeneration_analysis.ipynb
├── gtgmm/                  # Main package source code
│   ├── __init__.py         # High-level API
│   ├── data.py             # Data loading and preprocessing
│   ├── terrain.py          # GeneTerrain creation and manipulation
│   ├── gmm.py              # GMM optimization and analysis
│   ├── tda.py              # TerrainTDA functionality
│   ├── enrichment.py       # Pathway enrichment analysis
│   ├── visualization.py    # Plotting and visualization
│   └── utils.py            # Utility functions
├── tests/                  # Unit and integration tests
│   ├── test_data.py
│   ├── test_terrain.py
│   └── ...
├── .gitignore
├── LICENSE
├── README.md
├── requirements.txt
└── setup.py
```

---

## 3. Core Module Refactoring

### `gtgmm.terrain`

*   **`GeneTerrain` Class:** This will remain the central class for creating and managing terrains.
*   **Optimizations:**
    *   Default resolution will be set to `500x500`.
    *   `optimize_sigma` will be integrated as a standard method.
    *   The `create_terrain` method will be streamlined.

### `gtgmm.gmm`

*   **`GMMAnalyzer` Class:** A new class to handle GMM-related tasks.
*   **Optimizations:**
    *   The `GMMOptimizer` logic will be integrated into `GMMAnalyzer`.
    *   The `fit` method will automatically perform adaptive GMM optimization (BIC-based model selection).

### `gtgmm.tda`

*   **`TerrainTDA` Class:** This class will be enhanced and tightly integrated.
*   **Integration:**
    *   The `GeneTerrain` object will have a `.tda` property that returns a `TerrainTDA` instance, allowing for a more intuitive workflow (e.g., `gt.tda.get_summary()`).

### `gtgmm.enrichment`

*   **`EnrichmentAnalyzer` Class:** A new class for pathway enrichment analysis.
*   **Features:**
    *   Will support both standard and custom gene sets.
    *   Will include methods for Fisher's exact test and FDR correction.

### `gtgmm.visualization`

*   **`Visualizer` Class:** A dedicated class for generating all plots.
*   **Features:**
    *   Methods for plotting 2D and 3D terrains, component maps, and enrichment results.
    *   Publication-quality aesthetics by default (e.g., `jet` colormap, clear labels).

---

## 4. API Design

The goal is a clean and intuitive API. Here is a proposed high-level workflow:

```python
import gtgmm

# 1. Load data
genes, p1_expr, p8_expr = gtgmm.data.load_gse130699()

# 2. Create and analyze P1 terrain
gt_p1 = gtgmm.GeneTerrain(genes, p1_expr)
gt_p1.build()

# 3. Perform GMM and TDA
p1_gmm = gt_p1.run_gmm()
p1_tda = gt_p1.run_tda()

# 4. Run enrichment analysis
enrichment_analyzer = gtgmm.EnrichmentAnalyzer(pathways)
p1_enrichment = enrichment_analyzer.run(p1_gmm.get_components())

# 5. Visualize results
vis = gtgmm.Visualizer()
vis.plot_3d_terrain(gt_p1)
vis.plot_component_map(p1_gmm)
```

---

## 5. Implementation Plan

1.  **Create new repository structure:** Set up the directories and files as outlined above.
2.  **Migrate and refactor core modules:** Move the existing code into the new structure and refactor it according to the plan.
3.  **Integrate optimizations:** Incorporate the high-resolution terrains, adaptive GMM, and other enhancements.
4.  **Develop `EnrichmentAnalyzer`:** Create the new module for pathway enrichment.
5.  **Refine `TerrainTDA` integration:** Implement the seamless `.tda` property on the `GeneTerrain` object.
6.  **Write unit tests:** Develop a comprehensive test suite to ensure correctness and stability.
7.  **Create documentation:** Write user guides and API documentation.
8.  **Develop example notebooks:** Create new examples that demonstrate the refactored package's features.

This plan will guide the refactoring process, resulting in a robust and user-friendly `gtGMM` package that is ready for broader distribution and use in the scientific community.

