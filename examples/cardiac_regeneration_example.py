#!/usr/bin/env python3
"""
Cardiac Regeneration Analysis Example using gtGMM v2.0

This example demonstrates a complete analysis workflow for the GSE130699
cardiac regeneration dataset, comparing P1 (regenerative) and P8 (non-regenerative)
neonatal mouse hearts.

The analysis includes:
1. Loading and preprocessing data
2. Creating high-resolution GeneTerrains
3. Adaptive GMM decomposition
4. Topological data analysis (TDA)
5. Pathway enrichment analysis
6. Visualization of results
"""

import sys
import os
import numpy as np
import pandas as pd

# Add the package to the path if running from the examples directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import gtgmm

print("=" * 80)
print("Cardiac Regeneration Analysis using gtGMM v2.0")
print("=" * 80)

# ============================================================================
# STEP 1: Load Data
# ============================================================================
print("\n[STEP 1] Loading GSE130699 cardiac regeneration data...")

# For this example, we'll use a subset of genes for demonstration
# In a real analysis, you would use the full dataset via:
# genes, p1_expr, p8_expr = gtgmm.data.load_gse130699()

# Create a representative gene set for demonstration
pathways = gtgmm.data.get_cardiac_pathways()
demo_genes = []
for pathway_name, pathway_genes in pathways.items():
    demo_genes.extend(pathway_genes[:5])  # Take first 5 genes from each pathway

# Remove duplicates
demo_genes = list(set(demo_genes))
print(f"Using {len(demo_genes)} genes for demonstration")

# Create simulated expression data for demonstration
# In a real analysis, this would come from actual scRNA-seq data
np.random.seed(42)
p1_expr = {gene: np.random.uniform(1.0, 3.0) for gene in demo_genes}
p8_expr = {gene: np.random.uniform(0.5, 2.0) for gene in demo_genes}

print(f"✓ Loaded expression data for {len(demo_genes)} genes")

# ============================================================================
# STEP 2: Create GeneTerrains
# ============================================================================
print("\n[STEP 2] Creating high-resolution GeneTerrains...")

# Create P1 terrain
print("  Creating P1 (regenerative) terrain...")
gt_p1 = gtgmm.create_terrain(
    demo_genes, 
    p1_expr, 
    species=10090,  # Mouse
    resolution=500,  # High resolution
    score_threshold=400,
    optimize_sigma=True
)
print(f"  ✓ P1 terrain created: {gt_p1.terrain.shape}")

# Create P8 terrain
print("  Creating P8 (non-regenerative) terrain...")
gt_p8 = gtgmm.create_terrain(
    demo_genes, 
    p8_expr, 
    species=10090,
    resolution=500,
    score_threshold=400,
    optimize_sigma=True
)
print(f"  ✓ P8 terrain created: {gt_p8.terrain.shape}")

# ============================================================================
# STEP 3: GMM Decomposition
# ============================================================================
print("\n[STEP 3] Performing adaptive GMM decomposition...")

# P1 GMM analysis
print("  Analyzing P1 components...")
gmm_p1 = gtgmm.GMMOptimizer(max_components=8)
terrain_data_p1 = gt_p1.get_terrain_data()
n_components_p1 = gmm_p1.optimize(terrain_data_p1, min_components=3)
gmm_p1.fit(terrain_data_p1, n_components=n_components_p1)
print(f"  ✓ P1 optimal components: {n_components_p1}")

# P8 GMM analysis
print("  Analyzing P8 components...")
gmm_p8 = gtgmm.GMMOptimizer(max_components=8)
terrain_data_p8 = gt_p8.get_terrain_data()
n_components_p8 = gmm_p8.optimize(terrain_data_p8, min_components=3)
gmm_p8.fit(terrain_data_p8, n_components=n_components_p8)
print(f"  ✓ P8 optimal components: {n_components_p8}")

# Extract gene-to-component mappings
components_p1 = gmm_p1.extract_gene_components(gt_p1)
components_p8 = gmm_p8.extract_gene_components(gt_p8)

print(f"\n  Component sizes (P1):")
for comp_idx, comp_df in components_p1.items():
    print(f"    Component {comp_idx + 1}: {len(comp_df)} genes")

# ============================================================================
# STEP 4: Topological Data Analysis
# ============================================================================
print("\n[STEP 4] Computing topological data analysis (TDA) metrics...")

# P1 TDA
tda_p1 = gtgmm.TerrainTDA(gt_p1)
summary_p1 = tda_p1.get_tda_summary()
print("  ✓ P1 TDA complete")

# P8 TDA
tda_p8 = gtgmm.TerrainTDA(gt_p8)
summary_p8 = tda_p8.get_tda_summary()
print("  ✓ P8 TDA complete")

# Compare TDA features
comparison = gtgmm.compare_tda_features(summary_p1, summary_p8, 
                                       labels=['P1', 'P8'])
print("\n  TDA Comparison:")
print(comparison.to_string(index=False))

# ============================================================================
# STEP 5: Pathway Enrichment Analysis
# ============================================================================
print("\n[STEP 5] Performing pathway enrichment analysis...")

# Create enrichment analyzer
enrichment = gtgmm.EnrichmentAnalyzer(pathways, background_genes=demo_genes)

# Analyze both conditions
enrichment_results = enrichment.compare_conditions(
    components_p1, components_p8, 
    condition_a="P1", condition_b="P8"
)

print(f"  ✓ Enrichment analysis complete: {len(enrichment_results)} pathway-component pairs")

# Show top enrichments for each condition
print("\n  Top 5 enrichments in P1:")
p1_enrichments = enrichment_results[enrichment_results['condition'] == 'P1']
top_p1 = p1_enrichments.nsmallest(5, 'p_value')
for _, row in top_p1.iterrows():
    print(f"    {row['pathway']} (Component {row['component']}): "
          f"p={row['p_value']:.4e}, overlap={row['overlap']}/{row['pathway_size']}")

print("\n  Top 5 enrichments in P8:")
p8_enrichments = enrichment_results[enrichment_results['condition'] == 'P8']
top_p8 = p8_enrichments.nsmallest(5, 'p_value')
for _, row in top_p8.iterrows():
    print(f"    {row['pathway']} (Component {row['component']}): "
          f"p={row['p_value']:.4e}, overlap={row['overlap']}/{row['pathway_size']}")

# ============================================================================
# STEP 6: Visualization
# ============================================================================
print("\n[STEP 6] Generating visualizations...")

# Create visualizer
vis = gtgmm.TerrainVisualizer()

# Plot 3D terrains
print("  Creating 3D terrain plots...")
try:
    vis.plot_3d_terrain(gt_p1, title="P1 (Regenerative) GeneTerrain", 
                       save_path="p1_terrain_3d.png")
    vis.plot_3d_terrain(gt_p8, title="P8 (Non-Regenerative) GeneTerrain", 
                       save_path="p8_terrain_3d.png")
    print("  ✓ 3D terrain plots created")
except Exception as e:
    print(f"  Note: 3D plotting requires display: {e}")

# Plot component maps
print("  Creating component maps...")
try:
    vis.plot_gmm_components(gt_p1, gmm_p1, title="P1 GMM Components", 
                           save_path="p1_components.png")
    vis.plot_gmm_components(gt_p8, gmm_p8, title="P8 GMM Components", 
                           save_path="p8_components.png")
    print("  ✓ Component maps created")
except Exception as e:
    print(f"  Note: Component plotting requires display: {e}")

# ============================================================================
# Summary
# ============================================================================
print("\n" + "=" * 80)
print("Analysis Complete!")
print("=" * 80)
print(f"\nKey Findings:")
print(f"  - P1 has {n_components_p1} functional modules")
print(f"  - P8 has {n_components_p8} functional modules")
print(f"  - P1 mean Ricci curvature: {summary_p1['mean_ricci_curvature'].iloc[0]:.4f}")
print(f"  - P8 mean Ricci curvature: {summary_p8['mean_ricci_curvature'].iloc[0]:.4f}")
print(f"  - P1 Euler characteristic: {summary_p1['euler_characteristic'].iloc[0]}")
print(f"  - P8 Euler characteristic: {summary_p8['euler_characteristic'].iloc[0]}")

print("\nThis analysis demonstrates the power of gtGMM for unsupervised")
print("discovery of functional modules and topological differences in")
print("gene expression networks.")
print("=" * 80)

