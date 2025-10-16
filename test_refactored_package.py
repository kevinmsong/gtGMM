#!/usr/bin/env python3
"""
Test script for the refactored gtGMM package.

This script tests the core functionality of the refactored package,
including terrain creation, GMM analysis, TDA, and enrichment analysis.
"""

import sys
import os

# Add the package to the path
sys.path.insert(0, os.path.abspath('.'))

import numpy as np
import pandas as pd

print("=" * 80)
print("Testing Refactored gtGMM Package v2.0")
print("=" * 80)

# Test 1: Import all modules
print("\n[TEST 1] Importing modules...")
try:
    import gtgmm
    from gtgmm import GeneTerrain, GMMOptimizer, TerrainTDA, EnrichmentAnalyzer, TerrainVisualizer
    from gtgmm import data, utils
    print("✓ All modules imported successfully")
except Exception as e:
    print(f"✗ Import failed: {e}")
    sys.exit(1)

# Test 2: Load cardiac pathways
print("\n[TEST 2] Loading cardiac pathways...")
try:
    pathways = data.get_cardiac_pathways()
    print(f"✓ Loaded {len(pathways)} pathways")
    for name, genes in list(pathways.items())[:3]:
        print(f"  - {name}: {len(genes)} genes")
except Exception as e:
    print(f"✗ Pathway loading failed: {e}")
    sys.exit(1)

# Test 3: Create a simple GeneTerrain
print("\n[TEST 3] Creating a simple GeneTerrain...")
try:
    # Use a small subset of genes for testing
    test_genes = ['Yap1', 'Wwtr1', 'Tead1', 'Tead2', 'Ccnd1', 'Ccnd2', 
                  'Akt1', 'Akt2', 'Wnt1', 'Wnt3a', 'Notch1', 'Notch2']
    
    # Create random expression values
    np.random.seed(42)
    test_expr = {gene: np.random.uniform(0.5, 2.0) for gene in test_genes}
    
    # Create terrain with lower resolution for testing
    gt = GeneTerrain(test_genes, test_expr, species=10090, resolution=100)
    print(f"✓ GeneTerrain created with {len(test_genes)} genes")
    print(f"  Resolution: {gt.resolution}x{gt.resolution}")
except Exception as e:
    print(f"✗ GeneTerrain creation failed: {e}")
    sys.exit(1)

# Test 4: Build the terrain (fetch interactions, layout, create terrain)
print("\n[TEST 4] Building terrain (fetch, layout, create)...")
try:
    # Fetch interactions
    gt.fetch_string_interactions(score_threshold=400)
    print(f"✓ Fetched {len(gt.interactions) if gt.interactions is not None else 0} interactions")
    
    # Compute layout
    gt.compute_layout()
    print(f"✓ Computed layout for {len(gt.gene_data)} genes")
    
    # Create terrain
    gt.create_terrain(base_sigma=0.1, smooth_sigma=1.0)
    print(f"✓ Created terrain: {gt.terrain.shape}")
except Exception as e:
    print(f"✗ Terrain building failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 5: GMM Analysis
print("\n[TEST 5] Running GMM analysis...")
try:
    terrain_data = gt.get_terrain_data()
    gmm = GMMOptimizer(max_components=5)
    n_components = gmm.optimize(terrain_data, min_components=2)
    print(f"✓ Optimal number of components: {n_components}")
    
    gmm.fit(terrain_data, n_components=n_components)
    print(f"✓ GMM fitted with {n_components} components")
except Exception as e:
    print(f"✗ GMM analysis failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 6: TDA Analysis
print("\n[TEST 6] Running TDA analysis...")
try:
    tda = TerrainTDA(gt)
    summary = tda.get_tda_summary()
    print(f"✓ TDA summary computed")
    print(f"  Mean Ricci curvature: {summary['mean_ricci_curvature'].iloc[0]:.4f}")
    print(f"  Betti numbers: β₀={summary['median_beta_0'].iloc[0]}, β₁={summary['median_beta_1'].iloc[0]}")
    print(f"  Euler characteristic: {summary['euler_characteristic'].iloc[0]}")
except Exception as e:
    print(f"✗ TDA analysis failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 7: Enrichment Analysis
print("\n[TEST 7] Running enrichment analysis...")
try:
    enrichment = EnrichmentAnalyzer(pathways, background_genes=test_genes)
    
    # Extract components
    components = gmm.extract_gene_components(gt)
    print(f"✓ Extracted {len(components)} components")
    
    # Run enrichment
    results = enrichment.analyze_components(components, condition="Test")
    print(f"✓ Enrichment analysis complete: {len(results)} pathway-component pairs")
    
    # Show top enrichments
    if len(results) > 0:
        top_results = results.nsmallest(3, 'p_value')
        print("\n  Top 3 enrichments:")
        for _, row in top_results.iterrows():
            print(f"    - {row['pathway']} (Component {row['component']}): p={row['p_value']:.4e}")
except Exception as e:
    print(f"✗ Enrichment analysis failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Test 8: Convenience function
print("\n[TEST 8] Testing convenience function...")
try:
    gt_conv = gtgmm.create_terrain(test_genes, test_expr, species=10090, 
                                   resolution=100, optimize_sigma=False)
    print(f"✓ Convenience function works: terrain shape {gt_conv.terrain.shape}")
except Exception as e:
    print(f"✗ Convenience function failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("\n" + "=" * 80)
print("All tests passed! ✓")
print("=" * 80)

