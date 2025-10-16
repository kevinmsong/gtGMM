"""
gtGMM: GeneTerrain Gaussian Mixture Models

Topological Data Analysis for Gene Expression Networks

A Python library for unsupervised topological analysis of gene networks using
GeneTerrain spatial embeddings, Gaussian Mixture Model decomposition, and 
advanced topological measures.

Authors: Kevin Song, John Zhang, Lei Ye, Jianyi Zhang
Affiliation: Department of Biomedical Engineering, The University of Alabama at Birmingham

Key Features:
- High-resolution GeneTerrain generation (500x500 default)
- Adaptive GMM optimization with BIC model selection
- Comprehensive topological data analysis (TDA)
- Pathway enrichment analysis with FDR correction
- Publication-quality visualizations

Example Usage:
    >>> import gtgmm
    >>> 
    >>> # Load cardiac regeneration data
    >>> genes, p1_expr, p8_expr = gtgmm.data.load_gse130699()
    >>> 
    >>> # Create and analyze P1 terrain
    >>> gt_p1 = gtgmm.GeneTerrain(genes, p1_expr, resolution=500)
    >>> gt_p1.build()
    >>> 
    >>> # Run GMM and TDA
    >>> gmm_results = gt_p1.run_gmm()
    >>> tda_results = gt_p1.run_tda()
    >>> 
    >>> # Pathway enrichment
    >>> pathways = gtgmm.data.get_cardiac_pathways()
    >>> enrichment = gtgmm.EnrichmentAnalyzer(pathways)
    >>> results = enrichment.analyze_components(gmm_results.get_components())
"""

__version__ = "2.0.0"
__author__ = "Kevin Song, John Zhang, Lei Ye, Jianyi Zhang"
__email__ = "kmsong@uab.edu"
__institution__ = "The University of Alabama at Birmingham"

# Core classes
from .terrain import GeneTerrain
from .gmm import GMMOptimizer
from .tda import TerrainTDA, compare_tda_features
from .enrichment import EnrichmentAnalyzer
from .visualization import TerrainVisualizer

# Submodules
from . import data
from . import utils

__all__ = [
    # Core classes
    "GeneTerrain",
    "GMMOptimizer",
    "TerrainTDA",
    "EnrichmentAnalyzer",
    "TerrainVisualizer",
    
    # Functions
    "compare_tda_features",
    
    # Submodules
    "data",
    "utils",
]


# Convenience functions for common workflows
def create_terrain(genes, expression, species=10090, resolution=500, 
                  score_threshold=400, optimize_sigma=True):
    """
    Convenience function to create a GeneTerrain with optimal settings.
    
    Parameters
    ----------
    genes : list
        List of gene symbols
    expression : dict
        Gene expression values
    species : int
        NCBI taxonomy ID (10090 for mouse, 9606 for human)
    resolution : int
        Terrain resolution (default: 500)
    score_threshold : int
        STRING interaction score threshold (default: 400)
    optimize_sigma : bool
        Whether to optimize sigma parameter (default: True)
        
    Returns
    -------
    GeneTerrain
        Fully built terrain object
    """
    gt = GeneTerrain(genes, expression, species=species, resolution=resolution)
    gt.fetch_string_interactions(score_threshold=score_threshold)
    gt.compute_layout()
    
    if optimize_sigma:
        gt.optimize_sigma(n_samples=15)
    else:
        gt.create_terrain()
    
    return gt


def analyze_regeneration(genes_p1, expr_p1, genes_p8, expr_p8, 
                        pathways=None, species=10090):
    """
    Complete cardiac regeneration analysis workflow.
    
    Parameters
    ----------
    genes_p1 : list
        P1 gene symbols
    expr_p1 : dict
        P1 expression values
    genes_p8 : list
        P8 gene symbols
    expr_p8 : dict
        P8 expression values
    pathways : dict, optional
        Pathway gene sets. If None, uses cardiac pathways.
    species : int
        NCBI taxonomy ID
        
    Returns
    -------
    dict
        Complete analysis results including terrains, GMM, TDA, and enrichment
    """
    # Create terrains
    gt_p1 = create_terrain(genes_p1, expr_p1, species=species)
    gt_p8 = create_terrain(genes_p8, expr_p8, species=species)
    
    # GMM analysis
    gmm_p1 = GMMOptimizer(max_components=8)
    gmm_p1.optimize(gt_p1.get_terrain_data(), min_components=3)
    
    gmm_p8 = GMMOptimizer(max_components=8)
    gmm_p8.optimize(gt_p8.get_terrain_data(), min_components=3)
    
    # TDA analysis
    tda_p1 = TerrainTDA(gt_p1)
    tda_p8 = TerrainTDA(gt_p8)
    
    # Pathway enrichment
    if pathways is None:
        pathways = data.get_cardiac_pathways()
    
    enrichment = EnrichmentAnalyzer(pathways, background_genes=genes_p1)
    
    # Extract components
    comps_p1 = gmm_p1.extract_gene_components(gt_p1)
    comps_p8 = gmm_p8.extract_gene_components(gt_p8)
    
    enrich_results = enrichment.compare_conditions(comps_p1, comps_p8, 'P1', 'P8')
    
    return {
        'terrains': {'P1': gt_p1, 'P8': gt_p8},
        'gmm': {'P1': gmm_p1, 'P8': gmm_p8},
        'tda': {'P1': tda_p1, 'P8': tda_p8},
        'enrichment': enrich_results,
        'components': {'P1': comps_p1, 'P8': comps_p8}
    }

