"""
Utility functions for gtgmm package.
"""

import pandas as pd
import numpy as np


def export_gene_lists(component_genes_dict, output_dir='.', prefix='component'):
    """
    Export gene lists for each component to CSV files.
    
    Parameters
    ----------
    component_genes_dict : dict
        Dictionary mapping component index to gene DataFrame
    output_dir : str
        Output directory path
    prefix : str
        Prefix for output files
        
    Returns
    -------
    list
        List of created file paths
    """
    import os
    
    os.makedirs(output_dir, exist_ok=True)
    
    files = []
    for component_idx, genes_df in component_genes_dict.items():
        filename = os.path.join(output_dir, f'{prefix}_{component_idx}_genes.csv')
        genes_df.to_csv(filename, index=False)
        files.append(filename)
        print(f"Saved {len(genes_df)} genes to {filename}")
    
    return files


def create_summary_report(gene_terrain, gmm_optimizer, output_file='gtgmm_summary.txt'):
    """
    Create a text summary report of the analysis.
    
    Parameters
    ----------
    gene_terrain : GeneTerrain
        Fitted GeneTerrain object
    gmm_optimizer : GMMOptimizer
        Fitted GMMOptimizer object
    output_file : str
        Output file path
        
    Returns
    -------
    str
        Path to created file
    """
    with open(output_file, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("gtGMM Analysis Summary Report\n")
        f.write("=" * 70 + "\n\n")
        
        # GeneTerrain info
        f.write("GeneTerrain Information:\n")
        f.write("-" * 70 + "\n")
        f.write(f"Number of seed genes: {len(gene_terrain.seed_genes)}\n")
        if gene_terrain.network:
            f.write(f"Network nodes: {gene_terrain.network.number_of_nodes()}\n")
            f.write(f"Network edges: {gene_terrain.network.number_of_edges()}\n")
        f.write(f"Terrain resolution: {gene_terrain.resolution} x {gene_terrain.resolution}\n")
        if gene_terrain.optimal_sigma:
            f.write(f"Optimal sigma: {gene_terrain.optimal_sigma:.4f}\n")
        f.write("\n")
        
        # GMM info
        if gmm_optimizer.optimal_n_components:
            f.write("GMM Clustering Information:\n")
            f.write("-" * 70 + "\n")
            f.write(f"Optimal number of components: {gmm_optimizer.optimal_n_components}\n")
            
            # Component summary
            summary = gmm_optimizer.get_component_summary()
            f.write("\nComponent Summary:\n")
            for _, row in summary.iterrows():
                f.write(f"  Component {row['component']}: "
                       f"weight={row['weight']:.3f}, "
                       f"points={row['n_points']}, "
                       f"percentage={row['percentage']:.1f}%\n")
            f.write("\n")
        
        # Top genes by component
        f.write("Top Genes by Component:\n")
        f.write("-" * 70 + "\n")
        component_genes = gmm_optimizer.get_all_component_genes(gene_terrain, threshold=0.5)
        for comp_idx, genes_df in component_genes.items():
            f.write(f"\nComponent {comp_idx} ({len(genes_df)} genes):\n")
            if len(genes_df) > 0:
                top_genes = genes_df.head(10)['gene'].tolist()
                f.write(f"  Top genes: {', '.join(top_genes)}\n")
        
        f.write("\n" + "=" * 70 + "\n")
    
    print(f"Summary report saved to {output_file}")
    return output_file


def run_enrichment_analysis(gene_list, gene_sets='KEGG_2019_Human', organism='Human'):
    """
    Run enrichment analysis using gseapy (if installed).
    
    Parameters
    ----------
    gene_list : list
        List of gene symbols
    gene_sets : str
        Gene set library name
    organism : str
        Organism name
        
    Returns
    -------
    pd.DataFrame or None
        Enrichment results, or None if gseapy not available
    """
    try:
        import gseapy as gp
        
        print(f"Running enrichment analysis for {len(gene_list)} genes...")
        
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=None
        )
        
        if not enr.results.empty:
            print(f"Found {len(enr.results)} enriched terms")
            return enr.results
        else:
            print("No significant enrichment found")
            return None
            
    except ImportError:
        print("gseapy not installed. Install with: pip install gseapy")
        return None
    except Exception as e:
        print(f"Error running enrichment: {e}")
        return None


def get_component_statistics(gene_terrain, gmm_optimizer):
    """
    Calculate comprehensive statistics for each GMM component.
    
    Parameters
    ----------
    gene_terrain : GeneTerrain
        GeneTerrain object
    gmm_optimizer : GMMOptimizer
        Fitted GMMOptimizer
        
    Returns
    -------
    pd.DataFrame
        Statistics for each component
    """
    component_genes = gmm_optimizer.get_all_component_genes(gene_terrain, threshold=0.5)
    
    stats = []
    for comp_idx, genes_df in component_genes.items():
        if len(genes_df) > 0:
            stat = {
                'component': comp_idx,
                'n_genes': len(genes_df),
                'mean_expression': genes_df['expression'].mean(),
                'std_expression': genes_df['expression'].std(),
                'min_expression': genes_df['expression'].min(),
                'max_expression': genes_df['expression'].max(),
                'mean_probability': genes_df['probability'].mean(),
                'genes': ', '.join(genes_df.head(10)['gene'].tolist())
            }
        else:
            stat = {
                'component': comp_idx,
                'n_genes': 0,
                'mean_expression': 0,
                'std_expression': 0,
                'min_expression': 0,
                'max_expression': 0,
                'mean_probability': 0,
                'genes': ''
            }
        stats.append(stat)
    
    return pd.DataFrame(stats)


def quick_analysis(seed_genes, expression_values=None, optimize_sigma=True, 
                  optimize_gmm=True, plot=True):
    """
    Quick end-to-end analysis pipeline.
    
    Parameters
    ----------
    seed_genes : list
        List of gene symbols
    expression_values : dict or pd.Series, optional
        Gene expression values
    optimize_sigma : bool
        Whether to optimize sigma parameter
    optimize_gmm : bool
        Whether to optimize GMM components
    plot : bool
        Whether to generate plots
        
    Returns
    -------
    tuple
        (GeneTerrain, GMMOptimizer, TerrainVisualizer)
    """
    from .terrain import GeneTerrain
    from .gmm_optimizer import GMMOptimizer
    from .visualization import TerrainVisualizer
    
    print("Starting quick gtGMM analysis...")
    print("=" * 70)
    
    # Step 1: Create GeneTerrain
    print("\n1. Creating GeneTerrain...")
    gt = GeneTerrain(seed_genes, expression_values=expression_values)
    gt.fetch_string_interactions()
    gt.build_network()
    gt.compute_layout()
    
    # Step 2: Optimize and create terrain
    if optimize_sigma:
        print("\n2. Optimizing sigma...")
        gt.optimize_sigma(n_samples=20)
    else:
        print("\n2. Creating terrain with default sigma...")
        gt.create_terrain(base_sigma=0.15)
    
    # Step 3: Fit GMM
    print("\n3. Fitting GMM...")
    gmm = GMMOptimizer(max_components=8)
    terrain_data = gt.get_terrain_data()
    
    if optimize_gmm:
        gmm.optimize(terrain_data, min_components=2)
    else:
        gmm.fit(terrain_data, n_components=3)
    
    # Step 4: Visualize
    if plot:
        print("\n4. Creating visualizations...")
        viz = TerrainVisualizer()
        
        import matplotlib.pyplot as plt
        
        # Plot terrain
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        viz.plot_terrain_2d(gt, show_genes=True, show_labels=False, ax=axes[0])
        viz.plot_gmm_components(gmm_optimizer=gmm, gene_terrain=gt, ax=axes[1])
        plt.tight_layout()
        plt.show()
        
        # Plot optimization metrics
        if optimize_gmm:
            viz.plot_optimization_metrics(gmm, metrics='clustering')
            plt.show()
    else:
        viz = TerrainVisualizer()
    
    print("\n" + "=" * 70)
    print("Analysis complete!")
    
    return gt, gmm, viz
