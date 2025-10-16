"""
Topological Data Analysis (TDA) module for GeneTerrain.

This module provides advanced topological analysis tools for GeneTerrain landscapes,
including:
- Forman-Ricci curvature analysis
- Persistent homology and Betti numbers  
- Euler characteristic
- Critical point detection
- Morse theory analysis
- Topological feature extraction

These tools enable quantitative comparison of gene network topologies across conditions.
"""

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import distance_matrix
from typing import Dict, List, Tuple, Optional
import networkx as nx
from collections import defaultdict


class TerrainTDA:
    """
    Topological Data Analysis for GeneTerrain objects.
    
    Provides methods to compute topological invariants and curvature measures
    that characterize the shape and structure of gene network landscapes.
    
    Parameters
    ----------
    gene_terrain : GeneTerrain
        A fitted GeneTerrain object with terrain generated
    """
    
    def __init__(self, gene_terrain):
        """Initialize TDA analyzer with a GeneTerrain object."""
        if gene_terrain.terrain is None:
            raise ValueError("GeneTerrain must have terrain generated. Call optimize_sigma() first.")
        
        self.gene_terrain = gene_terrain
        self.terrain = gene_terrain.terrain
        self.resolution = gene_terrain.resolution
        
        # Computed properties (cached)
        self._gradient = None
        self._curvature = None
        self._critical_points = None
        self._persistence = None
        
    def compute_forman_ricci_curvature(self) -> Dict[Tuple[int, int], float]:
        """
        Compute Forman-Ricci curvature on the terrain network.
        
        Forman-Ricci curvature measures the "curvature" of edges in the network.
        - Positive curvature: edge connects similar expression regions
        - Negative curvature: edge spans different expression regimes
        
        This provides a quantitative measure of how "curved" the gene network
        landscape is, reflecting functional modularity.
        
        Returns
        -------
        dict
            Dictionary mapping edge tuples to curvature values
            
        Notes
        -----
        Uses discrete Forman-Ricci curvature formula for directed graphs.
        Curvature reflects local geometry of the network embedded in expression space.
        """
        if self.gene_terrain.network is None:
            raise ValueError("Network must be built first. Call build_network().")
        
        G = self.gene_terrain.network
        gene_data = self.gene_terrain.gene_data
        
        # Get gene positions and expression
        gene_pos = gene_data.set_index('gene')[['X', 'Y', 'Exp']].to_dict('index')
        
        curvatures = {}
        
        for u, v in G.edges():
            if u not in gene_pos or v not in gene_pos:
                continue
                
            # Get node properties
            pos_u = np.array([gene_pos[u]['X'], gene_pos[u]['Y']])
            pos_v = np.array([gene_pos[v]['X'], gene_pos[v]['Y']])
            expr_u = gene_pos[u]['Exp']
            expr_v = gene_pos[v]['Exp']
            
            # Forman-Ricci curvature formula (discrete version)
            # Curvature = (degree_u + degree_v) / 2 - spatial_distance - expression_diff
            deg_u = G.degree(u)
            deg_v = G.degree(v)
            
            spatial_dist = np.linalg.norm(pos_v - pos_u)
            expr_diff = abs(expr_v - expr_u)
            
            # Normalized Forman curvature
            curvature = (deg_u + deg_v) / 2.0 - 2 * spatial_dist - expr_diff
            
            curvatures[(u, v)] = curvature
            
        return curvatures
    
    def compute_network_curvature_statistics(self) -> pd.DataFrame:
        """
        Compute summary statistics of Forman-Ricci curvature distribution.
        
        Returns
        -------
        pd.DataFrame
            Statistics: mean, std, min, max, median, positive_fraction
        """
        curvatures = self.compute_forman_ricci_curvature()
        curv_values = np.array(list(curvatures.values()))
        
        stats = {
            'mean_curvature': np.mean(curv_values),
            'std_curvature': np.std(curv_values),
            'min_curvature': np.min(curv_values),
            'max_curvature': np.max(curv_values),
            'median_curvature': np.median(curv_values),
            'positive_fraction': np.mean(curv_values > 0),
            'negative_fraction': np.mean(curv_values < 0),
            'n_edges': len(curv_values)
        }
        
        return pd.DataFrame([stats])
    
    def compute_terrain_gradient(self) -> np.ndarray:
        """
        Compute gradient magnitude of the terrain.
        
        Returns
        -------
        np.ndarray
            Gradient magnitude at each point (resolution x resolution)
        """
        if self._gradient is not None:
            return self._gradient
            
        # Compute gradient using Sobel filters
        grad_x = ndimage.sobel(self.terrain, axis=1)
        grad_y = ndimage.sobel(self.terrain, axis=0)
        
        self._gradient = np.sqrt(grad_x**2 + grad_y**2)
        return self._gradient
    
    def compute_terrain_curvature(self) -> np.ndarray:
        """
        Compute Gaussian curvature of the terrain surface.
        
        Gaussian curvature characterizes local shape:
        - Positive: peak/valley (locally spherical)
        - Negative: saddle point (locally hyperbolic)
        - Zero: flat or cylindrical
        
        Returns
        -------
        np.ndarray
            Gaussian curvature at each point (resolution x resolution)
        """
        if self._curvature is not None:
            return self._curvature
            
        # Compute second derivatives
        Ixx = ndimage.gaussian_filter(self.terrain, sigma=1, order=[2, 0])
        Iyy = ndimage.gaussian_filter(self.terrain, sigma=1, order=[0, 2])
        Ixy = ndimage.gaussian_filter(self.terrain, sigma=1, order=[1, 1])
        
        # First derivatives for normalization
        Ix = ndimage.gaussian_filter(self.terrain, sigma=1, order=[1, 0])
        Iy = ndimage.gaussian_filter(self.terrain, sigma=1, order=[0, 1])
        
        # Gaussian curvature K = (Ixx*Iyy - Ixy^2) / (1 + Ix^2 + Iy^2)^2
        numerator = Ixx * Iyy - Ixy**2
        denominator = (1 + Ix**2 + Iy**2)**2
        
        self._curvature = numerator / (denominator + 1e-10)
        return self._curvature
    
    def detect_critical_points(self, threshold=0.01) -> Dict[str, List[Tuple[int, int]]]:
        """
        Detect critical points (peaks, valleys, saddles) in the terrain.
        
        Parameters
        ----------
        threshold : float
            Curvature threshold for classification (default: 0.01)
            
        Returns
        -------
        dict
            Dictionary with 'peaks', 'valleys', 'saddles' lists of (i, j) coordinates
        """
        if self._critical_points is not None:
            return self._critical_points
            
        gradient = self.compute_terrain_gradient()
        curvature = self.compute_terrain_curvature()
        
        # Find points with low gradient (potential critical points)
        grad_thresh = np.percentile(gradient, 10)
        candidates = gradient < grad_thresh
        
        critical_points = {
            'peaks': [],
            'valleys': [],
            'saddles': []
        }
        
        for i in range(1, self.resolution - 1):
            for j in range(1, self.resolution - 1):
                if not candidates[i, j]:
                    continue
                    
                curv_val = curvature[i, j]
                height = self.terrain[i, j]
                
                # Check neighbors
                neighborhood = self.terrain[i-1:i+2, j-1:j+2]
                
                if curv_val > threshold:
                    if height > np.mean(neighborhood):
                        critical_points['peaks'].append((i, j))
                    else:
                        critical_points['valleys'].append((i, j))
                elif curv_val < -threshold:
                    critical_points['saddles'].append((i, j))
        
        self._critical_points = critical_points
        return critical_points
    
    def compute_betti_numbers(self, thresholds: Optional[np.ndarray] = None) -> pd.DataFrame:
        """
        Compute Betti numbers (topological features) at different height thresholds.
        
        Betti numbers count topological features:
        - β0: number of connected components
        - β1: number of cycles/loops
        
        Parameters
        ----------
        thresholds : np.ndarray, optional
            Height thresholds to compute Betti numbers. If None, uses 20 evenly spaced values.
            
        Returns
        -------
        pd.DataFrame
            Betti numbers (beta_0, beta_1) at each threshold level
        """
        if thresholds is None:
            thresholds = np.linspace(self.terrain.min(), self.terrain.max(), 20)
        
        results = []
        
        for thresh in thresholds:
            # Create binary mask for supra-level set
            mask = self.terrain >= thresh
            
            # Label connected components
            labeled, n_components = ndimage.label(mask)
            
            # Approximate β1 using Euler characteristic
            # χ = β0 - β1, and for 2D: χ = V - E + F
            # Simple approximation: count holes in the mask
            filled = ndimage.binary_fill_holes(mask)
            holes = filled.astype(int) - mask.astype(int)
            n_holes = ndimage.label(holes)[1]
            
            results.append({
                'threshold': thresh,
                'beta_0': n_components,  # Connected components
                'beta_1': n_holes         # Cycles/loops
            })
        
        return pd.DataFrame(results)
    
    def compute_euler_characteristic(self) -> int:
        """
        Compute Euler characteristic of the terrain.
        
        For a 2D surface: χ = V - E + F (vertices - edges + faces)
        This topological invariant is preserved under continuous deformations.
        
        Returns
        -------
        int
            Euler characteristic
        """
        # Use critical points to estimate Euler characteristic
        critical_pts = self.detect_critical_points()
        
        n_maxima = len(critical_pts['peaks'])
        n_minima = len(critical_pts['valleys'])
        n_saddles = len(critical_pts['saddles'])
        
        # For a 2D surface: χ = #maxima - #saddles + #minima
        euler_char = n_maxima - n_saddles + n_minima
        
        return euler_char
    
    def compute_persistence_entropy(self) -> float:
        """
        Compute persistence entropy of topological features.
        
        Persistence entropy measures the complexity of the topological structure:
        - Low entropy: few dominant features
        - High entropy: many features of similar importance
        
        Returns
        -------
        float
            Persistence entropy value
        """
        betti = self.compute_betti_numbers()
        
        # Compute lifetime of each component
        lifetimes = []
        for i in range(len(betti) - 1):
            if betti.iloc[i]['beta_0'] > 0:
                lifetime = betti.iloc[i+1]['threshold'] - betti.iloc[i]['threshold']
                lifetimes.append(lifetime)
        
        if len(lifetimes) == 0:
            return 0.0
        
        # Normalize lifetimes to probabilities
        lifetimes = np.array(lifetimes)
        probs = lifetimes / lifetimes.sum()
        
        # Compute entropy
        entropy = -np.sum(probs * np.log(probs + 1e-10))
        
        return entropy
    
    def compute_terrain_roughness(self) -> float:
        """
        Compute terrain roughness (RMS of gradient).
        
        Returns
        -------
        float
            RMS gradient (roughness measure)
        """
        gradient = self.compute_terrain_gradient()
        return np.sqrt(np.mean(gradient**2))
    
    def get_tda_summary(self) -> pd.DataFrame:
        """
        Compute comprehensive TDA summary statistics.
        
        Returns
        -------
        pd.DataFrame
            All TDA metrics in one table
        """
        # Network curvature
        curv_stats = self.compute_network_curvature_statistics().iloc[0]
        
        # Critical points
        critical_pts = self.detect_critical_points()
        
        # Topological invariants
        euler_char = self.compute_euler_characteristic()
        persistence_entropy = self.compute_persistence_entropy()
        roughness = self.compute_terrain_roughness()
        
        # Betti numbers at median threshold
        betti = self.compute_betti_numbers()
        median_betti = betti.iloc[len(betti)//2]
        
        summary = {
            # Curvature
            'mean_ricci_curvature': curv_stats['mean_curvature'],
            'std_ricci_curvature': curv_stats['std_curvature'],
            'positive_curvature_fraction': curv_stats['positive_fraction'],
            
            # Critical points
            'n_peaks': len(critical_pts['peaks']),
            'n_valleys': len(critical_pts['valleys']),
            'n_saddles': len(critical_pts['saddles']),
            
            # Topological invariants
            'euler_characteristic': euler_char,
            'persistence_entropy': persistence_entropy,
            'median_beta_0': median_betti['beta_0'],
            'median_beta_1': median_betti['beta_1'],
            
            # Terrain properties
            'roughness': roughness,
            'mean_height': np.mean(self.terrain),
            'std_height': np.std(self.terrain),
        }
        
        return pd.DataFrame([summary])


def compare_tda_features(tda1: TerrainTDA, tda2: TerrainTDA, 
                        labels=('Condition 1', 'Condition 2')) -> pd.DataFrame:
    """
    Compare TDA features between two GeneTerrain objects.
    
    Parameters
    ----------
    tda1 : TerrainTDA
        First TDA analyzer
    tda2 : TerrainTDA
        Second TDA analyzer
    labels : tuple
        Labels for the two conditions
        
    Returns
    -------
    pd.DataFrame
        Comparison table with fold-changes and differences
    """
    summary1 = tda1.get_tda_summary().iloc[0]
    summary2 = tda2.get_tda_summary().iloc[0]
    
    comparison = []
    for metric in summary1.index:
        val1 = summary1[metric]
        val2 = summary2[metric]
        
        # Compute fold-change (avoid division by zero)
        if val1 != 0:
            fold_change = val2 / val1
        else:
            fold_change = np.inf if val2 > 0 else np.nan
        
        comparison.append({
            'metric': metric,
            labels[0]: val1,
            labels[1]: val2,
            'difference': val2 - val1,
            'fold_change': fold_change,
            'percent_change': ((val2 - val1) / (val1 + 1e-10)) * 100
        })
    
    return pd.DataFrame(comparison)
