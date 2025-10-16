"""
Visualization utilities for GeneTerrain and GMM results.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
from mpl_toolkits.mplot3d import Axes3D
import networkx as nx
from scipy.stats import multivariate_normal


class TerrainVisualizer:
    """
    Visualization toolkit for GeneTerrain analysis.
    
    Parameters
    ----------
    figsize : tuple
        Default figure size (default: (12, 10))
    dpi : int
        Default DPI for figures (default: 100)
    style : str
        Matplotlib style to use (default: 'default' for light theme)
    """
    
    def __init__(self, figsize=(12, 10), dpi=100, style='default'):
        self.figsize = figsize
        self.dpi = dpi
        
        # Set light theme with clean white background
        plt.style.use('default')
        plt.rcParams['figure.facecolor'] = 'white'
        plt.rcParams['axes.facecolor'] = 'white'
        plt.rcParams['savefig.facecolor'] = 'white'
        plt.rcParams['axes.edgecolor'] = 'black'
        plt.rcParams['axes.grid'] = False
        plt.rcParams['grid.alpha'] = 0.3
    
    def _compute_gaussian_contours(self, mean, covariance, resolution=100):
        """
        Compute probability density surface for a single Gaussian component.
        
        Parameters
        ----------
        mean : array-like, shape (2,)
            Mean of the Gaussian (X, Y coordinates)
        covariance : array-like, shape (2, 2)
            Covariance matrix
        resolution : int
            Grid resolution for evaluation
            
        Returns
        -------
        X, Y, Z : ndarrays
            Meshgrid coordinates and probability densities
        """
        # Create grid spanning [0, 1] x [0, 1]
        x = np.linspace(0, 1, resolution)
        y = np.linspace(0, 1, resolution)
        X, Y = np.meshgrid(x, y)
        
        # Stack coordinates for multivariate_normal
        pos = np.dstack((X, Y))
        
        # Compute probability density
        rv = multivariate_normal(mean, covariance, allow_singular=True)
        Z = rv.pdf(pos)
        
        return X, Y, Z
        
    def plot_terrain_2d(self, gene_terrain, show_genes=True, show_labels=False, 
                        show_edges=True, cmap='jet', ax=None, title=None,
                        label_fontsize=8):
        """
        Plot 2D heatmap of the terrain with network edges and gene labels.
        
        Parameters
        ----------
        gene_terrain : GeneTerrain
            GeneTerrain object with terrain data
        show_genes : bool
            Whether to show gene positions as points (default: True)
        show_labels : bool
            Whether to show gene labels with white outlines (default: False)
        show_edges : bool
            Whether to show network edges between genes (default: True)
        cmap : str
            Colormap name (default: 'jet' - red=upregulated, blue=downregulated)
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates new figure.
        title : str, optional
            Plot title
        label_fontsize : int
            Font size for gene labels (default: 8)
            
        Returns
        -------
        matplotlib.axes.Axes
        """
        if gene_terrain.terrain is None:
            raise ValueError("Terrain not created. Call create_terrain() first.")
        
        if ax is None:
            fig, ax = plt.subplots(figsize=self.figsize, dpi=self.dpi)
        
        # Plot high-resolution terrain with equal aspect ratio (square)
        im = ax.imshow(gene_terrain.terrain, cmap=cmap, aspect='equal',
                      extent=[0, 1, 0, 1], origin='lower', interpolation='bilinear')
        
        # Add colorbar
        plt.colorbar(im, ax=ax, label='Expression Level')
        
        # Plot network edges first (so they appear behind nodes)
        if show_edges and gene_terrain.network is not None and gene_terrain.gene_data is not None:
            gene_data = gene_terrain.gene_data
            for edge in gene_terrain.network.edges():
                gene1, gene2 = edge
                if gene1 in gene_data['gene'].values and gene2 in gene_data['gene'].values:
                    row1 = gene_data[gene_data['gene'] == gene1].iloc[0]
                    row2 = gene_data[gene_data['gene'] == gene2].iloc[0]
                    x1, y1 = row1['X'], row1['Y']
                    x2, y2 = row2['X'], row2['Y']
                    
                    # Draw thin black edge
                    ax.plot([x1, x2], [y1, y2], 'k-', linewidth=0.5, alpha=0.3, zorder=1)
        
        # Plot gene nodes
        if show_genes and gene_terrain.gene_data is not None:
            for _, gene in gene_terrain.gene_data.iterrows():
                ax.plot(gene['X'], gene['Y'], 'ko', markersize=6, 
                       markeredgecolor='white', markeredgewidth=1.5, zorder=2)
                
                if show_labels:
                    # Create text with white outline using path_effects
                    text = ax.text(gene['X'], gene['Y'], gene['gene'],
                                  fontsize=label_fontsize, ha='center', va='bottom',
                                  color='black', weight='bold', zorder=3)
                    text.set_path_effects([
                        path_effects.Stroke(linewidth=3, foreground='white'),
                        path_effects.Normal()
                    ])
        
        ax.set_xlabel('X Coordinate', fontsize=11)
        ax.set_ylabel('Y Coordinate', fontsize=11)
        ax.set_title(title or 'GeneTerrain 2D Heatmap', fontsize=13, weight='bold')
        ax.set_aspect('equal', adjustable='box')  # Square axes
        
        return ax
    
    def plot_terrain_3d(self, gene_terrain, cmap='jet', elev=30, azim=45,
                       title=None):
        """
        Plot 3D surface of the terrain.
        
        Parameters
        ----------
        gene_terrain : GeneTerrain
            GeneTerrain object
        cmap : str
            Colormap name
        elev : float
            Elevation angle for 3D view
        azim : float
            Azimuth angle for 3D view
        title : str, optional
            Plot title
            
        Returns
        -------
        matplotlib.axes.Axes
        """
        if gene_terrain.terrain is None:
            raise ValueError("Terrain not created.")
        
        fig = plt.figure(figsize=self.figsize, dpi=self.dpi)
        ax = fig.add_subplot(111, projection='3d')
        
        # Create coordinate meshes
        res = gene_terrain.resolution
        x = np.linspace(0, 1, res)
        y = np.linspace(0, 1, res)
        X, Y = np.meshgrid(x, y)
        
        # Plot surface
        surf = ax.plot_surface(X, Y, gene_terrain.terrain, cmap=cmap,
                              linewidth=0, antialiased=True, alpha=0.9)
        
        # Add colorbar
        fig.colorbar(surf, ax=ax, shrink=0.5, label='Expression Level')
        
        # Set view angle
        ax.view_init(elev=elev, azim=azim)
        
        ax.set_xlabel('X Coordinate')
        ax.set_ylabel('Y Coordinate')
        ax.set_zlabel('Expression')
        ax.set_title(title or 'GeneTerrain 3D Surface')
        
        return ax
    
    def plot_gmm_components(self, gene_terrain, gmm_optimizer, 
                           highlight_component=None,
                           cmap='jet', title=None,
                           contour_levels=[0.1, 0.3, 0.5, 0.7, 0.9]):
        """
        Plot terrain with TRUE GAUSSIAN probability density contours for GMM components.
        
        This function displays:
        - The complete GeneTerrain as a base layer
        - True Gaussian probability density contours (NOT argmax boundaries)
        - Multiple probability levels showing component shape and overlap
        - All genes with white-outlined black text labels
        - Network edges connecting genes (thin black lines)
        - Optional highlighting of a specific Gaussian component
        
        Parameters
        ----------
        gene_terrain : GeneTerrain
            GeneTerrain object with fitted terrain and network
        gmm_optimizer : GMMOptimizer
            Fitted GMMOptimizer with component assignments
        highlight_component : int, optional
            Specific GMM component index to highlight (0-based).
            If specified, this component is emphasized with thicker red contours
            while other components are shown with thinner dimmed contours.
        cmap : str
            Colormap for terrain (default: 'jet' - red=upregulated, blue=downregulated)
        title : str, optional
            Plot title
        contour_levels : list of float
            Relative probability levels for contours (default: [0.1, 0.3, 0.5, 0.7, 0.9])
            Higher values = tighter contours around component center
            
        Returns
        -------
        matplotlib.axes.Axes
            
        Notes
        -----
        This method plots TRUE Gaussian probability density contours, showing the actual
        shape and spread of each GMM component. This is scientifically more rigorous than
        showing argmax boundaries, as it:
        - Reveals component covariance structure (ellipse orientation/shape)
        - Shows overlapping probability regions between components
        - Enables identification of high-confidence vs uncertain regions
        - Displays the continuous probability surface rather than discrete boundaries
        """
        if gene_terrain.terrain is None or gmm_optimizer.labels is None:
            raise ValueError("Terrain and GMM must be fitted first.")
        
        fig, ax = plt.subplots(figsize=self.figsize, dpi=self.dpi)
        
        # Plot terrain as base with equal aspect ratio (square)
        im = ax.imshow(gene_terrain.terrain, cmap=cmap, aspect='equal',
                      extent=[0, 1, 0, 1], origin='lower', alpha=0.6)
        plt.colorbar(im, ax=ax, label='Expression Level')
        
        # Get GMM parameters
        gmm_model = gmm_optimizer.best_model
        n_components = gmm_optimizer.optimal_n_components
        
        # Create colormap for components
        colors = plt.cm.Set3(np.linspace(0, 1, n_components))
        
        # Plot TRUE GAUSSIAN PROBABILITY CONTOURS for each component
        for i in range(n_components):
            # Get mean and covariance for this component (3D: X, Y, Expression)
            mean_3d = gmm_model.means_[i]
            covariance_3d = gmm_model.covariances_[i]
            weight = gmm_model.weights_[i]
            
            # Extract 2D projection (X, Y coordinates only)
            mean_2d = mean_3d[:2]  # First two dimensions: X, Y
            covariance_2d = covariance_3d[:2, :2]  # Top-left 2x2 block
            
            # Compute probability density surface
            X, Y, Z = self._compute_gaussian_contours(mean_2d, covariance_2d, resolution=150)
            
            # Normalize Z to get relative probability levels (0 to 1)
            if Z.max() > 0:
                Z_norm = Z / Z.max()
            else:
                Z_norm = Z
            
            # Determine visualization style based on highlighting
            if highlight_component is not None:
                if i == highlight_component:
                    # PRIMARY HIGHLIGHT: Thick red contours with strong fill
                    contour_lines = ax.contour(X, Y, Z_norm, levels=contour_levels,
                                              colors='red', linewidths=3, 
                                              alpha=0.9, zorder=3)
                    ax.clabel(contour_lines, inline=True, fontsize=8, fmt='%0.1f')
                    
                    # Fill with component color
                    ax.contourf(X, Y, Z_norm, levels=contour_levels,
                               colors=[colors[i]], alpha=0.3, zorder=2)
                else:
                    # SECONDARY: Thin dimmed contours
                    ax.contour(X, Y, Z_norm, levels=contour_levels,
                              colors=[colors[i]], linewidths=1, 
                              alpha=0.3, zorder=2, linestyles='dashed')
            else:
                # ALL EQUAL: Medium thickness colored contours with fills
                contour_lines = ax.contour(X, Y, Z_norm, levels=contour_levels,
                                          colors=[colors[i]], linewidths=2,
                                          alpha=0.7, zorder=2)
                ax.contourf(X, Y, Z_norm, levels=contour_levels,
                           colors=[colors[i]], alpha=0.15, zorder=2)
        
        # Plot network edges first (behind genes)
        if gene_terrain.network is not None and gene_terrain.gene_data is not None:
            gene_pos = gene_terrain.gene_data.set_index('gene')[['X', 'Y']].to_dict('index')
            
            for edge in gene_terrain.network.edges():
                if edge[0] in gene_pos and edge[1] in gene_pos:
                    x_coords = [gene_pos[edge[0]]['X'], gene_pos[edge[1]]['X']]
                    y_coords = [gene_pos[edge[0]]['Y'], gene_pos[edge[1]]['Y']]
                    ax.plot(x_coords, y_coords, 'k-', linewidth=0.5, alpha=0.3, zorder=4)
        
        # Plot genes with white-outlined black text labels
        # Label genes within major component boundaries (probability > 0.5)
        if gene_terrain.gene_data is not None:
            for _, gene in gene_terrain.gene_data.iterrows():
                gene_pos = np.array([gene['X'], gene['Y']])
                
                # Determine if gene is within any component's boundary
                # Check probability at gene position for each component
                show_label = False
                max_prob = 0
                
                for i in range(n_components):
                    mean_3d = gmm_model.means_[i]
                    covariance_3d = gmm_model.covariances_[i]
                    mean_2d = mean_3d[:2]
                    covariance_2d = covariance_3d[:2, :2]
                    
                    # Compute probability at gene position
                    rv = multivariate_normal(mean_2d, covariance_2d, allow_singular=True)
                    prob = rv.pdf(gene_pos)
                    
                    # Get max probability for this component to normalize
                    max_prob_component = rv.pdf(mean_2d)
                    
                    # Normalize probability (0 to 1 relative to component center)
                    if max_prob_component > 0:
                        normalized_prob = prob / max_prob_component
                    else:
                        normalized_prob = 0
                    
                    # If gene is within 0.5 probability contour of any component, show label
                    if normalized_prob > 0.3:  # Show if within 30% probability threshold
                        show_label = True
                        max_prob = max(max_prob, normalized_prob)
                        break
                
                # When highlighting, always show labels for highlighted component genes
                if highlight_component is not None:
                    mean_3d = gmm_model.means_[highlight_component]
                    covariance_3d = gmm_model.covariances_[highlight_component]
                    mean_2d = mean_3d[:2]
                    covariance_2d = covariance_3d[:2, :2]
                    
                    rv = multivariate_normal(mean_2d, covariance_2d, allow_singular=True)
                    prob = rv.pdf(gene_pos)
                    max_prob_component = rv.pdf(mean_2d)
                    
                    if max_prob_component > 0:
                        normalized_prob = prob / max_prob_component
                        if normalized_prob > 0.2:  # Lower threshold for highlighted component
                            show_label = True
                
                # Plot gene node (always)
                ax.plot(gene['X'], gene['Y'], 'ko', markersize=6,
                       markeredgecolor='white', markeredgewidth=1, zorder=5)
                
                # Show label only if gene is within component boundary
                if show_label:
                    text = ax.text(gene['X'], gene['Y'], gene['gene'],
                                  fontsize=8, ha='center', va='bottom',
                                  color='black', weight='bold', zorder=6)
                    text.set_path_effects([
                        path_effects.Stroke(linewidth=3, foreground='white'),
                        path_effects.Normal()
                    ])
        
        ax.set_xlabel('X Coordinate', fontsize=11)
        ax.set_ylabel('Y Coordinate', fontsize=11)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_aspect('equal', adjustable='box')  # Square axes
        
        if title:
            ax.set_title(title, fontsize=13, weight='bold')
        elif highlight_component is not None:
            ax.set_title(f'GeneTerrain with Gaussian Component {highlight_component} (True Probability Contours)', 
                        fontsize=13, weight='bold')
        else:
            ax.set_title(f'GeneTerrain with {n_components} GMM Components (True Probability Contours)',
                        fontsize=13, weight='bold')
        
        # Add legend explaining contour levels
        legend_text = f'Contour levels: {contour_levels}\n(relative probability: 1.0 = component center)'
        ax.text(0.02, 0.98, legend_text, transform=ax.transAxes,
               fontsize=8, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        return ax
    
    def plot_optimization_metrics(self, gmm_optimizer, metrics='all'):
        """
        Plot GMM optimization metrics.
        
        Parameters
        ----------
        gmm_optimizer : GMMOptimizer
            GMMOptimizer with optimization results
        metrics : str or list
            Metrics to plot: 'all', 'clustering', 'information', or list of metric names
            
        Returns
        -------
        matplotlib.figure.Figure
        """
        if gmm_optimizer.optimization_results is None:
            raise ValueError("No optimization results available.")
        
        df = gmm_optimizer.optimization_results
        
        # Determine which metrics to plot
        if metrics == 'all':
            metric_names = ['silhouette', 'davies_bouldin', 'calinski_harabasz', 'bic', 'aic']
        elif metrics == 'clustering':
            metric_names = ['silhouette', 'davies_bouldin', 'calinski_harabasz']
        elif metrics == 'information':
            metric_names = ['bic', 'aic']
        else:
            metric_names = metrics
        
        # Create subplots
        n_metrics = len(metric_names)
        n_cols = min(3, n_metrics)
        n_rows = (n_metrics + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows), dpi=self.dpi)
        if n_metrics == 1:
            axes = [axes]
        else:
            axes = axes.flatten()
        
        # Plot each metric
        for i, metric in enumerate(metric_names):
            if metric not in df.columns:
                continue
            
            ax = axes[i]
            ax.plot(df['n_components'], df[metric], 'bo-', linewidth=2, markersize=8)
            
            # Mark optimal
            if gmm_optimizer.optimal_n_components is not None:
                optimal_val = df[df['n_components'] == gmm_optimizer.optimal_n_components][metric].values[0]
                ax.plot(gmm_optimizer.optimal_n_components, optimal_val, 
                       'r*', markersize=20, label='Optimal')
            
            ax.set_xlabel('Number of Components')
            ax.set_ylabel(metric.replace('_', ' ').title())
            ax.set_title(f'{metric.replace("_", " ").title()} vs Components')
            ax.grid(True, alpha=0.3)
            ax.legend()
        
        # Hide unused subplots
        for i in range(n_metrics, len(axes)):
            axes[i].axis('off')
        
        plt.tight_layout()
        return fig
    
    def plot_sigma_optimization(self, gene_terrain):
        """
        Plot sigma optimization results.
        
        Parameters
        ----------
        gene_terrain : GeneTerrain
            GeneTerrain with sigma_metrics
            
        Returns
        -------
        matplotlib.figure.Figure
        """
        if gene_terrain.sigma_metrics is None:
            raise ValueError("No sigma optimization results available.")
        
        df = gene_terrain.sigma_metrics
        
        fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=self.dpi)
        axes = axes.flatten()
        
        metrics = ['coverage', 'preservation', 'smoothness', 'dynamic_range', 'composite_score']
        
        for i, metric in enumerate(metrics):
            if metric not in df.columns:
                continue
            
            ax = axes[i]
            ax.plot(df['sigma'], df[metric], 'b-', linewidth=2)
            
            if gene_terrain.optimal_sigma is not None:
                optimal_val = df[df['sigma'] == gene_terrain.optimal_sigma][metric].values[0]
                ax.plot(gene_terrain.optimal_sigma, optimal_val, 
                       'r*', markersize=20, label=f'Optimal σ={gene_terrain.optimal_sigma:.3f}')
            
            ax.set_xlabel('Sigma')
            ax.set_ylabel(metric.replace('_', ' ').title())
            ax.set_title(f'{metric.replace("_", " ").title()} vs Sigma')
            ax.grid(True, alpha=0.3)
            ax.legend()
        
        # Plot optimal terrain in last subplot with square aspect and jet colormap
        if gene_terrain.terrain is not None:
            ax = axes[5]
            im = ax.imshow(gene_terrain.terrain, cmap='jet', aspect='equal',
                          extent=[0, 1, 0, 1], origin='lower', interpolation='bilinear')
            ax.set_title(f'Optimal Terrain (σ={gene_terrain.optimal_sigma:.3f})', 
                        fontsize=11, weight='bold')
            ax.set_xlabel('X Coordinate', fontsize=9)
            ax.set_ylabel('Y Coordinate', fontsize=9)
            ax.set_aspect('equal', adjustable='box')
            plt.colorbar(im, ax=ax, label='Expression')
        
        plt.tight_layout()
        return fig
    
    def plot_network(self, gene_terrain, node_size_by_expression=True,
                    layout='spring', title=None):
        """
        Plot the gene interaction network.
        
        Parameters
        ----------
        gene_terrain : GeneTerrain
            GeneTerrain with network
        node_size_by_expression : bool
            Scale node size by expression level
        layout : str
            Network layout: 'spring', 'kamada_kawai', or 'current'
        title : str, optional
            Plot title
            
        Returns
        -------
        matplotlib.axes.Axes
        """
        if gene_terrain.network is None:
            raise ValueError("Network not built.")
        
        fig, ax = plt.subplots(figsize=self.figsize, dpi=self.dpi)
        
        # Get positions
        if layout == 'current' and gene_terrain.gene_data is not None:
            pos = {row['gene']: (row['X'], row['Y']) 
                   for _, row in gene_terrain.gene_data.iterrows()}
        elif layout == 'spring':
            pos = nx.spring_layout(gene_terrain.network, seed=42)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(gene_terrain.network)
        else:
            pos = nx.spring_layout(gene_terrain.network, seed=42)
        
        # Node sizes
        if node_size_by_expression and gene_terrain.gene_data is not None:
            exp_dict = gene_terrain.gene_data.set_index('gene')['Exp'].to_dict()
            node_sizes = [abs(exp_dict.get(node, 1)) * 300 + 100 
                         for node in gene_terrain.network.nodes()]
        else:
            node_sizes = 300
        
        # Node colors
        if gene_terrain.gene_data is not None:
            exp_dict = gene_terrain.gene_data.set_index('gene')['Exp'].to_dict()
            node_colors = [exp_dict.get(node, 0) for node in gene_terrain.network.nodes()]
        else:
            node_colors = 'lightblue'
        
        # Draw network
        nx.draw_networkx_edges(gene_terrain.network, pos, ax=ax, alpha=0.3,
                              width=0.5, edge_color='gray')
        nx.draw_networkx_nodes(gene_terrain.network, pos, ax=ax,
                              node_size=node_sizes, node_color=node_colors,
                              cmap='jet', vmin=-2, vmax=2,
                              edgecolors='black', linewidths=0.5)
        nx.draw_networkx_labels(gene_terrain.network, pos, ax=ax, font_size=8)
        
        ax.set_title(title or 'Gene Interaction Network', fontsize=13, weight='bold')
        ax.axis('off')
        
        # Add colorbar with jet colormap (red=upregulated, blue=downregulated)
        sm = plt.cm.ScalarMappable(cmap='jet', 
                                   norm=mcolors.Normalize(vmin=-2, vmax=2))
        sm.set_array([])
        plt.colorbar(sm, ax=ax, label='Expression Level (log2FC)')
        
        return ax
