"""
GMM optimization module for GeneTerrain clustering.
"""

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
import warnings

warnings.filterwarnings('ignore', category=UserWarning)


class GMMOptimizer:
    """
    Optimize and fit Gaussian Mixture Models to GeneTerrain data.
    
    Parameters
    ----------
    max_components : int
        Maximum number of GMM components to test (default: 10)
    covariance_type : str
        Type of covariance: 'full', 'tied', 'diag', 'spherical' (default: 'full')
    random_state : int
        Random seed for reproducibility (default: 42)
    """
    
    def __init__(self, max_components=10, covariance_type='full', random_state=42):
        self.max_components = max_components
        self.covariance_type = covariance_type
        self.random_state = random_state
        
        self.optimization_results = None
        self.optimal_n_components = None
        self.best_model = None
        self.labels = None
        self.probabilities = None
        
    def optimize(self, data, min_components=2):
        """
        Find optimal number of GMM components using silhouette score.
        
        Parameters
        ----------
        data : np.ndarray
            Data array of shape (n_samples, n_features)
        min_components : int
            Minimum number of components to test (default: 2)
            
        Returns
        -------
        int
            Optimal number of components
        """
        print(f"\nOptimizing GMM components (range: {min_components}-{self.max_components})...")
        print("=" * 70)
        
        results = {
            'n_components': [],
            'silhouette': [],
            'davies_bouldin': [],
            'calinski_harabasz': [],
            'bic': [],
            'aic': []
        }
        
        best_silhouette = -1
        best_n = min_components
        models = {}
        
        for n in range(min_components, self.max_components + 1):
            print(f"\nTesting {n} components...")
            
            # Fit GMM
            gmm = GaussianMixture(
                n_components=n,
                covariance_type=self.covariance_type,
                random_state=self.random_state,
                n_init=5,
                max_iter=200
            )
            
            gmm.fit(data)
            models[n] = gmm
            
            # Get labels
            labels = gmm.predict(data)
            
            # Calculate metrics
            if len(np.unique(labels)) > 1:
                # Silhouette (primary metric)
                silhouette = silhouette_score(
                    data, labels, 
                    sample_size=min(5000, len(data))
                )
                
                # Other metrics
                davies_bouldin = davies_bouldin_score(data, labels)
                calinski_harabasz = calinski_harabasz_score(data, labels)
                
                # Track best
                if silhouette > best_silhouette:
                    best_silhouette = silhouette
                    best_n = n
            else:
                silhouette = -1
                davies_bouldin = float('inf')
                calinski_harabasz = 0
            
            # Store results
            results['n_components'].append(n)
            results['silhouette'].append(silhouette)
            results['davies_bouldin'].append(davies_bouldin)
            results['calinski_harabasz'].append(calinski_harabasz)
            results['bic'].append(gmm.bic(data))
            results['aic'].append(gmm.aic(data))
            
            # Print progress
            marker = " <-- BEST" if n == best_n else ""
            print(f"  Silhouette: {silhouette:.4f}{marker}")
            print(f"  Davies-Bouldin: {davies_bouldin:.4f} (lower is better)")
            print(f"  Calinski-Harabasz: {calinski_harabasz:.2f} (higher is better)")
            print(f"  Converged: {gmm.converged_}")
        
        # Store results
        self.optimization_results = pd.DataFrame(results)
        self.optimal_n_components = best_n
        self.best_model = models[best_n]
        
        # Fit final model
        self.labels = self.best_model.predict(data)
        self.probabilities = self.best_model.predict_proba(data)
        
        print(f"\n{'=' * 70}")
        print(f"Optimal number of components: {best_n}")
        print(f"Best silhouette score: {best_silhouette:.4f}")
        
        return best_n
    
    def fit(self, data, n_components=None):
        """
        Fit GMM with specified or optimal number of components.
        
        Parameters
        ----------
        data : np.ndarray
            Data to fit
        n_components : int, optional
            Number of components. If None, uses optimized value.
            
        Returns
        -------
        GMMOptimizer
            Self for method chaining
        """
        if n_components is None:
            if self.optimal_n_components is None:
                self.optimize(data)
            n_components = self.optimal_n_components
        
        print(f"\nFitting GMM with {n_components} components...")
        
        gmm = GaussianMixture(
            n_components=n_components,
            covariance_type=self.covariance_type,
            random_state=self.random_state,
            n_init=10,
            max_iter=300
        )
        
        gmm.fit(data)
        
        self.best_model = gmm
        self.optimal_n_components = n_components
        self.labels = gmm.predict(data)
        self.probabilities = gmm.predict_proba(data)
        
        print(f"GMM fitted successfully (converged: {gmm.converged_})")
        
        return self
    
    def get_component_genes(self, gene_terrain, component_idx, threshold=0.5):
        """
        Get genes belonging to a specific GMM component.
        
        Parameters
        ----------
        gene_terrain : GeneTerrain
            GeneTerrain object with gene_data
        component_idx : int
            Component index (0-based)
        threshold : float
            Probability threshold for assignment (default: 0.5)
            
        Returns
        -------
        pd.DataFrame
            DataFrame with genes, coordinates, expressions, and probabilities
        """
        if self.probabilities is None:
            raise ValueError("Model not fitted. Call fit() first.")
        
        if gene_terrain.gene_data is None:
            raise ValueError("GeneTerrain has no gene data.")
        
        # Get terrain data
        terrain_data = gene_terrain.get_terrain_data()
        
        # For each gene, find its probability in the component
        gene_assignments = []
        
        for _, gene_row in gene_terrain.gene_data.iterrows():
            # Find closest terrain point
            x_idx = int(gene_row['X'] * (gene_terrain.resolution - 1))
            y_idx = int(gene_row['Y'] * (gene_terrain.resolution - 1))
            flat_idx = y_idx * gene_terrain.resolution + x_idx
            
            # Get probability
            prob = self.probabilities[flat_idx, component_idx]
            
            if prob >= threshold:
                gene_assignments.append({
                    'gene': gene_row['gene'],
                    'X': gene_row['X'],
                    'Y': gene_row['Y'],
                    'expression': gene_row['Exp'],
                    'probability': prob,
                    'component': component_idx
                })
        
        df = pd.DataFrame(gene_assignments)
        
        if len(df) > 0:
            df = df.sort_values('probability', ascending=False)
        
        return df
    
    def get_all_component_genes(self, gene_terrain, threshold=0.5):
        """
        Get genes for all components.
        
        Parameters
        ----------
        gene_terrain : GeneTerrain
            GeneTerrain object with gene_data
        threshold : float
            Probability threshold for assignment
            
        Returns
        -------
        dict
            Dictionary mapping component_idx -> DataFrame of genes
        """
        if self.optimal_n_components is None:
            raise ValueError("Model not fitted.")
        
        result = {}
        for i in range(self.optimal_n_components):
            result[i] = self.get_component_genes(gene_terrain, i, threshold)
        
        return result
    
    def get_component_summary(self):
        """
        Get summary statistics for each component.
        
        Returns
        -------
        pd.DataFrame
            Summary with component weights and sizes
        """
        if self.best_model is None:
            raise ValueError("Model not fitted.")
        
        weights = self.best_model.weights_
        
        # Count assignments (hard assignment)
        unique, counts = np.unique(self.labels, return_counts=True)
        
        summary = []
        for i in range(self.optimal_n_components):
            count = counts[i] if i in unique else 0
            summary.append({
                'component': i,
                'weight': weights[i],
                'n_points': count,
                'percentage': count / len(self.labels) * 100
            })
        
        df = pd.DataFrame(summary)
        df = df.sort_values('weight', ascending=False)
        
        return df
