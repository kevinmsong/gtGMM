"""
Core GeneTerrain class for gene network analysis.
"""

import numpy as np
import pandas as pd
import requests
import networkx as nx
from sklearn.preprocessing import MinMaxScaler
from scipy.ndimage import gaussian_filter
import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning)


class GeneTerrain:
    """
    GeneTerrain: Generate topological landscapes from gene expression data.
    
    This class handles fetching protein-protein interactions from STRING,
    computing network layouts, and generating 3D expression terrains.
    
    Parameters
    ----------
    seed_genes : list of str
        List of gene symbols to analyze
    expression_values : dict or pd.Series, optional
        Gene expression values {gene: expression}. If None, uses ones.
    species : int, optional
        NCBI taxonomy ID (default: 9606 for human)
    resolution : int, optional
        Terrain grid resolution for high-quality rendering (default: 1000)
    """
    
    def __init__(self, seed_genes, expression_values=None, species=9606, resolution=1000):
        self.seed_genes = list(seed_genes)
        self.species = species
        self.resolution = resolution
        
        # Handle expression values
        if expression_values is None:
            self.expression_values = {gene: 1.0 for gene in seed_genes}
        elif isinstance(expression_values, dict):
            self.expression_values = expression_values
        elif isinstance(expression_values, pd.Series):
            self.expression_values = expression_values.to_dict()
        else:
            raise ValueError("expression_values must be dict or pd.Series")
        
        # Initialize attributes
        self.interactions = None
        self.network = None
        self.gene_data = None
        self.terrain = None
        self.optimal_sigma = None
        self.sigma_metrics = None
        
    def fetch_string_interactions(self, score_threshold=0.0):
        """
        Fetch protein-protein interactions from STRING database.
        
        Parameters
        ----------
        score_threshold : float
            Minimum confidence score (0-1000 scale)
            
        Returns
        -------
        pd.DataFrame
            Interactions with columns: gene_a, gene_b, score
        """
        print(f"Fetching STRING interactions for {len(self.seed_genes)} genes...")
        
        try:
            string_api_url = "https://string-db.org/api"
            identifiers = "%0d".join(self.seed_genes)
            request_url = (f"{string_api_url}/tsv/network?"
                          f"identifiers={identifiers}&species={self.species}")
            
            response = requests.get(request_url, timeout=30)
            response.raise_for_status()
            
            # Parse response
            lines = response.text.splitlines()
            if len(lines) < 2:
                print("Warning: No interactions found")
                return pd.DataFrame(columns=['gene_a', 'gene_b', 'score'])
            
            data = [line.split("\t") for line in lines[1:]]
            columns = lines[0].split("\t")
            
            df = pd.DataFrame(data, columns=columns)
            
            # Extract relevant columns
            df = df[["preferredName_A", "preferredName_B", "score"]].copy()
            df.columns = ['gene_a', 'gene_b', 'score']
            df['score'] = df['score'].astype(float)
            
            # Filter to seed genes only
            df = df[
                df['gene_a'].isin(self.seed_genes) & 
                df['gene_b'].isin(self.seed_genes)
            ]
            
            # Apply score threshold
            df = df[df['score'] >= score_threshold]
            
            self.interactions = df
            print(f"Retrieved {len(df)} interactions")
            return df
            
        except Exception as e:
            print(f"Error fetching STRING interactions: {e}")
            return pd.DataFrame(columns=['gene_a', 'gene_b', 'score'])
    
    def build_network(self):
        """
        Build NetworkX graph from interactions and compute layout.
        
        Returns
        -------
        nx.Graph
            Gene interaction network
        """
        if self.interactions is None or len(self.interactions) == 0:
            print("No interactions available. Using complete graph.")
            # Create a complete graph if no interactions
            self.network = nx.complete_graph(self.seed_genes)
        else:
            # Build network from interactions
            self.network = nx.Graph()
            
            # Add all seed genes as nodes
            self.network.add_nodes_from(self.seed_genes)
            
            # Add edges with weights
            for _, row in self.interactions.iterrows():
                self.network.add_edge(
                    row['gene_a'], 
                    row['gene_b'],
                    weight=row['score'] / 1000.0  # Normalize to [0,1]
                )
        
        print(f"Network: {self.network.number_of_nodes()} nodes, "
              f"{self.network.number_of_edges()} edges")
        
        return self.network
    
    def compute_layout(self, layout_type='spring', seed=42):
        """
        Compute 2D coordinates for genes using network layout algorithm.
        
        Parameters
        ----------
        layout_type : str
            Layout algorithm: 'spring', 'kamada_kawai', or 'spectral'
        seed : int
            Random seed for reproducibility
            
        Returns
        -------
        pd.DataFrame
            Gene data with columns: gene, X, Y, Exp
        """
        if self.network is None:
            self.build_network()
        
        print(f"Computing {layout_type} layout...")
        
        # Compute layout
        if layout_type == 'spring':
            pos = nx.spring_layout(self.network, seed=seed, iterations=50)
        elif layout_type == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(self.network)
        elif layout_type == 'spectral':
            pos = nx.spectral_layout(self.network)
        else:
            raise ValueError(f"Unknown layout type: {layout_type}")
        
        # Normalize coordinates to [0, 1]
        coords = np.array(list(pos.values()))
        scaler = MinMaxScaler()
        coords_norm = scaler.fit_transform(coords)
        
        # Create gene data DataFrame
        genes = list(pos.keys())
        data = []
        for i, gene in enumerate(genes):
            exp = self.expression_values.get(gene, 0.0)
            data.append({
                'gene': gene,
                'X': coords_norm[i, 0],
                'Y': coords_norm[i, 1],
                'Exp': exp
            })
        
        self.gene_data = pd.DataFrame(data)
        print(f"Layout computed for {len(self.gene_data)} genes")
        
        return self.gene_data
    
    def create_terrain(self, base_sigma=0.1, smooth_sigma=1.0):
        """
        Generate 3D terrain from gene positions and expression values.
        
        Parameters
        ----------
        base_sigma : float
            Base Gaussian kernel width for gene contributions
        smooth_sigma : float
            Additional Gaussian smoothing parameter
            
        Returns
        -------
        np.ndarray
            2D terrain array normalized to [-1, 1]
        """
        if self.gene_data is None:
            self.compute_layout()
        
        print(f"Creating terrain with sigma={base_sigma:.4f}...")
        
        res = self.resolution
        base_sigma = max(base_sigma, 1e-5)  # Prevent numerical issues
        
        # Create coordinate grids
        x_ = np.linspace(0, 1, res)
        y_ = np.linspace(0, 1, res)
        X, Y = np.meshgrid(x_, y_)
        
        terrain = np.zeros((res, res))
        
        # Add Gaussian contributions from each gene
        for _, row in self.gene_data.iterrows():
            x, y = row["X"], row["Y"]
            expression = row["Exp"]
            
            # Gaussian kernel
            dist_squared = ((X - x) ** 2 + (Y - y) ** 2) / (base_sigma ** 2)
            gaussian = np.exp(-0.5 * dist_squared)
            
            terrain += expression * gaussian
        
        # Apply additional smoothing
        terrain = gaussian_filter(terrain, sigma=smooth_sigma)
        
        # Normalize to [-1, 1]
        terrain_min, terrain_max = terrain.min(), terrain.max()
        if terrain_max > terrain_min:
            terrain = -1 + 2 * (terrain - terrain_min) / (terrain_max - terrain_min)
        else:
            terrain = np.zeros_like(terrain)
        
        self.terrain = terrain
        print("Terrain created successfully")
        
        return terrain
    
    def optimize_sigma(self, sigma_range=None, n_samples=30):
        """
        Find optimal base_sigma using composite scoring.
        
        Parameters
        ----------
        sigma_range : tuple, optional
            (min_sigma, max_sigma) range to search
        n_samples : int
            Number of sigma values to test
            
        Returns
        -------
        float
            Optimal sigma value
        """
        if self.gene_data is None:
            self.compute_layout()
        
        if sigma_range is None:
            sigma_range = (0.01, 0.5)
        
        print(f"\nOptimizing sigma in range {sigma_range}...")
        
        sigmas = np.linspace(sigma_range[0], sigma_range[1], n_samples)
        results = []
        
        for sigma in sigmas:
            terrain = self.create_terrain(base_sigma=sigma)
            metrics = self._calculate_terrain_metrics(terrain)
            metrics['sigma'] = sigma
            results.append(metrics)
            
            if len(results) % 10 == 0:
                print(f"  Tested {len(results)}/{n_samples} sigma values...")
        
        # Create results DataFrame
        df = pd.DataFrame(results)
        
        # Normalize metrics
        for col in df.columns:
            if col != 'sigma':
                col_min, col_max = df[col].min(), df[col].max()
                if col_max > col_min:
                    df[f'{col}_norm'] = (df[col] - col_min) / (col_max - col_min)
                else:
                    df[f'{col}_norm'] = 0.5
        
        # Composite score
        df['composite_score'] = (
            0.25 * df['coverage_norm'] +
            0.25 * df['preservation_norm'] +
            0.25 * df['smoothness_norm'] +
            0.25 * df['dynamic_range_norm']
        )
        
        # Find optimal
        optimal_idx = df['composite_score'].idxmax()
        optimal_sigma = df.loc[optimal_idx, 'sigma']
        
        self.optimal_sigma = optimal_sigma
        self.sigma_metrics = df
        
        print(f"\nOptimal sigma: {optimal_sigma:.4f}")
        
        # Regenerate terrain with optimal sigma
        self.create_terrain(base_sigma=optimal_sigma)
        
        return optimal_sigma
    
    def _calculate_terrain_metrics(self, terrain):
        """Calculate quality metrics for a terrain."""
        res = self.resolution
        metrics = {}
        
        # 1. Gene coverage
        covered = 0
        for _, row in self.gene_data.iterrows():
            x_idx = int(np.clip(row["X"] * (res - 1), 0, res - 1))
            y_idx = int(np.clip(row["Y"] * (res - 1), 0, res - 1))
            if np.abs(terrain[y_idx, x_idx]) > 0.1:
                covered += 1
        metrics['coverage'] = covered / len(self.gene_data)
        
        # 2. Expression preservation
        preservation_scores = []
        for _, row in self.gene_data.iterrows():
            x_idx = int(np.clip(row["X"] * (res - 1), 0, res - 1))
            y_idx = int(np.clip(row["Y"] * (res - 1), 0, res - 1))
            
            max_exp = self.gene_data["Exp"].abs().max()
            if max_exp > 0:
                norm_exp = row["Exp"] / max_exp
                norm_terrain = terrain[y_idx, x_idx]
                preservation_scores.append(1 - abs(norm_exp - norm_terrain))
        
        metrics['preservation'] = np.mean(preservation_scores) if preservation_scores else 0
        
        # 3. Smoothness
        grad_x, grad_y = np.gradient(terrain)
        gradient_mag = np.sqrt(grad_x**2 + grad_y**2)
        metrics['smoothness'] = 1.0 / (1.0 + np.mean(gradient_mag))
        
        # 4. Dynamic range
        metrics['dynamic_range'] = (terrain.max() - terrain.min()) / 2.0
        
        return metrics
    
    def get_terrain_data(self):
        """
        Get terrain data suitable for GMM fitting.
        
        Returns
        -------
        np.ndarray
            Array of shape (n_points, 3) with columns [X, Y, terrain_value]
        """
        if self.terrain is None:
            raise ValueError("Terrain not yet created. Call create_terrain() first.")
        
        res = self.resolution
        x_ = np.linspace(0, 1, res)
        y_ = np.linspace(0, 1, res)
        X, Y = np.meshgrid(x_, y_)
        
        # Flatten and stack
        data = np.column_stack([
            X.flatten(),
            Y.flatten(),
            self.terrain.flatten()
        ])
        
        return data
