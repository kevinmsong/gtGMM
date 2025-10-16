"""
Pathway Enrichment Analysis Module for gtGMM

This module provides tools for performing pathway enrichment analysis
on GMM components and gene sets.
"""

import numpy as np
import pandas as pd
from scipy import stats
from typing import Dict, List, Set, Tuple, Optional


class EnrichmentAnalyzer:
    """
    Pathway enrichment analyzer for gene sets and GMM components.
    
    This class performs statistical enrichment analysis using Fisher's exact test
    and provides FDR correction for multiple hypothesis testing.
    
    Parameters
    ----------
    pathways : dict
        Dictionary mapping pathway names to lists of gene symbols
    background_genes : list, optional
        List of all genes in the background set. If None, uses union of all pathway genes.
    """
    
    def __init__(self, pathways: Dict[str, List[str]], background_genes: Optional[List[str]] = None):
        self.pathways = pathways
        
        # Set background genes
        if background_genes is None:
            all_pathway_genes = set()
            for genes in pathways.values():
                all_pathway_genes.update(genes)
            self.background_genes = list(all_pathway_genes)
        else:
            self.background_genes = background_genes
        
        self.background_size = len(self.background_genes)
    
    def analyze_gene_set(self, gene_set: List[str], pathway_name: str = None) -> Dict:
        """
        Perform enrichment analysis for a single gene set against one or all pathways.
        
        Parameters
        ----------
        gene_set : list
            List of gene symbols to test for enrichment
        pathway_name : str, optional
            Specific pathway to test. If None, tests all pathways.
            
        Returns
        -------
        dict or pd.DataFrame
            Enrichment results with p-values and enrichment scores
        """
        gene_set = set(gene_set)
        
        if pathway_name is not None:
            # Test single pathway
            return self._test_pathway(gene_set, pathway_name)
        else:
            # Test all pathways
            results = []
            for pw_name in self.pathways:
                result = self._test_pathway(gene_set, pw_name)
                results.append(result)
            
            df = pd.DataFrame(results)
            
            # Apply FDR correction
            if len(df) > 0:
                df = self._apply_fdr_correction(df)
            
            return df.sort_values('p_value')
    
    def analyze_components(self, components: Dict[int, pd.DataFrame], 
                          condition: str = "Condition") -> pd.DataFrame:
        """
        Perform enrichment analysis for multiple GMM components.
        
        Parameters
        ----------
        components : dict
            Dictionary mapping component indices to DataFrames with 'gene' column
        condition : str
            Name of the condition (e.g., 'P1', 'P8')
            
        Returns
        -------
        pd.DataFrame
            Enrichment results for all component-pathway pairs
        """
        results = []
        
        for comp_idx, comp_df in components.items():
            comp_genes = set(comp_df['gene'].tolist())
            
            for pw_name, pw_genes in self.pathways.items():
                pw_set = set(pw_genes)
                overlap = len(comp_genes & pw_set)
                
                # Calculate enrichment score
                enrichment_score = overlap / len(pw_set) if len(pw_set) > 0 else 0
                
                # Fisher's exact test
                a = overlap  # in both
                b = len(comp_genes) - overlap  # in component, not pathway
                c = len(pw_set) - overlap  # in pathway, not component
                d = self.background_size - (a + b + c)  # in neither
                
                # Ensure all values are non-negative
                if a < 0 or b < 0 or c < 0 or d < 0:
                    p_value = 1.0
                elif a + b > 0 and c + d > 0:
                    _, p_value = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
                else:
                    p_value = 1.0
                
                results.append({
                    'condition': condition,
                    'component': comp_idx + 1,
                    'pathway': pw_name,
                    'overlap': overlap,
                    'pathway_size': len(pw_set),
                    'component_size': len(comp_genes),
                    'enrichment_score': enrichment_score,
                    'p_value': p_value
                })
        
        df = pd.DataFrame(results)
        
        # Apply FDR correction
        if len(df) > 0:
            df = self._apply_fdr_correction(df)
        
        return df
    
    def compare_conditions(self, components_a: Dict, components_b: Dict,
                          condition_a: str = "A", condition_b: str = "B") -> pd.DataFrame:
        """
        Compare pathway enrichment between two conditions.
        
        Parameters
        ----------
        components_a : dict
            Components from condition A
        components_b : dict
            Components from condition B
        condition_a : str
            Name of condition A
        condition_b : str
            Name of condition B
            
        Returns
        -------
        pd.DataFrame
            Combined enrichment results for both conditions
        """
        enrich_a = self.analyze_components(components_a, condition_a)
        enrich_b = self.analyze_components(components_b, condition_b)
        
        combined = pd.concat([enrich_a, enrich_b], ignore_index=True)
        
        return combined.sort_values(['pathway', 'condition', 'p_value'])
    
    def _test_pathway(self, gene_set: Set[str], pathway_name: str) -> Dict:
        """
        Test enrichment for a single pathway.
        
        Parameters
        ----------
        gene_set : set
            Set of genes to test
        pathway_name : str
            Name of the pathway
            
        Returns
        -------
        dict
            Enrichment statistics
        """
        pw_genes = set(self.pathways[pathway_name])
        overlap = len(gene_set & pw_genes)
        
        # Calculate enrichment score
        enrichment_score = overlap / len(pw_genes) if len(pw_genes) > 0 else 0
        
        # Fisher's exact test
        a = overlap
        b = len(gene_set) - overlap
        c = len(pw_genes) - overlap
        d = self.background_size - (a + b + c)
        
        if a + b > 0 and c + d > 0:
            _, p_value = stats.fisher_exact([[a, b], [c, d]], alternative='greater')
        else:
            p_value = 1.0
        
        return {
            'pathway': pathway_name,
            'overlap': overlap,
            'pathway_size': len(pw_genes),
            'gene_set_size': len(gene_set),
            'enrichment_score': enrichment_score,
            'p_value': p_value
        }
    
    def _apply_fdr_correction(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Apply Benjamini-Hochberg FDR correction.
        
        Parameters
        ----------
        df : pd.DataFrame
            DataFrame with 'p_value' column
            
        Returns
        -------
        pd.DataFrame
            DataFrame with added 'q_value' column
        """
        if 'p_value' not in df.columns or len(df) == 0:
            return df
        
        # Sort by p-value
        df_sorted = df.sort_values('p_value').copy()
        
        # Calculate q-values
        n = len(df_sorted)
        ranks = np.arange(1, n + 1)
        q_values = df_sorted['p_value'].values * n / ranks
        
        # Ensure monotonicity
        q_values = np.minimum.accumulate(q_values[::-1])[::-1]
        q_values = np.minimum(q_values, 1.0)
        
        df_sorted['q_value'] = q_values
        
        # Restore original order
        df_result = df_sorted.sort_index()
        
        return df_result
    
    def get_significant_enrichments(self, enrichment_df: pd.DataFrame, 
                                   p_threshold: float = 0.05,
                                   q_threshold: float = 0.1) -> pd.DataFrame:
        """
        Filter enrichment results to significant hits.
        
        Parameters
        ----------
        enrichment_df : pd.DataFrame
            Enrichment results from analyze_* methods
        p_threshold : float
            P-value threshold
        q_threshold : float
            Q-value (FDR) threshold
            
        Returns
        -------
        pd.DataFrame
            Filtered significant enrichments
        """
        if 'q_value' in enrichment_df.columns:
            return enrichment_df[
                (enrichment_df['p_value'] < p_threshold) & 
                (enrichment_df['q_value'] < q_threshold)
            ]
        else:
            return enrichment_df[enrichment_df['p_value'] < p_threshold]
    
    def get_pathway_genes(self, pathway_name: str) -> List[str]:
        """
        Get genes for a specific pathway.
        
        Parameters
        ----------
        pathway_name : str
            Name of the pathway
            
        Returns
        -------
        list
            List of gene symbols in the pathway
        """
        return self.pathways.get(pathway_name, [])
    
    def get_pathway_names(self) -> List[str]:
        """
        Get all pathway names.
        
        Returns
        -------
        list
            List of pathway names
        """
        return list(self.pathways.keys())

