"""
GEO (Gene Expression Omnibus) data integration module.

This module provides functions to fetch gene expression data from NCBI GEO
using the GEOparse library and direct API calls.
"""

import pandas as pd
import numpy as np
import requests
import re
from io import StringIO
import warnings


def fetch_geo_dataset(geo_id, platform=None, use_cached=True):
    """
    Fetch gene expression data from GEO database.
    
    Parameters
    ----------
    geo_id : str
        GEO dataset ID (e.g., 'GDS5093' or 'GSE48350')
    platform : str, optional
        Platform ID to use if multiple platforms available
    use_cached : bool
        Whether to use cached data if available
        
    Returns
    -------
    tuple
        (expression_df, metadata_df) where expression_df has genes as index
        and samples as columns
        
    Raises
    ------
    ImportError
        If GEOparse is not installed
    ValueError
        If dataset cannot be fetched
    """
    try:
        import GEOparse
    except ImportError:
        raise ImportError(
            "GEOparse is required for GEO integration. "
            "Install with: pip install GEOparse"
        )
    
    print(f"Fetching GEO dataset {geo_id}...")
    
    try:
        if geo_id.startswith('GDS'):
            # Dataset
            gds = GEOparse.get_GEO(geo=geo_id, destdir="./geo_cache" if use_cached else None)
            expression_df = gds.table
            metadata_df = gds.metadata
            
        elif geo_id.startswith('GSE'):
            # Series
            gse = GEOparse.get_GEO(geo=geo_id, destdir="./geo_cache" if use_cached else None)
            
            # Get platform
            if platform is None:
                platform = list(gse.gpls.keys())[0]
                print(f"Using platform: {platform}")
            
            # Extract expression data
            expression_dfs = []
            for gsm_name, gsm in gse.gsms.items():
                if gsm.metadata.get('platform_id', [None])[0] == platform:
                    expression_dfs.append(gsm.table.set_index('ID_REF'))
            
            expression_df = pd.concat(expression_dfs, axis=1)
            metadata_df = pd.DataFrame({k: gse.metadata[k] for k in gse.metadata})
            
        else:
            raise ValueError(f"Unsupported GEO ID format: {geo_id}")
        
        print(f"Successfully fetched {expression_df.shape[0]} features x {expression_df.shape[1]} samples")
        return expression_df, metadata_df
        
    except Exception as e:
        raise ValueError(f"Failed to fetch GEO dataset {geo_id}: {e}")


def fetch_geo_series_matrix(gse_id, use_cached=True):
    """
    Fetch GEO series matrix file (faster alternative to GEOparse).
    
    Parameters
    ----------
    gse_id : str
        GEO series ID (e.g., 'GSE48350')
    use_cached : bool
        Whether to cache downloaded files
        
    Returns
    -------
    tuple
        (expression_df, sample_metadata)
    """
    print(f"Fetching GEO series matrix for {gse_id}...")
    
    # Construct URL
    gse_short = gse_id.replace('GSE', '')
    base_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
    
    # Determine folder (e.g., GSE48350 -> GSE48nnn)
    folder = f"GSE{gse_short[:-3]}nnn/{gse_id}/matrix/"
    url = base_url + folder + f"{gse_id}_series_matrix.txt.gz"
    
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        
        # Parse series matrix
        import gzip
        content = gzip.decompress(response.content).decode('utf-8')
        
        # Extract sample info and expression data
        lines = content.split('\n')
        
        # Find data table start
        data_start = None
        for i, line in enumerate(lines):
            if line.startswith('!series_matrix_table_begin'):
                data_start = i + 1
                break
        
        if data_start is None:
            raise ValueError("Could not find data table in series matrix")
        
        # Find data table end
        data_end = None
        for i, line in enumerate(lines[data_start:], start=data_start):
            if line.startswith('!series_matrix_table_end'):
                data_end = i
                break
        
        # Read expression data
        data_lines = lines[data_start:data_end]
        expression_df = pd.read_csv(StringIO('\n'.join(data_lines)), sep='\t', index_col=0)
        
        # Extract metadata
        metadata = {}
        for line in lines[:data_start]:
            if line.startswith('!Sample_'):
                match = re.match(r'!Sample_(\w+)\s+"(.+)"', line)
                if match:
                    key, values = match.groups()
                    metadata[key] = values.split('\t')
        
        sample_metadata = pd.DataFrame(metadata)
        
        print(f"Successfully fetched {len(expression_df)} features x {len(expression_df.columns)} samples")
        return expression_df, sample_metadata
        
    except Exception as e:
        print(f"Failed to fetch series matrix: {e}")
        print("Falling back to GEOparse...")
        return fetch_geo_dataset(gse_id, use_cached=use_cached)


def extract_gene_expression(expression_df, gene_list, aggregate='mean'):
    """
    Extract and aggregate expression values for specific genes.
    
    Parameters
    ----------
    expression_df : pd.DataFrame
        Expression data with genes/probes as index
    gene_list : list
        List of gene symbols to extract
    aggregate : str
        How to aggregate multiple probes: 'mean', 'median', 'max', 'first'
        
    Returns
    -------
    pd.Series
        Expression values for requested genes
    """
    # Try direct gene symbol match
    matched_genes = {}
    
    for gene in gene_list:
        # Case-insensitive matching
        mask = expression_df.index.str.upper() == gene.upper()
        
        if mask.any():
            values = expression_df.loc[mask].values
            
            # Aggregate across samples
            sample_values = np.nanmean(values, axis=1)
            
            # Aggregate across probes
            if aggregate == 'mean':
                matched_genes[gene] = np.nanmean(sample_values)
            elif aggregate == 'median':
                matched_genes[gene] = np.nanmedian(sample_values)
            elif aggregate == 'max':
                matched_genes[gene] = np.nanmax(np.abs(sample_values))
            elif aggregate == 'first':
                matched_genes[gene] = sample_values[0]
            else:
                matched_genes[gene] = np.nanmean(sample_values)
    
    if len(matched_genes) == 0:
        warnings.warn(f"No genes matched in expression data")
    else:
        print(f"Matched {len(matched_genes)}/{len(gene_list)} genes")
    
    return pd.Series(matched_genes)


def compare_samples(expression_df, sample_indices, comparison='ratio'):
    """
    Compare expression between sample groups.
    
    Parameters
    ----------
    expression_df : pd.DataFrame
        Expression data
    sample_indices : tuple
        (control_indices, treatment_indices) as lists of column names or integers
    comparison : str
        Comparison method: 'ratio', 'difference', 'log2fc'
        
    Returns
    -------
    pd.Series
        Comparison values per gene
    """
    control_idx, treatment_idx = sample_indices
    
    control_data = expression_df.iloc[:, control_idx] if isinstance(control_idx[0], int) else expression_df[control_idx]
    treatment_data = expression_df.iloc[:, treatment_idx] if isinstance(treatment_idx[0], int) else expression_df[treatment_idx]
    
    control_mean = control_data.mean(axis=1)
    treatment_mean = treatment_data.mean(axis=1)
    
    if comparison == 'ratio':
        return treatment_mean / (control_mean + 1e-10)
    elif comparison == 'difference':
        return treatment_mean - control_mean
    elif comparison == 'log2fc':
        return np.log2((treatment_mean + 1) / (control_mean + 1))
    else:
        return treatment_mean - control_mean


def fetch_geo_for_geneterrain(geo_id, gene_list, sample_group=None, comparison=None):
    """
    Convenience function to fetch GEO data formatted for GeneTerrain.
    
    Parameters
    ----------
    geo_id : str
        GEO dataset ID
    gene_list : list
        List of genes to extract
    sample_group : int, list, or tuple, optional
        Sample indices to use. If tuple, compares (control, treatment)
    comparison : str, optional
        Comparison method if sample_group is tuple
        
    Returns
    -------
    dict
        Gene expression values ready for GeneTerrain
        
    Examples
    --------
    >>> # Single sample or average
    >>> expression = fetch_geo_for_geneterrain('GSE48350', genes)
    
    >>> # Specific samples
    >>> expression = fetch_geo_for_geneterrain('GSE48350', genes, sample_group=[0, 1, 2])
    
    >>> # Compare groups
    >>> expression = fetch_geo_for_geneterrain('GSE48350', genes, 
    ...     sample_group=([0,1,2], [3,4,5]), comparison='log2fc')
    """
    try:
        # Fetch data
        expression_df, metadata = fetch_geo_series_matrix(geo_id)
        
        # Handle sample selection/comparison
        if sample_group is not None:
            if isinstance(sample_group, tuple) and len(sample_group) == 2:
                # Comparison mode
                comparison = comparison or 'log2fc'
                expression_df = compare_samples(expression_df, sample_group, comparison)
                expression_df = expression_df.to_frame('comparison')
            elif isinstance(sample_group, (list, int)):
                # Subset mode
                if isinstance(sample_group, int):
                    sample_group = [sample_group]
                expression_df = expression_df.iloc[:, sample_group]
        
        # Extract genes
        gene_expression = extract_gene_expression(expression_df, gene_list)
        
        return gene_expression.to_dict()
        
    except Exception as e:
        print(f"Error fetching GEO data: {e}")
        print("Returning default expression values (all 1.0)")
        return {gene: 1.0 for gene in gene_list}


def list_geo_samples(geo_id):
    """
    List available samples in a GEO dataset.
    
    Parameters
    ----------
    geo_id : str
        GEO dataset ID
        
    Returns
    -------
    pd.DataFrame
        Sample metadata
    """
    try:
        _, metadata = fetch_geo_series_matrix(geo_id)
        return metadata
    except:
        print(f"Could not fetch metadata for {geo_id}")
        return None
