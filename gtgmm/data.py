"""
Data Loading and Preprocessing Module for gtGMM

This module provides utilities for loading gene expression data,
including the GSE130699 cardiac regeneration dataset.
"""

import numpy as np
import pandas as pd
from typing import Tuple, Dict, List, Optional
import warnings

warnings.filterwarnings('ignore')


def load_gse130699(data_path: Optional[str] = None) -> Tuple[List[str], Dict[str, float], Dict[str, float]]:
    """
    Load the GSE130699 cardiac regeneration dataset.
    
    This dataset contains single-cell RNA-seq data from neonatal mouse hearts
    at postnatal day 1 (P1, regenerative) and day 8 (P8, non-regenerative).
    
    Parameters
    ----------
    data_path : str, optional
        Path to the processed data file. If None, uses default location.
        
    Returns
    -------
    genes : list
        List of gene symbols
    p1_expression : dict
        P1 expression values {gene: expression}
    p8_expression : dict
        P8 expression values {gene: expression}
        
    Notes
    -----
    The data is preprocessed to include only genes with sufficient expression
    and variability across both conditions.
    """
    try:
        from .gse130699_data import GSE130699Data
        
        # Load data
        data_loader = GSE130699Data()
        genes, p1_expr, p8_expr = data_loader.load_data()
        
        return genes, p1_expr, p8_expr
        
    except ImportError:
        raise ImportError(
            "GSE130699 data module not found. Please ensure gse130699_data.py is available."
        )


def load_expression_matrix(file_path: str, 
                          gene_column: str = 'gene',
                          expression_columns: Optional[List[str]] = None) -> pd.DataFrame:
    """
    Load gene expression data from a CSV or TSV file.
    
    Parameters
    ----------
    file_path : str
        Path to the expression data file
    gene_column : str
        Name of the column containing gene symbols
    expression_columns : list, optional
        List of column names to use as expression values. If None, uses all numeric columns.
        
    Returns
    -------
    pd.DataFrame
        Expression matrix with genes as rows
    """
    # Detect file format
    if file_path.endswith('.csv'):
        df = pd.read_csv(file_path)
    elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
        df = pd.read_csv(file_path, sep='\t')
    else:
        raise ValueError("File must be .csv or .tsv format")
    
    # Set gene column as index
    if gene_column in df.columns:
        df = df.set_index(gene_column)
    
    # Select expression columns
    if expression_columns is not None:
        df = df[expression_columns]
    else:
        # Use all numeric columns
        df = df.select_dtypes(include=[np.number])
    
    return df


def aggregate_expression(expression_matrix: pd.DataFrame, 
                        method: str = 'mean') -> Dict[str, float]:
    """
    Aggregate expression values across multiple samples/cells.
    
    Parameters
    ----------
    expression_matrix : pd.DataFrame
        Expression matrix with genes as rows and samples as columns
    method : str
        Aggregation method: 'mean', 'median', or 'sum'
        
    Returns
    -------
    dict
        Aggregated expression values {gene: expression}
    """
    if method == 'mean':
        aggregated = expression_matrix.mean(axis=1)
    elif method == 'median':
        aggregated = expression_matrix.median(axis=1)
    elif method == 'sum':
        aggregated = expression_matrix.sum(axis=1)
    else:
        raise ValueError(f"Unknown aggregation method: {method}")
    
    return aggregated.to_dict()


def filter_genes(expression_dict: Dict[str, float],
                min_expression: float = 0.0,
                max_genes: Optional[int] = None) -> Dict[str, float]:
    """
    Filter genes based on expression criteria.
    
    Parameters
    ----------
    expression_dict : dict
        Gene expression values {gene: expression}
    min_expression : float
        Minimum expression threshold
    max_genes : int, optional
        Maximum number of genes to keep (selects top expressed)
        
    Returns
    -------
    dict
        Filtered expression values
    """
    # Filter by minimum expression
    filtered = {gene: expr for gene, expr in expression_dict.items() 
                if expr >= min_expression}
    
    # Select top genes if specified
    if max_genes is not None and len(filtered) > max_genes:
        sorted_genes = sorted(filtered.items(), key=lambda x: x[1], reverse=True)
        filtered = dict(sorted_genes[:max_genes])
    
    return filtered


def normalize_expression(expression_dict: Dict[str, float],
                        method: str = 'zscore') -> Dict[str, float]:
    """
    Normalize gene expression values.
    
    Parameters
    ----------
    expression_dict : dict
        Gene expression values {gene: expression}
    method : str
        Normalization method: 'zscore', 'minmax', or 'log2'
        
    Returns
    -------
    dict
        Normalized expression values
    """
    values = np.array(list(expression_dict.values()))
    genes = list(expression_dict.keys())
    
    if method == 'zscore':
        mean = values.mean()
        std = values.std()
        if std > 0:
            normalized = (values - mean) / std
        else:
            normalized = values - mean
    
    elif method == 'minmax':
        vmin = values.min()
        vmax = values.max()
        if vmax > vmin:
            normalized = (values - vmin) / (vmax - vmin)
        else:
            normalized = np.zeros_like(values)
    
    elif method == 'log2':
        normalized = np.log2(values + 1)
    
    else:
        raise ValueError(f"Unknown normalization method: {method}")
    
    return dict(zip(genes, normalized))


def get_cardiac_pathways() -> Dict[str, List[str]]:
    """
    Get curated cardiac regeneration pathway gene sets.
    
    Returns
    -------
    dict
        Dictionary mapping pathway names to gene lists
        
    Notes
    -----
    These pathways are curated for cardiac regeneration analysis and include:
    - Proliferation pathways (Hippo/YAP, AKT/PI3K, Cell Cycle)
    - Developmental signaling (WNT, Notch, FGF/VEGF)
    - Stress response (TGF-β/BMP, P38 MAPK)
    - Metabolism (Mitochondrial)
    - Epigenetic regulation
    - Cardiac structural genes
    """
    pathways = {
        'Hippo/YAP': ['Yap1', 'Wwtr1', 'Tead1', 'Tead2', 'Tead3', 'Tead4', 
                      'Lats1', 'Lats2', 'Stk3', 'Stk4', 'Sav1', 'Mob1a', 'Mob1b'],
        
        'AKT/PI3K': ['Akt1', 'Akt2', 'Akt3', 'Pik3ca', 'Pik3cb', 'Pik3r1', 
                     'Pik3r2', 'Pten', 'Pdk1', 'Pdk2', 'Mtor', 'Rptor', 'Rictor'],
        
        'Cell Cycle': ['Ccnd1', 'Ccnd2', 'Ccnd3', 'Cdk1', 'Cdk2', 'Cdk4', 
                       'Cdk6', 'Cdkn1a', 'Cdkn1b', 'Cdkn2a', 'Cdkn2b', 'Aurka', 'Aurkb'],
        
        'WNT': ['Wnt1', 'Wnt3a', 'Wnt5a', 'Wnt7a', 'Wnt11', 'Ctnnb1', 
                'Axin1', 'Axin2', 'Gsk3b', 'Lef1', 'Tcf7', 'Dkk1', 'Sfrp1'],
        
        'Notch': ['Notch1', 'Notch2', 'Dll1', 'Dll4', 'Jag1', 'Jag2', 
                  'Rbpj', 'Hes1', 'Hes5', 'Hey1', 'Hey2', 'Numb'],
        
        'FGF/VEGF': ['Fgf1', 'Fgf2', 'Fgf10', 'Fgfr1', 'Fgfr2', 'Vegfa', 
                     'Vegfb', 'Vegfc', 'Vegfr1', 'Vegfr2', 'Kdr', 'Flt1'],
        
        'TGF-β/BMP': ['Tgfb1', 'Tgfb2', 'Tgfb3', 'Tgfbr1', 'Tgfbr2', 
                      'Smad2', 'Smad3', 'Smad4', 'Smad7', 'Bmp2', 'Bmp4', 'Bmp7', 'Bmpr1a'],
        
        'P38 MAPK': ['Mapk14', 'Map2k3', 'Map2k6', 'Map3k5', 'Mapkapk2', 
                     'Atf2', 'Jnk1', 'Jnk2'],
        
        'Mitochondrial': ['Tfam', 'Nrf1', 'Nrf2', 'Ppargc1a', 'Ppargc1b', 
                         'Cox4i1', 'Cox5a', 'Atp5a1', 'Atp5b'],
        
        'Epigenetic': ['Ezh2', 'Dnmt1', 'Dnmt3a', 'Dnmt3b', 'Tet1', 
                       'Tet2', 'Kdm6a', 'Kdm6b', 'Dip2a'],
        
        'Cardiac Structural': ['Tnnt2', 'Tnni3', 'Myh6', 'Myh7', 'Actn2', 
                              'Ttn', 'Gata4', 'Nkx2-5', 'Mef2c', 'Tbx5', 'Hand1', 'Hand2'],
        
        'FST/Follistatin': ['Fst', 'Fstl1', 'Fstl3', 'Acvr1', 'Acvr2a', 
                           'Acvr2b', 'Inhba', 'Inhbb']
    }
    
    return pathways


def validate_gene_names(genes: List[str], species: str = 'mouse') -> Tuple[List[str], List[str]]:
    """
    Validate gene names and identify invalid entries.
    
    Parameters
    ----------
    genes : list
        List of gene symbols to validate
    species : str
        Species: 'mouse' or 'human'
        
    Returns
    -------
    valid_genes : list
        List of valid gene symbols
    invalid_genes : list
        List of invalid gene symbols
        
    Notes
    -----
    This is a basic validation. For comprehensive validation,
    consider using external databases like NCBI Gene or Ensembl.
    """
    # Basic validation: check for empty strings and special characters
    valid_genes = []
    invalid_genes = []
    
    for gene in genes:
        if isinstance(gene, str) and len(gene) > 0 and gene.replace('-', '').replace('_', '').isalnum():
            valid_genes.append(gene)
        else:
            invalid_genes.append(gene)
    
    return valid_genes, invalid_genes

