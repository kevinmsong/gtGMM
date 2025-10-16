# Changelog

All notable changes to the gtGMM project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2024-10-16

### Major Refactoring and Enhancements

This is a major refactoring release that incorporates all optimizations and enhancements developed during the cardiac regeneration analysis project.

### Added

- **New `enrichment.py` module**: Comprehensive pathway enrichment analysis with Fisher's exact test and FDR correction
- **New `data.py` module**: Data loading utilities including:
  - `load_gse130699()` function for cardiac regeneration dataset
  - `get_cardiac_pathways()` function with 12 curated pathways
  - Expression data filtering and normalization utilities
- **Convenience functions**: High-level API functions for common workflows:
  - `create_terrain()` for quick terrain creation with optimal settings
  - `analyze_regeneration()` for complete cardiac regeneration analysis
- **Enhanced documentation**: Full Sphinx documentation with API references
- **Comprehensive test suite**: Test script covering all major functionality

### Changed

- **Default resolution**: Increased from 1000x1000 to 500x500 for optimal balance of detail and performance
- **GMM optimization**: Now uses adaptive BIC-based model selection by default
- **Package structure**: Reorganized into a cleaner, more modular structure
- **API improvements**: More intuitive method names and clearer parameter descriptions
- **Version bump**: Updated to 2.0.0 to reflect major changes

### Enhanced

- **`GeneTerrain` class**: 
  - Improved `optimize_sigma()` method with better default parameters
  - Added `get_terrain_data()` method for easier GMM integration
- **`GMMOptimizer` class**:
  - Added `extract_gene_components()` method for component-to-gene mapping
  - Enhanced optimization metrics and reporting
- **`TerrainTDA` class**:
  - Tighter integration with `GeneTerrain`
  - Improved TDA metric calculations
- **`TerrainVisualizer` class**:
  - Publication-quality default settings
  - Better colormap handling (jet colormap for terrains)

### Fixed

- Fisher's exact test now properly handles edge cases with negative values
- Improved error handling throughout the package
- Better handling of empty interaction sets from STRING database

### Optimizations

- Terrain creation now uses optimized sigma values by default
- GMM analysis uses adaptive component selection
- TDA calculations are more efficient for high-resolution terrains

### Documentation

- Complete Sphinx documentation site
- Updated README with v2.0 features
- API reference for all modules
- Comprehensive docstrings throughout

## [1.0.0] - 2024-01-01

### Initial Release

- Basic GeneTerrain functionality
- GMM decomposition
- TDA analysis
- Visualization tools

