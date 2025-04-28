# sc_pipeline/modules/filtering.py

#### WIP NEEDS SOME TLC

import logging
from unicodedata import numeric
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation
from sc_pipeline.core.data_context import DataContext
from sc_pipeline.core.module import AnalysisModule

class FilterManager:
    def __init__(self, logger=None):
        self.logger = logger or logging.getLogger("FilterManager")
        self.filter_handlers = {
            'numeric': self._apply_numeric_filter,
            'boolean': self._apply_boolean_filter,
            'categorical': self._apply_categorical_filter,
            # Add more filter types as needed
        }
        self.filter_stats = {
            'cells_removed_by_filter': {},
            'filter_details': {}
        }
    
    def apply_filters(self, adata:AnnData, filters:dict, inplace=True, n_mads=3.0):
        """Apply all filters to the AnnData object
        
        Args:
            adata: AnnData object to filter. Should contain the metrics being used for filtering.
            filters: Dictionary of filter configurations
            inplace: Whether to modify the adata inplace
            n_mads: Number of Median Absolute Deviations for outlier detection

        Returns:
            Filtered AnnData object
        """
        if not inplace:
            adata = adata.copy()

            original_cells = adata.n_obs

        metric: str
        config: dict
        for metric, config in filters.items():
            filter_type = config.get('type')
            if filter_type in self.filter_handlers:
                adata = self.filter_handlers[filter_type](adata, metric, config, n_mads)
            else:
                self.logger.warning(f"Unknown filter type '{filter_type}' for metric '{metric}'")

        # Log overall results
        remaining_cells = adata.n_obs
        removed_cells = original_cells - remaining_cells
        self.logger.info(f"Filtering completed: {removed_cells}/{original_cells} cells removed ({remaining_cells} remaining)")

        return adata
    
    def _apply_numeric_filter(self, adata, metric, config, n_mads):
        """
        Handle numeric filters with min/max thresholds and outlier detection
        
        Args:
            adata: AnnData object
            metric: Name of the metric to filter on
            config: Filter configuration
            n_mads: Number of MADs for outlier detection
            
        Returns:
            Filtered AnnData object
        """
        # Check if the metric exists -- this allows metric name to be different than actual column name if ever needed for some reason?
        column = config.get('column', metric)
        
        if column not in adata.obs.columns:
            self.logger.warning(f"Numeric filter metric '{metric}' not found in data")
            return adata
        else:
            # Use the specified column name
            metric = column
        
        # Apply minimum threshold if specified
        if 'min' in config and config['min'] is not None:
            min_val = config['min']
            mask = adata.obs[metric] >= min_val
            n_removed = np.sum(~mask)
            if n_removed > 0:
                self.logger.info(f"Removing {n_removed} cells with {metric} < {min_val}")
                self.filter_stats['filter_details'][f'min_{metric}'] = min_val
                self.filter_stats['cells_removed_by_filter'][f'min_{metric}'] = int(n_removed)
                adata = adata[mask]
        
        # Apply maximum threshold if specified
        if 'max' in config and config['max'] is not None:
            max_val = config['max']
            mask = adata.obs[metric] <= max_val
            n_removed = np.sum(~mask)
            if n_removed > 0:
                self.logger.info(f"Removing {n_removed} cells with {metric} > {max_val}")
                self.filter_stats['filter_details'][f'max_{metric}'] = max_val
                self.filter_stats['cells_removed_by_filter'][f'max_{metric}'] = int(n_removed)
                adata = adata[mask]
        
        # Apply outlier detection if requested
        if config.get('outlier_detection', True):
            # Get MAD threshold (use global or metric-specific)
            mad_threshold = config.get('n_mads', n_mads)
            
            # Determine if filtering should be one-sided or two-sided
            two_sided = config.get('two_sided', True)
            
            # Detect outliers
            outlier_mask = self._detect_outliers(
                adata.obs[metric],
                mad_threshold,
                filter_high=True,
                filter_low=two_sided
            )

            # Create a column in adata.obs to mark outliers
            outlier_col = f"{metric}_outlier"
            adata.obs[outlier_col] = outlier_mask

            remove_outliers = config.get('remove_outliers', False)
            
            n_outliers = np.sum(outlier_mask)
            if n_outliers > 0:
                self.logger.info(f"Marked {n_outliers} cells as outliers in {metric} using {mad_threshold} MADs")
                self.filter_stats['filter_details'][f'outlier_{metric}'] = mad_threshold
                self.filter_stats['cells_removed_by_filter'][f'outlier_{metric}'] = int(n_outliers) if remove_outliers else 0

                if remove_outliers:
                    self.logger.info(f"Removing {n_outliers} cells marked as outliers in {metric}")
                    adata = adata[~outlier_mask]
        
        return adata
    
    def _apply_boolean_filter(self, adata, metric, config, n_mads):
        """
        Handle boolean filters
        
        Args:
            adata: AnnData object
            metric: Name of the metric or description
            config: Filter configuration
            n_mads: Not used for boolean filters, just here for flexible method calling
            
        Returns:
            Filtered AnnData object
        """
        # Get the column name (might be different from metric name)
        column = config.get('column', metric)
        
        # Check if the column exists
        if column not in adata.obs.columns:
            self.logger.warning(f"Boolean filter column '{column}' not found in data")
            return adata
        else: 
            metric = column
        
        # Get the expected value
        desired_value = config.get('keep', True)
        
        # Convert to boolean if needed
        if adata.obs[column].dtype != bool:
            try:
                # Create a new column with the boolean conversion
                bool_column = f"{column}_bool"
                
                if pd.api.types.is_numeric_dtype(adata.obs[column]):
                    # For numeric columns, treat >0 as True
                    adata.obs[bool_column] = adata.obs[column] > 0
                    self.logger.info(f"Converting numeric column '{column}' to boolean '{bool_column}'")
                elif pd.api.types.is_string_dtype(adata.obs[column]):
                    # For string columns, convert common truth values
                    adata.obs[bool_column] = adata.obs[column].str.lower().isin(['true', 'yes', '1', 'y'])
                    self.logger.info(f"Converting string column '{column}' to boolean '{bool_column}'")
                else:
                    self.logger.warning(f"Could not convert column '{column}' to boolean, skipping filter")
                    return adata
                
                # Use the new boolean column
                column = bool_column
                
            except Exception as e:
                self.logger.warning(f"Could not convert column '{column}' to boolean: {e}")
                return adata
        
        # Apply the filter
        mask = (adata.obs[column] == desired_value)
        n_removed = np.sum(~mask)
        
        if n_removed > 0:
            self.logger.info(f"Removing {n_removed} cells with {column} != {desired_value}")
            self.filter_stats['filter_details'][f'bool_{column}'] = desired_value
            self.filter_stats['cells_removed_by_filter'][f'bool_{column}'] = int(n_removed)
            adata = adata[mask]
        
        return adata
    
    def _apply_categorical_filter(self, adata, metric, config, n_mads):
        """
        Handle categorical filters with keep/drop lists
        
        Args:
            adata: AnnData object
            metric: Name of the metric or description
            config: Filter configuration
            n_mads: Not used for categorical filters
            
        Returns:
            Filtered AnnData object
        """
        # Get the column name (might be different from metric name)
        column = config.get('column', metric)
        
        # Check if the column exists
        if column not in adata.obs.columns:
            self.logger.warning(f"Categorical filter column '{column}' not found in data")
            return adata
        
        # Check for categories to keep
        if 'keep' in config:
            categories_to_keep = config['keep']
            if not isinstance(categories_to_keep, list):
                categories_to_keep = [categories_to_keep]
            
            mask = adata.obs[column].isin(categories_to_keep)
            n_removed = np.sum(~mask)
            
            if n_removed > 0:
                self.logger.info(f"Keeping only {len(mask) - n_removed} cells with {column} in {categories_to_keep}")
                self.filter_stats['filter_details'][f'keep_{column}'] = categories_to_keep
                self.filter_stats['cells_removed_by_filter'][f'keep_{column}'] = int(n_removed)
                adata = adata[mask]
        
        # Check for categories to drop
        if 'drop' in config:
            categories_to_drop = config['drop']
            if not isinstance(categories_to_drop, list):
                categories_to_drop = [categories_to_drop]
            
            mask = ~adata.obs[column].isin(categories_to_drop)
            n_removed = np.sum(~mask)
            
            if n_removed > 0:
                self.logger.info(f"Removing {n_removed} cells with {column} in {categories_to_drop}")
                self.filter_stats['filter_details'][f'drop_{column}'] = categories_to_drop
                self.filter_stats['cells_removed_by_filter'][f'drop_{column}'] = int(n_removed)
                adata = adata[mask]
        
        return adata

    def _detect_outliers(self, values, nmads, filter_high=True, filter_low=True):
        """
        Detect outliers in a Series using Median Absolute Deviation
        
        Args:
            values: Pandas Series of values
            nmads: Number of MADs to use as threshold
            filter_high: Whether to detect high outliers
            filter_low: Whether to detect low outliers
            
        Returns:
            Boolean mask of outliers
        """
        median = np.median(values)
        mad = median_abs_deviation(values, nan_policy='omit')
        
        # Initialize mask with no outliers
        outlier = np.zeros(len(values), dtype=bool)
        
        # Add low outliers if requested
        if filter_low:
            outlier = outlier | (values < median - nmads * mad)
            
        # Add high outliers if requested
        if filter_high:
            outlier = outlier | (values > median + nmads * mad)
            
        return outlier

class Filtering(AnalysisModule):
    """
    Module for filtering cells based on QC metrics.
    Supports both hard cutoffs and automatic outlier detection.
    """

    PARAMETER_SCHEMA = {
        'filters': {
            'type': dict,
            'default': {
                'n_genes_by_counts':{
                    'type': 'numeric',
                    'min': None,
                    'max': None,
                    'outlier_detection': True,
                    'remove_outliers': True,
                    'two_sided': True
                },
                'pct_counts_mt':{
                    'type': 'numeric',
                    'min': None,
                    'max': None,
                    'outlier_detection': True,
                    'remove_outliers': True,
                    'two_sided': False
                },
                'total_counts':{
                    'type': 'numeric',
                    'min': None,
                    'max': None,
                    'outlier_detection': True,
                    'remove_outliers': True,
                    'two_sided': True
                }
            },
            'description': 'Dictionary of filters to apply, each labeled according to their type and method.'
        },
        'outlier_detection': {
            'type': bool,
            'default': True,
            'description': 'Whether to use automatic outlier detection using Median Absolute Deviation method.'
        },
        'n_mads': {
            'type': float,
            'default': 3.0,
            'description': 'Number of MADs for outlier detection. Common threshold is 3.0'
        },
        'create_plots': {
            'type': bool,
            'default': True,
            'description': 'Generate plots of filtered cells'
        },
        'inplace': {
            'type': bool,
            'default': True,
            'description': 'Whether to filter the AnnData object in place (False not supported yet)'
        },
    }
    
    def __init__(self, name, params):
        super().__init__(name, params)
        self.logger = logging.getLogger(f"Module.{name}")
        
        # Define required inputs and outputs
        self.required_inputs = ["data"]
        self.outputs = ["data"]
        
    def is_outlier(self, adata, metric, nmads, high=True, low=True):
        """
        Detect outliers in a given metric using Median Absolute Deviation.
        Inspired by and borrowed from https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

        Args:
            adata: AnnData object
            metric: Metric to check for outliers
            nmads: Number of MADs to use as threshold
            high: Whether to detect high outliers
            low: Whether to detect low outliers
            
        Returns:
            Boolean mask of outliers
        """
        M = adata.obs[metric]
        median = np.median(M)
        mad = median_abs_deviation(M, nan_policy='omit')
        
        # Initialize mask with no outliers
        outlier = np.zeros(len(M), dtype=bool)
        
        # Add low outliers if requested
        if low:
            outlier = outlier | (M < median - nmads * mad)
            
        # Add high outliers if requested
        if high:
            outlier = outlier | (M > median + nmads * mad)
            
        return outlier
        
    def run(self, data_context: DataContext):
        """
        Filter cells based on QC metrics.
        
        Args:
            data_context: The shared data context
            
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            # Get the data
            adata = data_context.get("data")
            
            # Store original counts
            original_cells = adata.n_obs
            original_genes = adata.n_vars
            
            # Initialize the filter manager
            filter_manager = FilterManager(logger=self.logger)
            
            # Apply all configured filters
            adata = filter_manager.apply_filters(
                adata,
                self.params.get('filters', {}),
                inplace=self.params.get('inplace', True),
                n_mads=self.params.get('n_mads', 3.0)
            )
            
            # Get filtering statistics
            filter_stats = filter_manager.get_filter_stats()
            
            # Store stats in the AnnData object
            adata.uns['filtering'] = {
                'original_cells': int(original_cells),
                'remaining_cells': int(adata.n_obs),
                'percent_cells_kept': float(adata.n_obs / original_cells * 100) if original_cells > 0 else 0,
                'original_genes': int(original_genes),
                'remaining_genes': int(adata.n_vars),
                'filter_details': filter_stats['filter_details'],
                'cells_removed_by_filter': filter_stats['cells_removed_by_filter']
            }
            
            # Create plots if requested
            if self.params.get('create_plots', True):
                self._create_plots(adata, data_context)
            
            # Update the data context
            data_context.set("data", adata)
            
            return True
            
        except Exception as e:
            self.logger.error(f"Error filtering cells: {e}", exc_info=True)
            return False
            
    def _create_plots(self, adata, metrics, data_context):
        """Create plots showing the filtering results."""
        # Create a violin plot of the metrics after filtering
        if len(metrics) > 0:
            available_metrics = [m for m in metrics if m in adata.obs.columns]
            if available_metrics:
                fig, axs = plt.subplots(1, len(available_metrics), figsize=(4*len(available_metrics), 5))
                if len(available_metrics) == 1:
                    axs = [axs]
                
                for i, metric in enumerate(available_metrics):
                    sc.pl.violin(adata, metric, ax=axs[i], show=False)
                    axs[i].set_title(f"{metric} after filtering")
                
                plt.tight_layout()
                
                # Save the figure
                img_path = self.save_figure(self.name, fig, name="qc_after_filtering")
                
                # Add to report
                data_context.add_figure(
                    module_name=self.name,
                    title="QC Metrics After Filtering",
                    description="Distribution of quality control metrics after cell filtering",
                    image_path=img_path,
                    caption=f"Filtered from {adata.uns['filtering']['original_cells']} to {adata.uns['filtering']['remaining_cells']} cells ({adata.uns['filtering']['percent_cells_kept']:.1f}% kept)"
                )
                
                # Create a bar chart of removed cells by filter
                cells_removed = adata.uns['filtering'].get('cells_removed_by_filter', {})
                if cells_removed:
                    fig, ax = plt.subplots(figsize=(8, 4))
                    filters = list(cells_removed.keys())
                    counts = list(cells_removed.values())
                    
                    # Sort by number of cells removed
                    idx = np.argsort(counts)[::-1]
                    filters = [filters[i] for i in idx]
                    counts = [counts[i] for i in idx]
                    
                    # Create bar chart
                    ax.bar(filters, counts)
                    ax.set_ylabel('Cells Removed')
                    ax.set_title('Cells Removed by Filter Type')
                    plt.xticks(rotation=45, ha='right')
                    plt.tight_layout()
                    
                    # Save the figure
                    img_path = self.save_figure(self.name, fig, name="cells_removed_by_filter")
                    
                    # Add to report
                    data_context.add_figure(
                        module_name=self.name,
                        title="Cells Removed by Filter Type",
                        description="Number of cells removed by each filter criterion",
                        image_path=img_path
                    )