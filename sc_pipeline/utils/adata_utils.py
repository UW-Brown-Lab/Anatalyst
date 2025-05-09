def save_layer(adata, name, data=None, make_active=False):
    """
    Save data as a layer in AnnData with clear naming.
    
    Args:
        adata: AnnData object
        name: Name for the layer
        data: Data to store (if None, uses current X)
        make_active: Whether to make this the active layer (X)
        
    Returns:
        adata: The updated AnnData object
    """
    # Store data as a layer
    if data is None:
        adata.layers[name] = adata.X.copy()
    else:
        adata.layers[name] = data
    
    # Make active if requested
    if make_active:
        adata.X = adata.layers[name].copy()
        
        # Track the active layer in uns
        if 'active_layer' not in adata.uns:
            adata.uns['active_layer'] = {}
        
        adata.uns['active_layer']['current'] = name
    
    return adata

def set_active_layer(adata, layer_name):
    """
    Set a specific layer as the active layer (X) in AnnData.
    
    Args:
        adata: AnnData object
        layer_name: Name of the layer to make active
        
    Returns:
        adata: The updated AnnData object
    """
    if layer_name not in adata.layers:
        raise ValueError(f"Layer '{layer_name}' not found in AnnData object")
    
    # Set the specified layer as X
    adata.X = adata.layers[layer_name].copy()
    
    # Track the active layer in uns
    if 'active_layer' not in adata.uns:
        adata.uns['active_layer'] = {}
    
    adata.uns['active_layer']['current'] = layer_name
    
    return adata