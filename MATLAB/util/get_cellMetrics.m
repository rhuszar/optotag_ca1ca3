function cell_metrics = get_cellMetrics( basepath )

    basename = bz_BasenameFromBasepath( basepath );
    load( fullfile( basepath, [basename '.cell_metrics.cellinfo.mat' ] ) )
    
end