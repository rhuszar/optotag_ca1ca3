function tracking = get_tracking( basepath )
% tracking = get_tracking( basepath )
% get tracking.behavior.mat file
    
    basename = bz_BasenameFromBasepath( basepath );
    load( fullfile( basepath, [basename '.Tracking.Behavior.mat' ] ) )

end