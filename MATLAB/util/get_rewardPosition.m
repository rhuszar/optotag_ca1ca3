function rewardCoordinates = get_rewardPosition( basepath, varargin )
%{

return reward coordinates
either maze locations- 'maze'- or locations where animal drank- 'animal'

%}
    
    p=inputParser;
    addParameter(p,'coordinates','maze',@ischar);     % 'maze' or 'animal'      
    parse(p,varargin{:});
    coordinates = p.Results.coordinates;

    basename = bz_BasenameFromBasepath( basepath );
    
    % maze coordinates- recorded based on video
    if strcmp( 'maze', coordinates )
        load( fullfile( basepath, 'rewardCoordinates.mat' ) )
    else
        load( fullfile( basepath, [ basename '.Tracking.Behavior.mat' ] ) )
        rewT = cell2mat( tracking.events.rewardTime );
        rewO = cell2mat( tracking.events.rewardOrder );
        rewPos = [ interp1( tracking.timestamps, tracking.position.x, rewT ) ...
                    interp1( tracking.timestamps, tracking.position.y, rewT) ];
        rewardCoordinates = [ mean( rewPos(rewO == 1,:) ); mean( rewPos(rewO == 2,:) ) ];
    end
    
    
end