function indicators = getCellIndicators_cheese( basepath, varargin )
% load cell indicators 

p=inputParser;
addParameter(p,'savemat',true, @islogical );
addParameter(p,'forceReload',false, @islogical );
addParameter(p,'extract_tagging',true, @islogical );
parse(p,varargin{:});

savemat = p.Results.savemat;
forceReload = p.Results.forceReload;
extract_tagging = p.Results.extract_tagging;

basename = bz_BasenameFromBasepath( basepath );
filepath = fullfile( basepath, [basename '.cellIndicators.mat'] );

% check if file exists
if exist( filepath ) & ~forceReload
    load( filepath )
    return
end

disp('Loading metrics and computing indicators...')
load( fullfile( basepath, [basename '.cell_metrics.cellinfo.mat'] ) )

unknown_cell = find( cellfun(@isempty, cell_metrics.brainRegion ) );
if ~isempty( unknown_cell )
    cell_metrics.brainRegion( unknown_cell ) = {''};
end

% output struct 
indicators = struct();

% cell types
indicators.pyr_id = cellfun(@(x) any( regexp( x, 'Pyr' ) ), cell_metrics.putativeCellType )';
indicators.int_id = cellfun(@(x) any( regexp( x, 'Int' ) ), cell_metrics.putativeCellType )'; 

% brain regions
indicators.ca1_id = cellfun(@(x) strcmp( x, 'CA1'), cell_metrics.brainRegion )';
indicators.ca3_id = cellfun(@(x) strcmp( x, 'CA3'), cell_metrics.brainRegion )';

% hpc pyramidal cells
indicators.ca1_pyr = indicators.pyr_id & indicators.ca1_id;
indicators.ca3_pyr = indicators.pyr_id & indicators.ca3_id;
indicators.ca1_int = indicators.int_id & indicators.ca1_id;
indicators.ca3_int = indicators.int_id & indicators.ca3_id;

indicators.hpc_pyr = indicators.ca1_pyr | indicators.ca3_pyr;
indicators.hpc_int = indicators.ca1_int | indicators.ca3_int;

if isfield( 'probeID', cell_metrics )
    probe_id = cellfun(@(x) strcmp(x, 'P1'), cell_metrics.probeID ) + 2.*cellfun(@(x) strcmp(x, 'P2'), cell_metrics.probeID ); 
    indicators.probe_id = probe_id(:);
end
indicators.maxWaveformCh1 = cell_metrics.maxWaveformCh1( : );

if isfield( cell_metrics, 'optoTag' ) && extract_tagging
    % get light responsive neurons
    [indicators.ca1_tag_sameShank , indicators.ca3_tag_sameShank] = get_taggedCells(cell_metrics, 'same_probe', true);
    [indicators.ca1_tag_all, indicators.ca3_tag_all] = get_taggedCells(cell_metrics, 'same_probe', false);
    % neurons that are exclusively driven by stimulating a different shank
    indicators.ca1_tag_diffShank = indicators.ca1_tag_all; indicators.ca1_tag_diffShank( indicators.ca1_tag_sameShank ) = false;
    indicators.ca3_tag_diffShank = indicators.ca3_tag_all; indicators.ca3_tag_diffShank( indicators.ca3_tag_sameShank ) = false;
end

% save
if savemat
    save( filepath, 'indicators' )
end
    

end