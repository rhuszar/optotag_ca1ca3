function [ca1_tag, ca3_tag] = get_taggedCells( cell_metrics, varargin )

p=inputParser;
addParameter(p,'celltype','pyr',@ischar);
addParameter(p,'same_probe',true,@islogical);
addParameter(p,'remove_longResponse',true,@islogical);

parse(p,varargin{:});

celltype = p.Results.celltype;
same_probe = p.Results.same_probe;
remove_longResponse = p.Results.remove_longResponse;

% cell types 
if strcmp(celltype, 'pyr')
    cell_id = cellfun(@(x) any(regexp(x, 'Pyr')), cell_metrics.putativeCellType )';
elseif strcmp(celltype, 'int')
    cell_id = cellfun(@(x) any(regexp(x, 'Int')), cell_metrics.putativeCellType )';
else
    cell_id = true( cell_metrics.general.cellCount, 1 );
end

cell_region = cell_metrics.brainRegion;
cell_shankID = arrayfun(@(x) [ 'S' num2str( cell_metrics.shankID(x) ) ], 1:cell_metrics.general.cellCount , 'UniformOutput', false)';
opto_id = cell_metrics.optoTag.opto_id;

uniqLeds = cell_metrics.optoTag.lightSource(1,:);

ca1_tag = false( cell_metrics.general.cellCount, 1 );
ca3_tag = false( cell_metrics.general.cellCount, 1 );
for jp = 1:cell_metrics.general.cellCount
    
    % get led sources for this cell- specified by region
    if same_probe
        thiscell_ledsources = cellfun(@(x) any( regexp( x, cell_shankID{jp})), uniqLeds );
    else
        thiscell_ledsources = ~cellfun(@(x) any(regexp(x, 'ignor')), uniqLeds );
    end

    if any( opto_id(jp, thiscell_ledsources) )
        % CA1 light responsive
        if strcmp( cell_region{jp}, 'CA1' )
            ca1_tag( jp ) = true;
         % CA3 light responsive
        elseif strcmp( cell_region{jp}, 'CA3' )       
            ca3_tag( jp ) = true;
        end
    end
end

ca1_tag = ca1_tag & cell_id;
if remove_longResponse && isfield(cell_metrics.optoTag, 'long_response' )
    ca1_tag( cell_metrics.optoTag.long_response ) = false;
end
ca3_tag = ca3_tag & cell_id;


