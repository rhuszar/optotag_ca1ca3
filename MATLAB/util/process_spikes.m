function [X, ts] = process_spikes( spikes, varargin )
%
% Process spikes- binning / smoothing
% 
% inputs:
%     spikes - either a matrix or cell array of spike times
% 
%     varargin:
%         intervals - for binning (events mode), or specifying relevant time windows (window mode)
%         interval_mode :
%              'events' - intervals are treated as 'event' bins
%              'window' - intervals are treated as relevant windows (e.g.,
%                         behavior)- bin, smooth and correlate
%         nanpad ('window' mode) - number of padding samples
%         split_flag ('window' mode) - indicator whether to split the outputs
%         binsize  - binsize, by default the same as fwhm, in which case the
%                    spikes are just binned
%         fwhm     - smoothing window
%         unit     - 'count', 'rate', or 'z'
%         inclusion - struct
%                val   - values that are the basis for inclusion
%                ts     - timestamps
%                thr    - positive threshold on values
% 
% outputs:
%     
%     X  - binned data
%     ts - intervals for binning ('events' mode) or midpoints of timebins ('window' mode) 
%        
p=inputParser;
addParameter(p,'interval_mode', 'events', @(x) strcmp(x, 'events') || strcmp(x, 'window')  )
addParameter(p,'intervals',[],@isnumeric);  
addParameter(p,'fwhm',0.025,@isnumeric);
addParameter(p,'binsize',0.025,@isnumeric);
addParameter(p,'inclusion',struct(), @isstruct )
addParameter(p,'nanpad',5000,@isnumeric);
addParameter(p,'unit','count', @(x) strcmp(x, 'count') || strcmp(x, 'rate') || strcmp(x, 'z')  )
addParameter(p,'split_flag',false,@islogical);


parse(p,varargin{:});

interval_mode = p.Results.interval_mode;
intervals = p.Results.intervals;
inclusion = p.Results.inclusion;
fwhm = p.Results.fwhm;
binsize = p.Results.binsize;
nanpad = p.Results.nanpad;
split_flag = p.Results.split_flag;
unit = p.Results.unit;

%% check inputs
  
if ~iscell( spikes )
    error('Supplied spikes need to be a cell array .. ')
end

N = length( spikes );

% no intervals- consider all data, run in window mode
if isempty( intervals )
    intervals = [ 0 ceil( max( cellfun(@max, spikes ) ) ) ];
    interval_mode = 'window';
end


%%

    % intervals are events- e.g., ripples
    if strcmp( interval_mode, 'events' )

        X = get_spkcountMat( spikes, intervals );
        % convert to rate
        if strcmp( unit, 'rate' )
           X = X ./ diff( intervals, [], 2);
        end
        ts = intervals;

    % intervals are time windows- e.g., times on a maze 
    elseif strcmp( interval_mode, 'window' )

        % Get the behavior data binned at dt
        X = [];
        ts = [];
        for n = 1:size( intervals, 1 )
            X_tmp = bz_SpktToSpkmat( spikes, 'dt',  binsize, 'win', intervals(n,:) );
            X = [ X ;  nan( nanpad, N ) ; double( X_tmp.data ) ];  
            ts = [ ts ; nan( nanpad,1 ) ; double( X_tmp.timestamps ) ];
        end

        % Smoothing if mismatch between fwhm and binning
        if length( fwhm ) ~= length(binsize) || fwhm ~= binsize
            X_all = cell( length( fwhm ), 1 );
            for ip = 1:length( fwhm )

                stddev = fwhm(ip) ./ (2.*sqrt(2*log(2)));
                kernelx = [-fliplr(binsize:binsize:10*stddev) 0 binsize:binsize:10*stddev];
                kernel = Gauss(kernelx,0,stddev);
                kernel = kernel./Gauss(0,0,stddev);
                X_all{ip} = cell2mat( arrayfun(@(x) nanconv( X(:,x)', kernel,'same', 'nanout')', 1:N, 'UniformOutput', false ) );
                
            end
            X = X_all;
        end
        
        % ensure we're working with cell arrays 
        if ~iscell( X )
            X = {X};
        end

        % remove padded nans
        nanlogical = isnan( ts );
        for ip = 1:length(X)
           X{ip}( nanlogical,: ) = []; 
        end
        ts( nanlogical , : ) = [];
        
        % convert to rate
        if strcmp( unit, 'rate' )
           for ip = 1:length( fwhm )
               X{ip} = X{ip} ./ fwhm(ip);
           end
        elseif strcmp( unit, 'z' )
           for ip = 1:length( fwhm )
               X{ip} = zscore( X{ip} );
           end
        end
        
        
        % only include data that satisfy conditions- e.g., speed threshold,
        % theta power etc.
        if 0 ~= length( fieldnames(inclusion) )
            include = true( length( ts ), 1 );
            % go over each exclusion condition
            for c = 1:length( inclusion )
                val_interp = interp1( inclusion( c ).ts, inclusion( c ).val, ts );
                include = include & val_interp >= inclusion( c ).thr;
            end
            
            % only keep relevant data
            for ip = 1:length(X)
                X{ip} = X{ip}( include,: );
            end
            ts = ts( include,: );
        end



        if split_flag
            for ip = 1:length(X)
                % some flag if to split by interval
                X{ip} = arrayfun(@(x) X{ip}( InIntervals( ts, intervals(x,:) ), : ), 1:size(intervals,1) , 'UniformOutput', false );
            end
            ts = arrayfun(@(x) ts( InIntervals( ts, intervals(x,:) ), : ), 1:size(intervals,1) , 'UniformOutput', false );
        end
        
        % If we a
        if length(X) == 1
            X = X{1};
        end
        
    end

end