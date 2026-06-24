function out = interpToBouts( basepath, ts_target, varargin )
% interpToBouts( basepath, ts_target, varargin )
%     interpolate cheeseboard bout data to target timestamps  

p=inputParser;
addParameter(p,'subsess_n',2,@isnumeric);

parse(p,varargin{:});
subsess_n = p.Results.subsess_n;

load( fullfile( basepath, 'run_bouts_5cms.mat' ) )

% running bouts during behavior
trial_bout_ints = int_cat( subsess_cat == subsess_n,: );
trial_bout_rewID = rewID_cat( subsess_cat == subsess_n );
trial_bout_rewB = rewBout_cat( subsess_cat == subsess_n );
trial_n = epoch_cat( subsess_cat == subsess_n );

trial_v_all = [];
trial_ts_all = [];
trial_pos_all = [];
trial_rewID_all = [];
trial_rewBout_all = [];
trial_n_all = [];
bout_id_all = [];
% for each 
for jp = 1:size( trial_bout_ints, 1 )
    
    ipp = InIntervals( ts_all, trial_bout_ints( jp,: ) );
    trial_v_all = [ trial_v_all ; v_all(ipp) ];
    trial_ts_all = [ trial_ts_all ; ts_all(ipp) ];
    trial_pos_all = [trial_pos_all ; pos_all(ipp,:) ];
    trial_rewID_all = [ trial_rewID_all ; repmat( trial_bout_rewID(jp), sum(ipp), 1 )];
    trial_rewBout_all = [ trial_rewBout_all ; repmat( trial_bout_rewB(jp), sum(ipp), 1 )];
    trial_n_all = [ trial_n_all ; repmat( trial_n(jp), sum(ipp), 1 ) ];
    bout_id_all = add_to_mat( bout_id_all, repmat( jp, sum(ipp), 1 ) );
end

% interpolate bout-specific data to target timestamps
out = struct();
xpos = interp1( trial_ts_all, trial_pos_all(:,1), ts_target );
ypos = interp1( trial_ts_all, trial_pos_all(:,2), ts_target );
out.pos = [ xpos ypos ];
out.v = interp1( trial_ts_all, trial_v_all, ts_target );
out.rewID = round( interp1( trial_ts_all, trial_rewID_all, ts_target ) );
out.rewBout = round( interp1( trial_ts_all, trial_rewBout_all, ts_target ) );
out.trialn = round( interp1( trial_ts_all, trial_n_all, ts_target ) );
out.boutid = round( interp1( trial_ts_all, bout_id_all, ts_target ) );
