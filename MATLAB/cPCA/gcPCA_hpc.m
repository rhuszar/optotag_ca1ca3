function gcPCA_hpc( index )
%{

- bin into ripples, then CA1 PBEs, then CA3 PBE...
- run_gcpca - ensures that variance explained is balanced across conditions  

%}

reg = 'all'; % 'CA1' or 'CA3'
speed_thresh = 5;

% number of shuffling for gcPCA
nshuff = 100;

% filter out cells that do not emit enough spikes in ripples
min_n_spk = 5;
filter_low_fr = true;

% bin sizes for various epochs
ripple_dt = 0.01;            % 10 ms ripple bins
run_dt = 1/30;             % camera dt

params = struct();
params.reg = reg;
params.speed_thresh = speed_thresh;
params.nshuff = nshuff;
params.ripple_dt = ripple_dt;
params.run_dt = run_dt;
params.min_n_spk = min_n_spk;
params.filter_low_fr = filter_low_fr;

% path to sessions folder on the HPC
datpath = '...';

index = str2double(index);
% load the session paths
load('basepaths_ca1ca3.mat')
basepath = basepaths_all{ index };

% Basepath on HPC
path_parts = strsplit(basepath,'\');
basepath_hpc = fullfile( datpath,path_parts{end-1},path_parts{end} );
basename = bz_BasenameFromBasepath( basepath_hpc );
disp(basepath_hpc)

% Define output paths
st = dbstack;
funcname = st.name;
if strcmp( reg, 'CA1' )
    outpath = fullfile(datpath, [funcname '_CA1_OUT']);
elseif strcmp( reg, 'CA3' )
    outpath = fullfile(datpath, [funcname '_CA3_OUT']);
else
    outpath = fullfile(datpath, [funcname '_OUT']);
end

if ~exist(outpath, 'dir')
    mkdir(outpath);
end
fprintf('Running %s\n', funcname)

% Output folder for this basepath
outfile = fullfile( outpath, [ basename '.mat'] );
outfile_err = fullfile( outpath, [ basename '_ERROR.mat'] );


%% Initialize the outputs struct

e = [];
details = struct();
details.description = [ 'bin data into high-synchrony events in rest periods; run contrastive PCA' ];
details.err_msg = cell(1,4);
contr_assemblies = struct();
try

    % Prepare the data
    load( fullfile( basepath_hpc, [basename '.cell_metrics.cellinfo.mat'] ) )
    load( fullfile( basepath_hpc, [basename '.Tracking.Behavior.mat'] ) )
    load( fullfile( basepath_hpc, [basename '.ripples.events.mat'] ) )
    
    indicators = getCellIndicators_cheese( basepath_hpc );
    if isfield( indicators, 'noisy_id' )
        noisy_id = indicators.noisy_id;
    else
        noisy_id = false( length( indicators.hpc_pyr ), 1 );
    end

    if strcmp( reg, 'CA1' )
        pyrid = indicators.pyr_id &  indicators.ca1_id;
    elseif strcmp( reg, 'CA3' )
        pyrid = indicators.pyr_id &  indicators.ca3_id;
    else
        pyrid = indicators.pyr_id & ( indicators.ca1_id | indicators.ca3_id );
    end
    % get rid of noisy units
    pyrid = pyrid & ~noisy_id;
    pyrid = find( pyrid );
    
    % behavior intervals
    contr_assemblies.probe1_ints = [ tracking.events.onMaze_startTime{1} tracking.events.openBox_ints{1}(:,2) ];
    contr_assemblies.learn_trial_ints = [ tracking.events.onMaze_startTime{2} cellfun( @(x) x(end), tracking.events.rewardTime ) ];
    contr_assemblies.probe2_ints = [ tracking.events.onMaze_startTime{3} tracking.events.openBox_ints{3}(:,2) ];
    
    % inclusion criteria for running process_spikes
    inclusion = struct();
    inclusion(1).val = tracking.position.v_diff_smooth;
    inclusion(1).ts = tracking.timestamps;
    inclusion(1).thr = speed_thresh;
    inclusion(2).val = double( ~isnan(tracking.position.x) & ~isnan(tracking.position.y) );
    inclusion(2).ts = tracking.timestamps;
    inclusion(2).thr = 0.9;
    
    % rest 1- take ripples and bin them finely
    rest1_ts = [ tracking.events.subSessions(1,2) tracking.events.subSessions(2,1) ];
    ripples_rest1_ts = ripples.timestamps( InIntervals( ripples.peaks, rest1_ts ), : );
    hsync_rest1_ts_binned = arrayfun(@(x) ripples_rest1_ts(x, 1) - ripple_dt/2 : ripple_dt : ripples_rest1_ts(x, 2) + ripple_dt/2, 1:size(ripples_rest1_ts, 1) , 'UniformOutput', false);
    hsync_rest1_ts_binned = cellfun(@(x) [ x(1:end-1)' x(2:end)' ], hsync_rest1_ts_binned , 'UniformOutput', false)';
    contr_assemblies.hsync_rest1_id = cell2mat( arrayfun(@(X) repmat(X, length( hsync_rest1_ts_binned{X} ), 1), 1:length( hsync_rest1_ts_binned ), 'UniformOutput', false)' );
    hsync_rest1_ts_binned = cell2mat( hsync_rest1_ts_binned );

    % rest 2- take ripples and bin them finely
    rest2_ts = [ tracking.events.subSessions(2,2) tracking.events.subSessions(3,1) ];
    ripples_rest2_ts = ripples.timestamps( InIntervals( ripples.peaks, rest2_ts ), : );
    hsync_rest2_ts_binned = arrayfun(@(x) ripples_rest2_ts(x, 1) - ripple_dt/2 : ripple_dt : ripples_rest2_ts(x, 2) + ripple_dt/2, 1:size(ripples_rest2_ts, 1) , 'UniformOutput', false);
    hsync_rest2_ts_binned = cellfun(@(x) [ x(1:end-1)' x(2:end)' ], hsync_rest2_ts_binned , 'UniformOutput', false)';
    contr_assemblies.hsync_rest2_id = cell2mat( arrayfun(@(X) repmat(X, length( hsync_rest2_ts_binned{X} ), 1), 1:length( hsync_rest2_ts_binned ), 'UniformOutput', false)' );
    hsync_rest2_ts_binned = cell2mat( hsync_rest2_ts_binned );

    % bin the data
    [ Z_probe1, contr_assemblies.probe1_ts  ] = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', contr_assemblies.probe1_ints, 'interval_mode', 'window', 'fwhm', run_dt, 'binsize', run_dt, 'inclusion', inclusion  );
    Z_rest1 = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', hsync_rest1_ts_binned, 'interval_mode', 'events' );
    [ Z_learn, contr_assemblies.learn_ts  ] = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', contr_assemblies.learn_trial_ints, 'interval_mode', 'window', 'fwhm', run_dt, 'binsize', run_dt, 'inclusion', inclusion  );
    Z_rest2 = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', hsync_rest2_ts_binned, 'interval_mode', 'events' );
    [ Z_probe2, contr_assemblies.probe2_ts  ] = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals',  contr_assemblies.probe2_ints, 'interval_mode', 'window', 'fwhm', run_dt, 'binsize', run_dt, 'inclusion', inclusion  );

    % consider removing low firing rate cells
    if params.filter_low_fr
        low_r_cells = any( [ sum(Z_rest1) ; sum(Z_rest2) ] < params.min_n_spk );
        Z_rest1( :, low_r_cells ) = []; Z_rest2( :, low_r_cells ) = [];
    end

    % run contrastive PCA
    [ assemblies, err_msg ] = run_gcpca( Z_rest2, Z_rest1, 'gcPCAversion', 4.1, 'nshuffle', nshuff );
    details.synchronyEvents = 'ripples';
    details.err_msg{1} = err_msg;


    %  succeeded in obtaining gcPCs with removal of low rate cells  
    if ~isempty( fieldnames( assemblies ) ) && params.filter_low_fr
        pyrid( low_r_cells ) = [];
        Z_probe1( :, low_r_cells ) = []; Z_learn( :, low_r_cells ) = []; Z_probe2( :, low_r_cells ) = [];
    end
        
    % variance explained by each gcPC per subsession
    disp('explained variance per subsession')
    nc = size( assemblies.gcPCA.X, 2 );
    zProbe1_cov = cov( zscore( Z_probe1 ) );
    zLearn_cov = cov( zscore( Z_learn ) );
    zProbe2_cov = cov( zscore( Z_probe2 ) );
    % probe1
    contr_assemblies.pcVar.d_probe1_gcpca = arrayfun(@(v) assemblies.gcPCA.X(:,v)' * zProbe1_cov * assemblies.gcPCA.X(:,v), 1:nc );
    [~,~, contr_assemblies.pcVar.d_probe1_pca] = pca( zscore( Z_probe1 ) );
    % learn
    contr_assemblies.pcVar.d_learn_gcpca = arrayfun(@(v) assemblies.gcPCA.X(:,v)' * zLearn_cov * assemblies.gcPCA.X(:,v), 1:nc );
    [~,~, contr_assemblies.pcVar.d_learn_pca] = pca( zscore( Z_learn ) );
    % probe2
    contr_assemblies.pcVar.d_probe2_gcpca = arrayfun(@(v) assemblies.gcPCA.X(:,v)' * zProbe2_cov * assemblies.gcPCA.X(:,v), 1:nc );
    [~,~, contr_assemblies.pcVar.d_probe2_pca] = pca( zscore( Z_probe2 ) );

    % store the remaining data
    contr_assemblies.assemblies = assemblies;
    contr_assemblies.hsync_rest1_ts_binned = hsync_rest1_ts_binned;
    contr_assemblies.hsync_rest2_ts_binned = hsync_rest2_ts_binned;
     
    save( outfile, 'contr_assemblies', 'params', 'outpath', 'details', 'e', 'pyrid', '-v7.3' );

catch e
    disp('Function returned due to error in loop ')
    save(outfile_err, 'e', 'basepath')
    
end

end