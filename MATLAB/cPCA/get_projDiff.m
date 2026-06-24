


%% variance explained version of the analysis
% extract Fisher z-scores capturing the relative variance capture by the
% two subspaces

dt_run = 1/30;
fwhm = .1;
speed_thresh = 5;
rewRadius = 10;
mindim = 5;
distbins = 0:5:30;

outfile = 'projDiff.mat';
mindim_sess = false( length( basepaths_all ), 1 );
for bp = 1:length( basepaths_all )
    
    [ basepath, basename ] = get_basepath( basepaths_all, bp );
    animal = get_animal( basepath );
    clear pyrid
    basename
    c_a_path = fullfile( basepath, 'contrastive_assemblies.mat' );
    if ~exist( c_a_path )
        continue
    end
    load( c_a_path )
    load( fullfile( basepath, [basename '.cell_metrics.cellinfo.mat'] ) )
    load( fullfile( basepath, [basename '.tracking.behavior.mat'] ) )

    indicators = getCellIndicators_cheese( basepath );
    sessints = getSession_ints( tracking );
    
    % only include data when animal runs
    inclusion = struct();
    inclusion.val = tracking.position.v_diff_smooth;
    inclusion.ts = tracking.timestamps;
    inclusion.thr = speed_thresh;
    if ~exist( 'pyrid', 'var' )
        pyrid = find( indicators.hpc_pyr );
    end
    ncells = length( pyrid );
    ndim = min( [ contr_assemblies.assemblies.nc_A, contr_assemblies.assemblies.nc_B ] );
    % dimensionality of shuffled data too low... skip
    if ndim < 2
        ndim = mindim;
        sigA = 1:ndim; sigA = sigA( : );
        sigB = ( ncells-ndim+1 ):ncells; sigB = sigB( : );

        mindim_sess( bp ) = true;
    else
        % significant components
        sigA = find( contr_assemblies.assemblies.sigA_gpc_select ); 
        sigB = flip( find( contr_assemblies.assemblies.sigB_gpc_select ) );
        % ensure equal dimensionality of the two spaces
        sigA = sigA( 1:ndim );
        sigB = sigB( 1:ndim );

    end

    rewPos = get_rewardPosition( basepath, 'coordinates', 'animal' );
    % bin / smooth the data
    [ X_pre, ts_pre ] = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', sessints.pre_mazeints(1,:), ...
                                    'interval_mode', 'window', 'binsize', dt_run, 'fwhm', fwhm, 'inclusion', inclusion );
    [ X_learn, ts_learn ] = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', sessints.learn_mazeints, ...
                                    'interval_mode', 'window', 'binsize', dt_run, 'fwhm', fwhm, 'inclusion', inclusion );
    [ X_post, ts_post ] = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', sessints.post_mazeints(1,:), ...
                                    'interval_mode', 'window', 'binsize', dt_run, 'fwhm', fwhm, 'inclusion', inclusion );
    
    % interpolate to bout data
    bout_data = interpToBouts( basepath, ts_learn );
    % pre / post position data
    pos_pre = [ interp1( tracking.timestamps, tracking.position.x, ts_pre ) interp1( tracking.timestamps, tracking.position.y, ts_pre ) ];
    pos_post = [ interp1( tracking.timestamps, tracking.position.x, ts_post ) interp1( tracking.timestamps, tracking.position.y, ts_post ) ];

    % distances away from reward
    rew_dist_pre = [ vecnorm( pos_pre- rewPos(1,:), 2, 2 )  vecnorm( pos_pre- rewPos(2,:), 2, 2 ) ];
    rew_dist_learn = [ vecnorm( bout_data.pos - rewPos(1,:), 2, 2 )  vecnorm( bout_data.pos - rewPos(2,:), 2, 2 ) ];
    rew_dist_post = [ vecnorm( pos_post- rewPos(1,:), 2, 2 )  vecnorm( pos_post- rewPos(2,:), 2, 2 ) ];

    % bin distances, separate bouts that end in different rewards
    [~, ~, rew1_bininds] = histcounts( rew_dist_learn( :, 1 ), distbins );
    rew1_bininds(~( bout_data.rewID == 1 )) = 0;  % remove data that are not behavior ending with reward1 
    [~, ~, rew2_bininds] = histcounts( rew_dist_learn( :, 2 ), distbins );
    rew2_bininds(~( bout_data.rewID == 2 )) = 0;

    rew_bininds_learn = rew1_bininds + rew2_bininds;
    rew_proximity_pre = InIntervals( rew_dist_pre(:,1), [ 0 rewRadius ] ) | InIntervals( rew_dist_pre(:,2), [ 0 rewRadius ] ) ;
    rew_proximity_post = InIntervals( rew_dist_post(:,1), [ 0 rewRadius ] ) | InIntervals( rew_dist_post(:,2), [ 0 rewRadius ] ) ;

    % mean / sd of running data
    m = nanmean( [ X_pre ; X_learn ; X_post ] );
    sd = nanstd( [ X_pre ; X_learn ; X_post ] );
    
    % z score running data
    X_pre_z = ( X_pre - m ) ./ sd;
    X_learn_z = ( X_learn - m ) ./ sd;
    X_post_z = ( X_post - m ) ./ sd;
    % whatever is nan, set to zero- won't contribute to projections
    X_pre_z( isnan( X_pre_z ) ) = 0;
    X_learn_z( isnan( X_learn_z ) ) = 0;
    X_post_z( isnan( X_post_z ) ) = 0;


    V_learn = contr_assemblies.assemblies.gcPCA.X( :, sigA );
    V_base = contr_assemblies.assemblies.gcPCA.X( :, sigB );
    
    % PRE PROBE
    % reconstruct activity based on subspace projection
    X_l_hat = ( X_pre_z * V_learn ) * V_learn' ;
    X_b_hat = ( X_pre_z * V_base ) * V_base' ;

    % reconstruction error
    recErr_l = vecnorm( X_pre_z - X_l_hat, 2, 2 );
    recErr_b = vecnorm( X_pre_z - X_b_hat, 2, 2 );
    
    % r squuarred variance explained by each subspace
    r2_l_pre = 1 - ( recErr_l.^2 ) ./ vecnorm(X_pre_z, 2, 2).^2;
    r2_b_pre = 1 - ( recErr_b.^2 ) ./ vecnorm(X_pre_z, 2, 2).^2;
    % square root of r squared
    r2_l_sq = sqrt( r2_l_pre );
    r2_b_sq = sqrt( r2_b_pre );
    z_pre = nan( length( r2_l_sq ), 1 );
    % fisher z transform
    for k = 1:length( r2_l_sq )
        [~, z_pre(k)] = compare_correlation_coefficients( r2_l_sq(k), r2_b_sq(k), ncells, ncells);
    end
        
    % POST PROBE
    % reconstruct activity based on subspace projection
    X_l_hat = ( X_post_z * V_learn ) * V_learn' ;
    X_b_hat = ( X_post_z * V_base ) * V_base' ;

    % reconstruction error
    recErr_l = vecnorm( X_post_z - X_l_hat, 2, 2 );
    recErr_b = vecnorm( X_post_z - X_b_hat, 2, 2 );
    
    % variance explained by each subspace
    r2_l_post = 1 - ( recErr_l.^2 ) ./ vecnorm(X_post_z, 2, 2).^2;   % in the denominator we don't subtract anything since 
    r2_b_post = 1 - ( recErr_b.^2 ) ./ vecnorm(X_post_z, 2, 2).^2;
    % square root of r squared- correlation coefficient
    r2_l_sq = sqrt( r2_l_post );
    r2_b_sq = sqrt( r2_b_post );
    z_post = nan( length( r2_l_sq ), 1 );
    % fisher z transform
    for k = 1:length( r2_l_sq )
        [~, z_post(k)] = compare_correlation_coefficients( r2_l_sq(k), r2_b_sq(k), ncells, ncells);
    end

    % LEARNING
    % reconstruct activity based on subspace projection
    X_learn_hat = ( X_learn_z * V_learn ) * V_learn' ;
    X_base_hat = ( X_learn_z * V_base ) * V_base' ;

    % reconstruction error
    recErr_l = vecnorm( X_learn_z - X_learn_hat, 2, 2 );
    recErr_b = vecnorm( X_learn_z - X_base_hat, 2, 2 );
    
    % variance explained by each subspace
    r2_l_learn = 1 - ( recErr_l.^2 ) ./ vecnorm(X_learn_z, 2, 2).^2;
    r2_b_learn = 1 - ( recErr_b.^2 ) ./ vecnorm(X_learn_z, 2, 2).^2;
    % square root of r squared
    r2_l_sq = sqrt( r2_l_learn );
    r2_b_sq = sqrt( r2_b_learn );
    z_learn = nan( length( r2_l_sq ), 1 );
    % fisher z transform
    for k = 1:length( r2_l_sq )
        [~, z_learn(k)] = compare_correlation_coefficients( r2_l_sq(k), r2_b_sq(k), size( X_learn_z, 2), size( X_learn_z, 2) );
    end

    % average per trial for learning subsession
    ntrials = size( sessints.learn_mazeints, 1 );
    % trial ids 
    trial_ip = arrayfun(@(x) InIntervals( ts_learn, sessints.learn_mazeints(x,:) ), 1:ntrials , 'UniformOutput', false);
    z_learn_trial = cellfun(@(x) nanmean( z_learn(x) ), trial_ip );
    r2_l_learn_trial = cellfun(@(x) nanmean( r2_l_learn(x) ), trial_ip );
    r2_b_learn_trial = cellfun(@(x) nanmean( r2_b_learn(x) ), trial_ip );
    
    % store the outputs
    projDiff = struct();
    
    % pre subsession
    projDiff.pre.z = z_pre;
    projDiff.pre.ts = ts_pre;
    projDiff.pre.rewDist = rew_proximity_pre;

    % learning subsession
    % indicators
    projDiff.learn.trialn = bout_data.trialn;
    projDiff.learn.rewID = bout_data.rewID;
    projDiff.learn.rewBout = bout_data.rewBout;
    projDiff.learn.rewDist = rew_bininds_learn;
    % projections
    projDiff.learn.z = z_learn;
    projDiff.learn.z_trial = z_learn_trial;
    projDiff.learn.r2_l_trial = r2_l_learn_trial;
    projDiff.learn.r2_b_trial = r2_b_learn_trial;
    projDiff.learn.ts = ts_learn;
    
    % post subsession
    projDiff.post.z = z_post;
    projDiff.post.ts = ts_post;
    projDiff.post.rewDist = rew_proximity_post;

    % save 
    save( fullfile(basepath, outfile), 'projDiff', 'ndim', 'ncells', 'pyrid', '-v7.3' )
    clear pyrid

end