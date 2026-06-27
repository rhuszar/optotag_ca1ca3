function decodePos_gpr_hpc( index )
%
% decode position during the learning phase using Gaussian process regression
%
% decode from activity in contrastive subspaces 
%
% this function is written to run for each recording session on a high
% performance computing cluster, and requires that contrastive PCA has been
% run on all of these recording sessions

params = struct();
params.run_dt = 1 / 30;
params.speed_thresh = 5;
params.fwhm_run = 0.1;
params.N = 5000;
params.mindim = 5;

% path to folders containing outputs of contrastive PCA ( i.e., gcPCA_hpc.m )
infolder = '...';
% path to session folder on the HPC
datpath = '/gpfs/scratch/rh2618/rh_data';

index = str2double( index );
load('basepaths_ca1ca3.mat')
basepath = basepaths_all{ index };

% Basepepath on HPC
path_parts = strsplit(basepath,'\');
basepath_hpc = fullfile( datpath,path_parts{end-1},path_parts{end} );
basename = bz_BasenameFromBasepath( basepath_hpc );
disp(basepath_hpc)

% Define
st = dbstack;
funcname = st.name;
outpath = fullfile(datpath, [funcname '_OUT1']);
if ~exist(outpath, 'dir')
    mkdir(outpath);
end
fprintf('Running %s\n', funcname)

% Output folder for this basepath
outfile = fullfile( outpath, [ basename '.mat'] );
outfile_err = fullfile( outpath, [ basename '_ERROR.mat'] );


% output for contrast-defined assemblies
try

load( fullfile( infolder, [ basename '.mat'] ) )
cell_metrics = get_cellMetrics( basepath_hpc );
tracking = get_tracking( basepath_hpc );

indicators = getCellIndicators_cheese( basepath_hpc );
rewpos = get_rewardPosition( basepath_hpc, 'coordinates', 'animal' );

% inclusion criteria for running process_spikes
inclusion = struct();
inclusion(1).val = tracking.position.v_diff_smooth;
inclusion(1).ts = tracking.timestamps;
inclusion(1).thr = params.speed_thresh;
inclusion(2).val = double( ~isnan(tracking.position.x) & ~isnan(tracking.position.y) );
inclusion(2).ts = tracking.timestamps;
inclusion(2).thr = 0.9;

%% bin the data


if ~exist( 'pyrid', 'var')
    pyrid = find( indicators.hpc_pyr );
end

[ X, ts  ] = process_spikes( cell_metrics.spikes.times( pyrid ), 'intervals', contr_assemblies.learn_trial_ints, 'interval_mode', 'window', 'fwhm', params.fwhm_run, 'binsize', params.run_dt, 'inclusion', inclusion  );
X_z = zscore( X );

xpos = interp1( tracking.timestamps, tracking.position.x, ts  );
ypos = interp1( tracking.timestamps, tracking.position.y, ts  );

rew1_dist = vecnorm( [ xpos ypos ] - rewpos(1,:), 2, 2) ;
rew2_dist = vecnorm( [ xpos ypos ] - rewpos(2,:), 2, 2) ;

%% fit the GP model

sigA_inds = find( contr_assemblies.assemblies.sigA_gpc_select );
sigB_inds = flip( find( contr_assemblies.assemblies.sigB_gpc_select ) );
dim = min( [ length(sigA_inds), length(sigB_inds) ] ); 
if dim >= params.mindim
    
    % new subspace
    V_new = contr_assemblies.assemblies.gcPCA.X( :,  sigA_inds(1:dim) ) ;
    % old subspace
    V_old = contr_assemblies.assemblies.gcPCA.X( :,  sigB_inds(1:dim)  ) ;

else
    
    V_new = contr_assemblies.assemblies.gcPCA.X( :,  1:params.mindim ) ;
    V_old = contr_assemblies.assemblies.gcPCA.X( :,  ( end-params.mindim+1 ):end  ) ;
    dim = params.mindim;
end

% reconstructed neural activity
Z_learn_rec_new = ( X_z * V_new ) * V_new';
Z_learn_rec_old = ( X_z * V_old ) * V_old';

rng('shuffle')
ip = randperm( size( Z_learn_rec_new, 1 ) );
train_ids = sort( ip( 1:params.N ) );
test_ids = sort( ip( params.N+1:end ) );

% distance to reward
rew1_dist_test = rew1_dist( test_ids );
rew2_dist_test = rew2_dist( test_ids );

% train a separate model for each dimension

% new subspace
gpr_rec_new_x = fitrgp( Z_learn_rec_new( train_ids,: ), xpos( train_ids ) );
gpr_rec_new_y = fitrgp( Z_learn_rec_new( train_ids,: ), ypos( train_ids ) );

% old subspace
gpr_rec_old_x = fitrgp( Z_learn_rec_old( train_ids,: ), xpos( train_ids ) );
gpr_rec_old_y = fitrgp( Z_learn_rec_old( train_ids,: ), ypos( train_ids ) );


%% assess the model 

% new subspace
xpred_new_train = gpr_rec_new_x.predict( Z_learn_rec_new( train_ids,: ) );  % train
ypred_new_train = gpr_rec_new_y.predict( Z_learn_rec_new( train_ids,: ) );

xpred_new_test = gpr_rec_new_x.predict( Z_learn_rec_new( test_ids,: ) );    % test
ypred_new_test = gpr_rec_new_y.predict( Z_learn_rec_new( test_ids,: ) );

err_new_train = vecnorm( [ xpos( train_ids )  ypos( train_ids ) ] - [xpred_new_train ypred_new_train], 2, 2 );
err_new_test = vecnorm(  [ xpos( test_ids  )  ypos( test_ids ) ] -  [xpred_new_test  ypred_new_test], 2, 2 );

% old subspace
xpred_old_train = gpr_rec_old_x.predict( Z_learn_rec_old( train_ids,: ) );        % train
ypred_old_train = gpr_rec_old_y.predict( Z_learn_rec_old( train_ids,: ) );

xpred_old_test = gpr_rec_old_x.predict( Z_learn_rec_old( test_ids,: ) );     % test
ypred_old_test = gpr_rec_old_y.predict( Z_learn_rec_old( test_ids,: ) );

err_old_train = vecnorm( [ xpos( train_ids )  ypos( train_ids ) ] -  [ xpred_old_train ypred_old_train ], 2, 2 );
err_old_test = vecnorm( [ xpos( test_ids )  ypos( test_ids ) ] - [ xpred_old_test  ypred_old_test  ], 2, 2 );

% error as a function of trial
[~, trial_id] = InIntervals( ts( test_ids ), contr_assemblies.learn_trial_ints );
ntrials = length( contr_assemblies.learn_trial_ints  );

err_old_test_trial = [ arrayfun(@(x) mean( err_old_test( trial_id == x ) ), 1:ntrials )' arrayfun(@(x) sem( err_old_test( trial_id == x ) ), 1:ntrials )' ] ;
err_new_test_trial =  [ arrayfun(@(x) mean( err_new_test( trial_id == x ) ), 1:ntrials )' arrayfun(@(x) sem( err_new_test( trial_id == x ) ), 1:ntrials )' ];

save( outfile,'pyrid', 'train_ids', 'test_ids', 'Z_learn_rec_new', 'Z_learn_rec_old', 'gpr_rec_new_x', 'gpr_rec_new_y', ...
    'gpr_rec_old_x', 'gpr_rec_old_y','err_old_train', 'err_old_test', 'err_new_train', 'err_new_test', 'ts', 'trial_id', ...
    'err_old_test_trial', 'err_new_test_trial', 'dim','infolder','params', '-v7.3')


catch e
    disp('Function returned due to error in loop ')
    save(outfile_err, 'e', 'basepath')
    
end