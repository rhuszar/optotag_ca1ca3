
function getAttribution_hpc( index, gpid )
%
% obtain attributions to a gaussian process position decoder
% written for HPC, since evaluating the attribution value of each cell at each timepoint is time consuming 
%
% this function is written to run on the cluster
% - index- indexes into recording session
% - gpid - indexes a group of timebins over which attributions are computed

% choice of reference baseline- zero was used in the published study
baseline = 'zero';  %  'zero' or 'mean'

% gpsize connections processed in a single call
gpsize = 160;
dat_ids = ( (gpid-1) * gpsize + 1) : gpid * gpsize;

% points to  files that store GPR decoders per each session (i.e., decodePos_gpr_hpc)
inpath = '...';
suff = strsplit( inpath, '_' ); suff = suff{ end };

% Path to sessions folder on the HPC
datpath = '...';

index = str2double(index);
load( 'basepaths_ca1ca3.mat' )


basepath = basepaths_all{ index };

% Basepepath on HPC
path_parts = strsplit(basepath,'\');
basepath_hpc = fullfile( datpath,path_parts{end-1},path_parts{end} );
basename = bz_BasenameFromBasepath( basepath_hpc );
disp(basepath_hpc)

% output path for this function
st = dbstack;
funcname = st.name;
outpath = fullfile(datpath, [funcname sprintf('_%sBase_%s', baseline, suff) ]);
if ~exist(outpath, 'dir')
    mkdir(outpath);
end
fprintf('Running %s\n', funcname)

% output for this specific session
outpath_session = fullfile( outpath, basename );
if ~exist(outpath_session, 'dir')
    mkdir(outpath_session);
end

% output folder for this basepath
outfile = fullfile( outpath_session, sprintf( '%d.mat', gpid ) );
outfile_err = fullfile( outpath_session, sprintf( '%d_ERROR.mat', gpid ) );

%% Initialize the outputs struct

% output for contrast-defined assemblies
out_struct = struct();
try
    
    load( fullfile( inpath, [basename '.mat'] ) )
    
    if strcmp( baseline, 'mean' )
        % note: data is zscored, so numerically, 0 is the mean
        x0_learn = zeros( 1, size( Z_learn_rec_learn, 2) );
        x0_pre = zeros( 1, size( Z_learn_rec_pre, 2) );
    elseif strcmp( baseline, 'zero' )
        x0_learn = min( Z_learn_rec_learn );
        x0_pre = min( Z_learn_rec_pre );
    end
    
    % store the expected values and variances of attributions, which are themselves Gaussians 
    eAttr_learn_x = [];
    vAttr_learn_x = [];
    eAttr_learn_y = [];
    vAttr_learn_y = [];

    eAttr_pre_x = [];
    vAttr_pre_x = [];
    eAttr_pre_y = [];
    vAttr_pre_y = [];

    train_data_learn = Z_learn_rec_learn( train_ids,: );
    test_data_learn =  Z_learn_rec_learn( test_ids,: );

    train_data_pre = Z_learn_rec_pre( train_ids,: );
    test_data_pre =  Z_learn_rec_pre( test_ids,: );
    
    % check if any datapoints to process
    if dat_ids(1) > length( test_ids )
        disp('First index exceeds number of test data...')
        return
    end

    tic
    prevtoc = toc;
    for dat_id = 1:length( dat_ids )
        
        if dat_ids( dat_id ) > length( test_ids )
            break
        end
        xtest_learn = test_data_learn( dat_ids( dat_id ), : );
        xtest_pre =   test_data_pre(   dat_ids( dat_id ), : );
        % get_attr(gp,Xt,xp,x0)
        % rest2 subspace- single neuron attributions
        [e,v] = get_attr( gpr_rec_learn_x, train_data_learn, xtest_learn, x0_learn );
        eAttr_learn_x = [ eAttr_learn_x e(:) ];
        vAttr_learn_x = [ vAttr_learn_x v(:) ];

        [e,v] = get_attr( gpr_rec_learn_y, train_data_learn, xtest_learn, x0_learn );
        eAttr_learn_y = [ eAttr_learn_y e(:) ];
        vAttr_learn_y = [ vAttr_learn_y v(:) ];

        % rest1 subspace- single neuron attributions
        [e,v] = get_attr( gpr_rec_pre_x, train_data_pre, xtest_pre, x0_pre );
        eAttr_pre_x = [ eAttr_pre_x e(:) ];
        vAttr_pre_x = [ vAttr_pre_x v(:) ];

        [e,v] = get_attr( gpr_rec_pre_y, train_data_pre, xtest_pre, x0_pre );
        eAttr_pre_y = [ eAttr_pre_y e(:) ];
        vAttr_pre_y = [ vAttr_pre_y v(:) ];
        
        % user feedback
        if mod(dat_id, 5) == 0
            thistoc = toc;
            fprintf('%d ; %.2fs elapsed\n', dat_id, thistoc-prevtoc)
            prevtoc = thistoc;
        end
    end
    t = toc;

    save( outfile, 'eAttr_learn_x', 'vAttr_learn_x', 'eAttr_learn_y', 'vAttr_learn_y', 'eAttr_pre_x',...
                   'vAttr_pre_x', 'eAttr_pre_y', 'vAttr_pre_y', 't', 'inpath' );

catch e
    disp('Function returned due to error in loop ')
    save(outfile_err, 'e', 'basepath')
    
end
