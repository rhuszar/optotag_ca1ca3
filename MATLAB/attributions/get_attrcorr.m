
%% 

gpr_out = 'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\tbr2\processed\decodePos_gpr_hpc_OUT'; % 'C:\Users\Roman\Documents\DATA\decodePos_gpr_hpc_OUT';
%gpr_out = 'C:\Users\Roman\Documents\DATA\decodePos_gpr_hpc_OUT_1';
pca_out = 'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\tbr2\processed\decodePos_gpr_hpc_OUT_1PCA'; %'C:\Users\Roman\Documents\DATA\decodePos_gpr_hpc_OUT_1PCA';
attrfile = 'attribution_ts.mat';
projd = 'projDiff.mat';


%%
basepaths_all = basepaths_saver_multiRegion;
basepaths_all_local = basepaths_saver_multiRegion( 'islocal', true );

saveflag = false;

% percentile bins
prcbin = 10;
nprcbins = 100 / prcbin;

% variables for storing data
z_dist = [];

% attribution correlations
ccAttr_learnSplit_all = [];   % split by subspace expressed
ccAttr_baseSplit_all = [];
ccAttr_learn_all = [];        % no splitting
ccAttr_base_all = [];

% convergence on interneurons
converg_new_all = [];   
converg_old_all  = [];
converg_all = [];
 
% subspace-reconstructed spiking activity
cc_learnSplit_all = [];  
cc_baseSplit_all = [];
cc_learn_all = [];
cc_base_all = [];

% decoding error, split by subspace expression
decErr_learnSplit = []; 
decErr_baseSplit = [];
decErr_learn = [];    % decoding error, all datapoints
decErr_base = [];
rewDist_all = [];
trialn_all = []; 

probe_id_all = [];
dist_all = [];
sbd_id_all = [];
dbd_id_all = [];
lat_class_all = [];
reg_all = [];
animal_all_pairs = [];
basename_all_pairs = [];

animal_all_testdata = [];
basename_all_testdata = [];

outfile = '...\attrcorr.mat';

for bp = 1:length( basepaths_all )

[ basepath, basename ] = get_basepath( basepaths_all, bp );
basename
clear pyrid
animal = get_animal( basepath );

c_a_path = fullfile( basepath, projd );
attr_path = fullfile( basepath, attrfile );
if ~exist( c_a_path ) || ~exist( attr_path )
    continue
end
load( c_a_path )


tracking = get_tracking( basepaths_all{bp} );
load( fullfile( basepath, 'chanMap.mat') )
load( fullfile( basepath, attrfile ) )
load( fullfile( basepath, [basename '.mono_res.cellinfo.mat'] ) )
load( fullfile( gpr_out, [basename '.mat']), 'err_learn_test', 'err_pre_test', 'Z_learn_rec_learn', 'Z_learn_rec_pre' )

% data reconstructed based on projections into the two subspaces
X_learn_rec = Z_learn_rec_learn( id_test, : );
X_pre_rec = Z_learn_rec_pre( id_test, : );

% session indices and indicators
sessints = getSession_ints( tracking );
cellindicators = getCellIndicators_cheese( basepath );

% new assembly ints
% new_int = diff( int_r')' > 0 & rs < .05;
% old_int = diff( int_r')' < 0 & rs < .05;
new_int = cellindicators.int_mod == 1;
old_int = cellindicators.int_mod == -1;
pv_int = cellindicators.classes.pv_id ;
pv_new_int = new_int & pv_int;
sst_int = cellindicators.classes.sst_id ;
sst_new_int = new_int & sst_int;

if exist( 'pyrid', 'var' )
    cellindicators.hpc_pyr = pyrid;
end
% int_id = find( cellindicators.hpc_int );
int_id = find( cellindicators.classes.int_id  );
new_int = new_int( int_id);
old_int = old_int ( int_id  );
pv_int = pv_int ( int_id  );
pv_new_int = pv_new_int( int_id );
sst_int = sst_int ( int_id  );
sst_new_int = sst_new_int( int_id );

% opto id
if strcmp( animal, 'iuvE14_cagCre_ef1aDIOchr2_8' )
    opto_id = cellindicators.ca1_tag_all | cellindicators.ca3_tag_all;
elseif strcmp( animal, 'iuvE15_camkiiCre_cagDIOchr2_1' )
    opto_id = false( length( cellindicators.maxWaveformCh1 ), 1 );
else
    opto_id = cellindicators.ca1_tag_sameShank | cellindicators.ca3_tag_sameShank;
end

% cell locations and distances between the neurons
cell_locations = [ xcoords( cellindicators.maxWaveformCh1( cellindicators.hpc_pyr  ) )' ycoords( cellindicators.maxWaveformCh1( cellindicators.hpc_pyr ) )' ];
distmat = get_estimatedCellDist( cell_locations );
converg = get_convergenceIndex( mono_res.pyr2int, cellindicators.hpc_pyr );

brainreg = double( cellindicators.ca1_id + 2.*cellindicators.ca3_id );

z_learn = projDiff.learn.z(id_test);
z_dist = add_to_mat(z_dist, {z_learn});

% split the data by subspace expression
edges = [ -inf prctile( z_learn, [ prcbin:prcbin:( 100-prcbin ) ] ) inf];
[~,~, ip] = histcounts( z_learn, edges );

nbins = max( unique(ip) );
ccAttr_learn = {}; ccAttr_base = {};
cc_learn = {}; cc_base = {};
decErr_l = {}; decErr_b = {}; decErr_pc = {};
% compute correlations conditioned on subspace expression
for b = 1:nbins
    
    % attribution correlation
    this_cc = get_cc( A_learn( b == ip, : ), 'subset', opto_id( cellindicators.hpc_pyr  ), 'uid', cellindicators.hpc_pyr, ...
                'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old} );
    ccAttr_learn{b} = this_cc.cc;
    this_cc = get_cc( A_base( b == ip, : ), 'subset', opto_id( cellindicators.hpc_pyr  ), 'uid', cellindicators.hpc_pyr,...
                'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old} );
    ccAttr_base{b} = this_cc.cc;
    % correlations of subspace-reconstructed spiking activity
    this_cc = get_cc( X_learn_rec( b == ip, : ), 'subset', opto_id( cellindicators.hpc_pyr  ), 'uid', cellindicators.hpc_pyr,...
                'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old } );
    cc_learn{b} = this_cc.cc;
    this_cc = get_cc( X_pre_rec( b == ip, : ), 'subset', opto_id( cellindicators.hpc_pyr  ), 'uid', cellindicators.hpc_pyr,...
                'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old } );
    cc_base{b} = this_cc.cc;
    % decoding error
    decErr_l{b} = err_learn_test( b == ip );
    decErr_b{b} = err_pre_test( b == ip );

end

% attribution correlations- subspace-conditioned 
ccAttr_learn = cell2mat( ccAttr_learn );
ccAttr_base = cell2mat( ccAttr_base );
ccAttr_learnSplit_all = add_to_mat( ccAttr_learnSplit_all, ccAttr_learn );
ccAttr_baseSplit_all = add_to_mat( ccAttr_baseSplit_all, ccAttr_base );

% correlations of subspace-reconstructed data
cc_learn = cell2mat( cc_learn );
cc_base = cell2mat( cc_base );
cc_learnSplit_all = add_to_mat( cc_learnSplit_all, cc_learn );
cc_baseSplit_all = add_to_mat( cc_baseSplit_all, cc_base );

% decoding error
decErr_learn = add_to_mat( decErr_learn, err_learn_test );     % all datapoints
decErr_base = add_to_mat( decErr_base, err_pre_test );
decErr_learnSplit =add_to_mat( decErr_learnSplit, decErr_l );  % subspace conditioned
decErr_baseSplit = add_to_mat( decErr_baseSplit, decErr_b );
% store animals / basenames for samples
animal_all_testdata = add_to_mat( animal_all_testdata, repmat({animal}, length( id_test ), 1) );
basename_all_testdata = add_to_mat( basename_all_testdata, repmat({basename}, length( id_test ), 1) );

% attribution correlation of all datapoints
ccAttr_learn = get_cc( A_learn, 'subset', opto_id( cellindicators.hpc_pyr  ), 'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old converg} );
ccAttr_base = get_cc( A_base, 'subset', opto_id( cellindicators.hpc_pyr  ), 'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old converg} );
ccAttr_learn_all = add_to_mat( ccAttr_learn_all, ccAttr_learn.cc );
ccAttr_base_all = add_to_mat( ccAttr_base_all, ccAttr_base.cc );
% correlation of subspace reconstructed activity
cc_learn = get_cc( X_learn_rec, 'subset', opto_id( cellindicators.hpc_pyr  ), 'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old converg ...
                                                                                                        converg_pv converg_sst converg_pv_new converg_sst_new} );
cc_base = get_cc( X_pre_rec, 'subset', opto_id( cellindicators.hpc_pyr  ), 'other_indicators', { brainreg( cellindicators.hpc_pyr ) cellindicators.probe_id( cellindicators.hpc_pyr ) distmat converg_new converg_old converg} );
cc_learn_all = add_to_mat( cc_learn_all, cc_learn.cc );
cc_base_all = add_to_mat( cc_base_all, cc_base.cc );

% convergence on to old and new 
converg_all = cat( 1, converg_new_all, cc_learn.other_indicator_pairs{ 6 } );
converg_new_all = cat( 1, converg_new_all, cc_learn.other_indicator_pairs{ 4 } );
converg_old_all = cat( 1, converg_old_all, cc_learn.other_indicator_pairs{ 5 } );
converg_pv_all = cat( 1, converg_pv_all, cc_learn.other_indicator_pairs{ 7 } );
converg_sst_all = cat( 1, converg_sst_all, cc_learn.other_indicator_pairs{ 8 } );
converg_pv_new_all = cat( 1, converg_pv_new_all, cc_learn.other_indicator_pairs{ 9 } );
converg_sst_new_all = cat( 1, converg_sst_new_all, cc_learn.other_indicator_pairs{ 10 } );


rewDist_all = cat( 1, rewDist_all, projDiff.learn.rewDist( id_test ) );
trialn_all = cat( 1, trialn_all, projDiff.learn.trialn( id_test ) );


% probe
probeid =  sum( this_cc.other_indicator_pairs{2}, 2);
probe_id_all = add_to_mat( probe_id_all, probeid );
% pair distance
dist_all = cat( 1, dist_all, this_cc.other_indicator_pairs{3} );
% region
regid = sum( this_cc.other_indicator_pairs{1}, 2 );
reg_all = cat( 1, reg_all, regid );

% add latency class for a more stringent analysis
if isfield( cellindicators, 'opto_latency_class' )
    lat_class_all = cat( 1, lat_class_all, cellindicators.opto_latency_class( this_cc.uid_pair ) );
else
    lat_class_all = cat( 1, lat_class_all, nan( size( this_cc.uid_pair ) ) );
end

% same birthdate
s = sum( this_cc.subset_pair_id, 2);
sbd_id_all = cat( 1, sbd_id_all, s == 2 );
dbd_id_all = cat( 1, dbd_id_all, s == 1 );

npairs = length( this_cc.cc );
animal_all_pairs = cat( 1, animal_all_pairs, repmat({animal}, npairs, 1) );
basename_all_pairs = cat( 1, basename_all_pairs, repmat({basename}, npairs, 1) );

end


%%
if saveflag
    save( outfile, 'basename_all_pairs', 'animal_all_pairs', 'reg_all', 'dbd_id_all', 'sbd_id_all', 'dist_all', 'probe_id_all', ...
        'decErr_learnSplit', 'decErr_baseSplit', 'decErr_learn', 'decErr_base', 'rewDist_all', 'trialn_all', ...
        'cc_learnSplit_all', 'cc_baseSplit_all', 'cc_learn_all', 'cc_base_all', 'converg_new_all', 'converg_old_all', 'converg_all', ...
          'ccAttr_learnSplit_all', 'ccAttr_baseSplit_all', 'ccAttr_learn_all', 'ccAttr_base_all', 'z_dist', 'animal_all_testdata', 'basename_all_testdata', 'description' )
end

%%