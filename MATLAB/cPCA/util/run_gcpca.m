function [ out, err_msg ] = run_gcpca(Za, Zb, varargin)
% Za- target data (nsamples X nfeatures)
% Zb- background data
%
% isolates components that pass the shuffle + higher variance explained in
% one versus the other task condition
% 
% only looks in some prctile of gcPCs for the significantly discriminating subspace
%
p=inputParser;
addParameter(p,'nshuffle',100,@isnumeric);  
addParameter(p,'sigthresh',2,@isnumeric);   % number of sd above mean
addParameter(p,'gcPCAversion',1,@isnumeric);% version 1 is just cPCA with alpha = 1
addParameter(p,'cumVarRat_thresh',1.03,@isnumeric);
addParameter(p,'prctl',20,@isnumeric);   % percentile of components among which to search for the significant ones

parse(p,varargin{:});
nshuffle = p.Results.nshuffle;
gcPCAversion = p.Results.gcPCAversion;
sigthresh = p.Results.sigthresh;
cumVarRat_thresh = p.Results.cumVarRat_thresh;
prctl = p.Results.prctl;

out = struct();
err_msg = [];

nfeatures = size( Za, 2 );
nf_subset = ceil( nfeatures / 100 * prctl );

% Run gcPCA
[B, S, C] = gcPCA( Za, Zb, gcPCAversion, 'Nshuffle', nshuffle );

% store gcPCA output
out.gcPCA.B = B;
out.gcPCA.S = S;
out.gcPCA.X = C;
if ~isfield( S, 'a_shuf' )
    err_msg = 'gcPCA cannot shuffle';
    warning( err_msg )
    return
end

cumvarRat_ab = abs( sum( S.a ) ./ sum( S.b ) );
cumvarRat_ba = abs( sum( S.b ) ./ sum( S.a ) );
cumvarRat = max( [cumvarRat_ab , cumvarRat_ba] );
if cumvarRat > cumVarRat_thresh
    err_msg = sprintf('Absolute cumulative variance is %.2f- exceeds threshold of %d\n', cumvarRat, cumVarRat_thresh );
    warning( err_msg )
    return
end

out.nf_subset = nf_subset;


%% A-specific components

% significant components for the A condition
out.sA_mu = mean( S.a_shuf, 2 );
out.sA_sd = std( S.a_shuf,[], 2 );
sigA_gpc = ( S.a - out.sA_mu) ./ out.sA_sd > sigthresh;

% select components that pass the shuffle and higher variance in A than B
sigA_gpc_select = sigA_gpc & S.a > S.b;
sigA_gpc_select( nf_subset+1:end ) = false;
out.sigA_gpc_select = sigA_gpc_select;

% number of significant components in the A condition
out.nc_A = sum( sigA_gpc_select );


%% B-specific components

out.sB_mu = mean( S.b_shuf, 2 );
out.sB_sd = std( S.b_shuf,[], 2 );
sigB_gpc = ( S.b - out.sB_mu) ./ out.sB_sd > sigthresh;

% select components that pass the shuffle and higher variance in B than A
sigB_gpc_select = sigB_gpc & S.b > S.a;
sigB_gpc_select( 1:(nfeatures-nf_subset) ) = false;
out.sigB_gpc_select = sigB_gpc_select;

% number of significant components in the B condition
out.nc_B = sum( sigB_gpc_select );


end