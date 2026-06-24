
%% concatenate and store the attribution timeseries 

outfile = 'attribution_ts.mat';

gpr_path = '...\decodePos_gpr_hpc_OUT';
attr_path = '...\getAttribution_hpc_zeroBase_OUT';
description = 'attribution timeseries per pyramidal neuron ; test data of GPR model';
for bp = 1:length( basepaths_all )

    [ basepath, basename ] = get_basepath( basepaths_all, bp );

    attr_path_full = fullfile( attr_path, basename );
    if ~exist( attr_path_full )
        continue
    end
    v = load( fullfile( gpr_path, [basename '.mat'] ) );
    fils = dir( attr_path_full );
    fils = { fils.name };
    fils(1:2) = [];
    eAttr_learn_x_all = []; eAttr_pre_x_all = [];
    eAttr_learn_y_all = []; eAttr_pre_y_all = [];
    % combine all the batched attribution outputs to produce a single time series
    for k = 1:length( fils )
        load( fullfile( attr_path_full, sprintf('_%d.mat', k) ) )
        eAttr_learn_x_all = add_to_mat( eAttr_learn_x_all, eAttr_learn_x' );
        eAttr_learn_y_all = add_to_mat( eAttr_learn_y_all, eAttr_learn_y' );
        eAttr_pre_x_all = add_to_mat( eAttr_pre_x_all, eAttr_pre_x' );
        eAttr_pre_y_all = add_to_mat( eAttr_pre_y_all, eAttr_pre_y' );
    end
    A_learn = abs( eAttr_learn_x_all + eAttr_learn_y_all ); 
    A_base = abs( eAttr_pre_x_all + eAttr_pre_y_all );
    ts_test = v.ts(v.test_ids);
    id_test = v.test_ids;
    
    % save datas
    save( fullfile( basepath, outfile ), 'A_learn', 'A_base', 'ts_test', 'gpr_path', 'attr_path', 'id_test', 'description', '-v7.3')

end