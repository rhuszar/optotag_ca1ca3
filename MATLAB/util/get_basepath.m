function [basepath, basename] = get_basepath( basepaths, ip )

    basepath = basepaths{ip};
    basename = bz_BasenameFromBasepath( basepath );

end