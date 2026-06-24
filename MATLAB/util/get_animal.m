function animal = get_animal( basepath )
% animal = get_animal( basepath )
% assumes cheeseboard basepath

    parts = strsplit( basepath, filesep );
    animal = parts{end-1};

end