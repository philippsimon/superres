function dirs = subfolders(dir_name)
    dirs = dir(dir_name);
    dirs = dirs([dirs(:).isdir]);
    dirs = {dirs(:).name};
    dirs(ismember(dirs,{'.','..'})) = [];
end