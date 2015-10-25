function dirs = get_all_tests()
    init();
    
    dirs = {};
    
    for camera_dir_name = subfolders('test')
        for test_dir_name = subfolders(['test/',camera_dir_name{1}])
            dirs{end+1} = ['test/',camera_dir_name{1},'/',test_dir_name{1}];
        end
    end
end

