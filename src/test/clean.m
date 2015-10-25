clear

init()

for test_dir = get_all_tests()
    try
        rmdir([test_dir{1},'/results'],'s');
    catch
    end
end