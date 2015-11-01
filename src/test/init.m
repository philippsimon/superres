% go to main path if necessary and run there init
if regexpi(pwd(),'[\/\\]test$') % folder
    cd('..');
end

init();