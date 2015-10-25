% go to main path if necessary
if regexpi(pwd(),'[\/\\]test$') % folder
    cd('..');
    init();
end