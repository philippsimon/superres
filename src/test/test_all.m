% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

init();

global SAVE_FIGURES;
SAVE_FIGURES = 1;

global SHOW_FIGURES;
SHOW_FIGURES = 0;

% run Image tests
test_Image();

% run SuperResolution tests for all camera pictures
for folder = get_all_tests()
    if exist([folder{1},'/results'],'file') == 7 % folder
        warning([folder{1},' already processed'])
    else
        test_SuperResolution(folder{1});
    end
end