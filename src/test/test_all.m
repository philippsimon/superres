% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

init();

global SAVE_FIGURES;
SAVE_FIGURES = 1;

global SHOW_FIGURES;
SHOW_FIGURES = 0;

folders = {...
    'test/Panasonic_LX100/exposure_series_1',...
    'test/Sony_A300/exposure_series_1',...
    'test/Sony_A300/exposure_series_2',...
    'test/Raspberry_Pi_NoIR_Camera/exposure_series_1',...
    'test/Raspberry_Pi_NoIR_Camera/exposure_series_2',...
    'test/Samsung_GT-S7580/nesquik',...
    'test/Panasonic_LX100/amount_images_needed',...
    'test/Basler_acA3800/exposure_series_1',...
};

test_Image();

for folder = folders
    if exist([folder{1},'/results'],'file') == 7 % folder
        warning([folder{1},' already processed'])
    else
        test_SuperResolution(folder{1});
    end
end