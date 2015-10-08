% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

function test_SuperResolution(folder)
    %opengl('save','software');

    init();
    
    global RESULTS_FOLDER;
    RESULTS_FOLDER = [folder,'/results'];

    close all;
    super_res = SuperResolution();
    
    function run()
        disp(['create clean results folder ',RESULTS_FOLDER]);
        try
            rmdir(RESULTS_FOLDER,'s');
        catch err
        end
        mkdir(RESULTS_FOLDER);
        
        %super_res.images{1}.show('red');
        %super_res.undistortImages();
        super_res.extractFeatures();
        super_res.matchFeatures();
        [arranged_imgs, weights_arranged_imgs, arranged_imgs_orig] = ...
            super_res.adaptImagesAndCalcWeights('red', true);
        joined_image = SuperResolution.joinArrangedImages(...
            arranged_imgs, weights_arranged_imgs);
        super_res.writeImages(joined_image, arranged_imgs,...
            weights_arranged_imgs, arranged_imgs_orig);
    end

    % SuperResolution.generateTestImages(...
    %    '../misc/ISO_12233-reschart.png','test/',10,[600 800]);
    
    switch folder
        % default test behavior from RAW camera data
        case {...
            'test/Panasonic_LX100/exposure_series_1',...
            'test/Sony_A300/exposure_series_1',...
            'test/Sony_A300/exposure_series_2',...
            'test/Raspberry_Pi_NoIR_Camera/exposure_series_1'...
            'test/Raspberry_Pi_NoIR_Camera/exposure_series_2'...
        }
            super_res.read({[folder,'/*.*']});
            run();
            
        case 'test/Samsung_GT-S7580/nesquik'
            % super res test - no changing expsosure time
            % or gbrg?
            super_res.read({[folder,'/*.jpg']}, 'bggr');
            run();
            
        case 'test/Panasonic_LX100/amount_images_needed'
            % amount images needed for super resolution
            super_res.read({[folder,'/P1020520.RW2']});
            
            original_RESULTS_FOLDER = RESULTS_FOLDER;
            
            for i = 1:9
                RESULTS_FOLDER = [original_RESULTS_FOLDER,'/',num2str(i+1)];
                image_idx = length(super_res.images)+1;
                super_res.images{image_idx} = Image.read(...
                    [folder,'/P102052',num2str(i),'.RW2']);
                run();
                movefile([RESULTS_FOLDER,'/1_joined_image_mean.png'],...
                    [original_RESULTS_FOLDER,'/joined_image_from_',...
                    num2str(i+1),'_images.png']);
            end
            
        otherwise
            error(['Test for folder ',folder,' not found'])
    end
end