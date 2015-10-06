function SuperResolution_test(folder)
    %opengl('save','software');

    init();
    
    global results_folder;
    results_folder = [folder,'/results'];

    close all;
    super_res = SuperResolution();
    
    function run()
        disp(['create clean results folder ',results_folder]);
        try
            rmdir(results_folder,'s');
        catch err
        end
        mkdir(results_folder);
        
        %super_res.images{1}.show('red');
        %super_res.undistortImages();
        super_res.extractFeatures();
        super_res.matchFeatures();
        [joined_image, arranged_images] = super_res.joinImages();
        super_res.writeImages(joined_image, arranged_images);
    end

    % SuperResolution.generateTestImages('../misc/ISO_12233-reschart.png','test/',10,[600 800]);
    
    switch folder
        % default test behavior from RAW camera data
        case {...
            'test/Panasonic_LX100/exposure_series_1',...
            'test/Sony_A300/exposure_series_1',...
            'test/Sony_A300/exposure_series_2',...
            'test/RPiTelecine/exposure_series_1'...
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
            
            original_results_folder = results_folder;
            
            for i = 1:9
                results_folder = [original_results_folder,'/',num2str(i+1)];
                image_idx = length(super_res.images)+1;
                super_res.images{image_idx} = Image.read(...
                    [folder,'/P102052',num2str(i),'.RW2']);
                run();
                movefile([results_folder,'/1_joined_image_mean.png'],...
                    [original_results_folder,'/joined_image_from_',num2str(i+1),'_images.png']);
            end
            
        otherwise
            error(['Test for folder ',folder,' not found'])
    end
end