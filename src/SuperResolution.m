% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

classdef SuperResolution < handle
    %SUPERRESOLUTION Summary of this class goes here
    %   Detailed explanation goes here
    
    % TODO:
    % 1. remove lens distorion
    % 2. use dwt to merge images
    % 3. use optical flow to place pictures (areas) more accurate over each other
    % 4. remove errors from lens & sensor (i.e. dead pixels)
    
    properties
        sensorAlignment
        images = {}
        imageTransformations
        cameraParameters
    end
    
    methods
        function super_res = SuperResolution(sensorAlignment)
            if nargin == 1
                super_res.sensorAlignment = sensorAlignment;
            end
        end
        
        function estimateLensDistortion(super_res, varargin)
        % http://de.mathworks.com/help/vision/examples/evaluating-the-accuracy-of-single-camera-calibration.html
        %
        % doc about camera calibration and webcam support in Matlab
        % http://de.mathworks.com/help/vision/ug/single-camera-calibrator-app.html
        % 
        % get pattern pdf with "open checkerboardPattern.pdf"
        % (squareSize 23mm)
        
            disp('superResolution.estimateLensDistortion(...)');
            
            % TODO: only done for red channel at the moment
            data_prop = 'CDataRed';
            
            checkboard_images = Image.read(varargin{1}, super_res.sensorAlignment).(data_prop);
            for i = 2:length(varargin)
                image_data = Image.read(varargin{i}, super_res.sensorAlignment).(data_prop);
                checkboard_images = cat(4, checkboard_images, image_data);
            end
            
            % Detect calibration pattern.
            [imagePoints, boardSize] = detectCheckerboardPoints(checkboard_images);

            % Generate world coordinates of the corners of the squares.
            squareSizeInMM = 23;
            worldPoints = generateCheckerboardPoints(boardSize, squareSizeInMM);

            % Calibrate the camera.
            super_res.cameraParameters ...
                = estimateCameraParameters(...
                    imagePoints, worldPoints,...
                    'NumRadialDistortionCoefficients', 3,...
                    'EstimateTangentialDistortion', true);
            
            super_res.cameraParameters
                
            figure;
            showReprojectionErrors(super_res.cameraParameters);
            
            figure;
            for i = 1:size(checkboard_images,4)
                subplot(3, 2, i);
                imshow(checkboard_images(:,:,1,i));
                hold on;
                plot(imagePoints(:,1,i), imagePoints(:,2,i),'go');
                plot(super_res.cameraParameters.ReprojectedPoints(:,1,i), super_res.cameraParameters.ReprojectedPoints(:,2,i),'r+');
                legend('Detected Points','ReprojectedPoints');
                hold off;
            end

            figure;
            showExtrinsics(super_res.cameraParameters, 'CameraCentric');
        end
        
        function undistortImages(super_res)
            for image = super_res.images
                image{1}.undistortImage();
            end
        end
        
        function read(super_res, fileSelectors, sensorAlignment, filename_ext, mask_filename)
            disp('superResolution.read(...)');
            
            tic
            
            if ischar(fileSelectors)
                fileSelectors = {fileSelectors};
            end
            
            if nargin >= 5
                NaN_mask = imread(mask_filename);
                mask_true = max(NaN_mask(:));
                NaN_mask = NaN_mask ~= mask_true;
            end
            
            for i = 1:length(fileSelectors)
                fileSelector = fileSelectors{i};
                pathstr = fileparts(fileSelector);
                files = dir(fileSelector);
                for j = 1:length(files)
                    file = files(j);
                    filename = fullfile(pathstr,file.name);
                    if exist(filename,'file') == 2
                        image_idx = length(super_res.images)+1;
                        if nargin >= 5
                            super_res.images{image_idx} = Image.read(filename, sensorAlignment, filename_ext, NaN_mask);
                        elseif nargin >= 4
                            super_res.images{image_idx} = Image.read(filename, sensorAlignment, filename_ext);
                        elseif nargin >= 3
                            super_res.images{image_idx} = Image.read(filename, sensorAlignment);
                        else
                            super_res.images{image_idx} = Image.read(filename);
                        end
                    end
                end
            end
            
            toc
        end
        
        function extractFeatures(super_res, method)
            if nargin == 1
                method = Image.defaultValue('extractFeatures_method');
            end
            
            fprintf('superResolution.extractFeatures(%s)\n', method);
            
            tic
            
            for image = super_res.images
                image{1}.extractFeatures(method);
            end
            
            toc
        end
        
        function matchFeatures(super_res)
        %%
        % Match features by using their descriptors.
            disp('superResolution.matchFeatures()');
            
            tic
            
            length_images = length(super_res.images);
            
            super_res.imageTransformations = cell(length_images, length_images);
            
            %for i = 1:length_images-1
                i = 1;
                for j = i+1:length_images
                    super_res.imageTransformations{i,j} = ...
                        Image.getTransformation( ...
                            super_res.images{i}, super_res.images{j});
                end
            %end
            
            toc
        end
        
        function [joined_img,arranged_imgs] = joinImages(super_res, channel)
            disp('superResolution.joinImages(...)');
            
            if nargin < 2
                channel = 'red';
            end
            
            % guessed variables
            PSF = fspecial('gauss',[3 3],0.5);
            NOISE_VARIANCE = 0.0002;
            
            arranged_imgs = super_res.arrangeImages(...
                channel,PSF,NOISE_VARIANCE);
            
            arranged_imgs(:,:,1) = arranged_imgs(:,:,1);
            
            % adapte values relative to first image (if images were
            % recorded with different exposures)
            for i = 2:size(arranged_imgs,3)
                % TODO: do with every picture as picture 1 and 2 could be
                % both overexposed and then this function doesn't see any
                % problem.
%                 [arranged_imgs(:,:,1), arranged_imgs(:,:,i), correlation_curves{i}] = ...
%                     Image.adaptRelative(...
%                         arranged_imgs(:,:,1), arranged_imgs(:,:,i));
                CData1 = arranged_imgs(:,:,1);
                CData2 = arranged_imgs(:,:,i);
                idx_finite = isfinite(CData1) & isfinite(CData2);
                samples = [CData1(idx_finite) CData2(idx_finite)];
                [overexposure_CData1(i), overexposure_CData2(i), p_1_to_2(i,:), correlation_curves{i}] = ...
                    Image.findOverexposure(samples(:,1),samples(:,2));
            end
            
            figure
            hold on
            for i = 2:size(arranged_imgs,3)
                plot(correlation_curves{i}(:,1),(correlation_curves{i}(:,2) - p_1_to_2(i,2)) / p_1_to_2(i,1))
            end
            
%             amount_images = size(arranged_images,3);
%             
%             figure;
%             variance = var(arranged_images,0,3,'omitnan');
%             imshow(variance / max(variance(:)));
%             imwrite(im2uint16(variance),['test_data/results/variance-',int2str(amount_images),'-arranged_images.png']);
%             
%             figure;
%             histogram(variance(:));

            joined_img = SuperResolution.joinArrangedImagesMean(arranged_imgs);
        end
        
        function arranged_images = arrangeImages(super_res, channel, PSF, noise_variance)
        % scale relative to original image
        % TODO: use estimated_nsr only from reference image
            
            fprintf('SuperResolution.arrangeImages(%s)\n',channel);
            tic
            
            amount_images = length(super_res.images);
            reference_image = super_res.images{1};
            
            scale = round(amount_images^0.5);
            scale = min(2,scale);
            if ~isempty(reference_image.CDataBayer)
                % use half resolution for bayer data
                scale = scale / 2;
            end
            fprintf('  scale: %d\n', scale);
            
            image_size = scale * [reference_image.height reference_image.width];
            images_ref = imref2d(image_size);
            
            % Bayer: double internal scale for red and blue channel as they only
            % appear once in the bayer pattern
            if ~isempty(reference_image.CDataBayer) && ~strcmp(channel,'green')
                scale = 2 * scale;
            end
            
            % init arranged_images
            % Construction: [I1_data I1_alpha I2_data ....]
            % every second layer is the alpha channel for the previous data
            arranged_images = nan(...
                image_size(1),...
                image_size(2),...
                amount_images);

            for i = 1:amount_images
                image = super_res.images{i};
                
                switch channel
                    case {'red', 'blue'}
                        image_data = image.get('CData','channel',channel);
                        
                    case 'green'
                        if isempty(image.CDataBayer)
                            image_data = image.CDataGreen;
                            
                        else
                            image_data_size = 2 * size(image.CDataGreen1);
                            image_data = nan(image_data_size);

                            if nargin >= 3
                                % deblur channels green1, green2 using Wiener filter
                                % TODO: don't apply Wiener filter independent for
                                % green1 and green 2
                                green1 = image.CDataGreen1;
                                green2 = image.CDataGreen2;
                                estimated_nsr = noise_variance / var([image.CDataGreen1(:);image.CDataGreen2(:)]);
                                green1 = deconvwnr(green1, PSF, estimated_nsr);
                                green2 = deconvwnr(green2, PSF, estimated_nsr);
                            end

                            switch image.sensorAlignment
                                case {'bggr', 'rggb'}
                                %      B  G1   R  G1
                                %      G2 R    G2 B
                                    image_data(1:2:end, 2:2:end) = green1;
                                    image_data(2:2:end, 1:2:end) = green2;

                                case {'gbrg', 'grbg'}
                                %      G1 B    G1 R
                                %      R  G2   B  G2
                                    image_data(1:2:end, 1:2:end) = green1;
                                    image_data(2:2:end, 2:2:end) = green2;
                            end

                            % interpolate as imwarp doesn't support NaN values
                            % image_data = inpaint_nans(image_data);
                            % image_data
                        end
                end
                
                % not possible as a mask can be applied to the image during
                % reading and deconvwnr only allows finite values
%                 if strcmp(channel,'red') || strcmp(channel,'blue') || isempty(image.CDataBayer)
%                     if nargin >= 3
%                         % deblur image using Wiener filter
%                         estimated_nsr = noise_variance / var(image_data(:));
%                         image_data = deconvwnr(image_data, PSF, estimated_nsr);
%                     end
%                 end
                
                % scale image
                image_data = imresize(image_data, scale, 'cubic');
                
                % transpose all images relative to first image
                if i > 1
                    tform = super_res.imageTransformations{1,i};
                    
                    % transpose part of affine/projective matrix
                    if isa(tform,'projective2d')
                        % TODO: doesn't work yet
                        [1 0 0; 0 1 0; scale scale 1] * tform.T; % for test
                        tform.T(3,1:2) = [scale 0 0; 0 scale 0; 0 0 1] * tform.T(3,1:2);
                    else
                        tform.T(3,1:2) = scale * tform.T(3,1:2);
                    end
                    image_data = imwarp(image_data,tform,'Interp','cubic',...
                        'OutputView',images_ref,'FillValues',NaN);
                end

                arranged_images(:,:,i) = image_data;
            end
            
            toc
        end
        
        function writeImages(super_res, joined_image, arranged_images)
            disp('SuperResolution.writeImages(...)');
            
            global results_folder;
            
            % scale arranged images to 0..1
            [arranged_images,min_arranged_images,max_arranged_images] = ...
                normalize_img(arranged_images);
            
            % gamma value for writing results
            gamma = Image.defaultValue('gamma');
            
            reference_image = arranged_images(:,:,1);
            
            disp(['  write ',results_folder,'/2_reference_image.png']);
            imwrite(im2uint16(reference_image .^ (1/gamma)),[results_folder,'/2_reference_image.png']);
            
            original_image = super_res.images{1}.get('CData','channel','red');
            original_image = normalize_img(original_image,...
                min_arranged_images,max_arranged_images);
            original_image = imresize(original_image, size(reference_image));
            disp(['  write ',results_folder,'/3_original_image.png']);
            imwrite(im2uint16(original_image .^ (1/gamma)),[results_folder,'/3_original_image.png']);
            
            mkdir(results_folder,'/arranged_images/');
            for i = 1:size(arranged_images,3)
                disp(['  write ',results_folder,'/arranged_images/',int2str(i),'.png']);
                
                img_overexp_marked = image_mark_NaN(arranged_images(:,:,i) .^ (1/gamma));
                
                imwrite(im2uint16(img_overexp_marked),[results_folder,'/arranged_images/',int2str(i),'.png']);
            end
            
            disp(['  write ',results_folder,'/4_used_images.png']);
            used_images = sum(isfinite(arranged_images),3) / length(super_res.images);
            imwrite(im2uint16(used_images),[results_folder,'/4_used_images.png']);
            
            joined_image = normalize_img(joined_image,...
                min_arranged_images, max_arranged_images);
            
            disp(['  write ',results_folder,'/1_joined_image_mean.png']);
            imwrite(im2uint16(joined_image .^ (1/gamma)),[results_folder,'/1_joined_image_mean.png']);
        end
    end  
    
    methods(Static)
        function joined_image = joinArrangedImagesMean(arranged_images)
            disp('SuperResolution.joinArrangedImagesMean(...)');
            tic
            
            joined_image = nanmean(arranged_images, 3);
            
            toc
        end
        
        function joined_image = joinArrangedImagesMedian(arranged_images)
            disp('SuperResolution.joinArrangedImagesMedian(...)');
            tic
            
            joined_image = median(arranged_images,3,'omitnan');
            
            toc
        end
        
        function joined_image = joinArrangedImagesWavelet(arranged_images)
            disp('SuperResolution.joinArrangedImageWavelets(...)');
            tic
            
            WAVELET_LEVELS = 4;
            WAVELET = 'db4'; % see wfilters
            
            for i = size(arranged_images,3):-1:1
                [C(:,i),S] = wavedec2(arranged_images(:,:,i),WAVELET_LEVELS,WAVELET);
            end
            
            %C_joined = median(C,2,'omitnan');
            C_joined = nanmean(C,2);
            
            % reconstruct image from wavelets
            joined_image = waverec2(C_joined,S,WAVELET);
            
            toc
        end
        
        function generateTestImages(original_image_filename,path,amount_images,image_size)
            disp('SuperResolution.generateTestImages(...)');
            
            MAX_SCALE       = 0.5;% percent
            MAX_ROTATION    = 0.5;% degree
            MAX_TRANSLATION = 0.5;% percent
            
            NOISE_VARIANCE  = 0.0002;
            
            % read image and load in GPU
            original_image = im2double(imread(original_image_filename));
            
            % scale image so pure white and black will have also noise
            original_image = 0.9 * original_image + 0.05;
            
            original_image_size = size(original_image);
            
            reference_picture_ref = imref2d(image_size);
            
            default_scale = max(image_size(1:2) ./ original_image_size(1:2));
            default_translation = (image_size - default_scale * original_image_size(1:2)) / 2;
            
            for i = 1:amount_images
                if i == 1
                    scale = default_scale;
                    theta = 0;
                    translation = default_translation;
                else
                    scale = MAX_SCALE / 100 * (2 * rand(1) - 1) * default_scale + default_scale;
                    theta = MAX_ROTATION * (2 * rand(1) - 1);
                    translation = default_translation + (MAX_TRANSLATION / 100 * (2 * rand(1,2) - ones(1,2)) .* original_image_size(1:2));
                end
                
                % apply point spread function (PSF)
                image = imfilter(original_image,fspecial('disk',1/scale));
                
                % transform image relative to reference picture 1
                tform = affine2d(...
                    [scale*cosd(theta) -scale*sind(theta) 0;
                     scale*sind(theta) scale*cosd(theta)  0;
                     translation(2)    translation(1)     1]);
                image = imwarp(image,tform,'Interp','cubic','OutputView',reference_picture_ref);
                
                % add errors to image
                image = imnoise(image,'gaussian',0,NOISE_VARIANCE);
                
                imwrite(im2uint16(image),[path,'test_image_',int2str(i),'.png']);
            end
        end
    end
end

