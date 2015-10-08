% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

classdef Image < handle
    %IMAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sensorAlignment
        
        metaInfo
        width
        height
        
        colorSpace
        
        cameraParameters
        
        CDataBayer
        CData
        
        CDataRed
        CDataRedNormalized
        CDataGreen  % RGB green channel image
        CDataGreenNormalized
        CDataGreen1 % Bayer first green pixel image
        CDataGreen1Normalized
        CDataGreen2 % Bayer second green pixel image
        CDataGreen2Normalized
        CDataBlue
        CDataBlueNormalized
        
        featuresRed
        featuresGreen
        featuresGreen1
        featuresGreen2
        featuresBlue
        
        validPtsRed
        validPtsGreen
        validPtsGreen1
        validPtsGreen2
        validPtsBlue
    end
    
    methods
        function img = Image(CData, sensorAlignment)
            if nargin >= 1
                image_size = size(CData);
                img.height = image_size(1);
                img.width  = image_size(2);
                if nargin == 1
                    img.CData = CData;
                else
                    img.CDataBayer = CData;
                    img.sensorAlignment = sensorAlignment;
                    img.separateBayerChannels();
                end
            end
        end
        
        function undistortImage(img)
            if isempty(img.cameraParameters)
                error('image.cameraParameters not set');
            end
            
            img.CDataRed = undistortImage(img.CDataRed, img.cameraParameters);
        end
        
        function extractFeatures(img, method)
            if nargin == 1
                method = Image.defaultValue('extractFeatures_method');
            end
            
            fprintf('image.extractFeatures(%s)\n', method);
            
            map = containers.Map();
            if isempty(img.CDataBayer)
                map('red')    = {'featuresRed',    'validPtsRed'};
                map('green')  = {'featuresGreen',  'validPtsGreen'};
                map('blue')   = {'featuresBlue',   'validPtsBlue'};
            else
                map('red')    = {'featuresRed',    'validPtsRed'};
                map('green1') = {'featuresGreen1', 'validPtsGreen1'};
                map('green2') = {'featuresGreen2', 'validPtsGreen2'};
                map('blue')   = {'featuresBlue',   'validPtsBlue'};
            end
            for prop = map.keys
                feature_pts = map(prop{1});
                img_CData = img.get('CData',...
                    'channel',prop{1},...
                    'normalized',true);

                switch method
                    case 'SURF'
                        points = detectSURFFeatures(img_CData);
                        [features_, valid_points] = extractFeatures(...
                            img_CData, points);

                    case 'MSER'
                        points = detectMSERFeatures(img_CData);
                        [features_, valid_points] = extractFeatures(...
                            img_CData, points, 'Upright', true);

                    otherwise
                        error('method not supported');
                end

                if isempty(features_)
                    error('no features found');
                end

%                     figure; imshow(image_arr); hold on;
%                     plot(valid_points);

                img.(feature_pts{1}) = features_;
                img.(feature_pts{2}) = valid_points;

            end
        end
        
        function separateBayerChannels(img)
        % Split Bayer in separate Channels
        
            img.CDataRed    = nan(img.height/2,img.width/2);
            img.CDataGreen1 = nan(img.height/2,img.width/2);
            img.CDataGreen2 = nan(img.height/2,img.width/2);
            img.CDataBlue   = nan(img.height/2,img.width/2);
            
            switch img.sensorAlignment
            case 'gbrg'
                %       j%2=0 j%2=1
                % i%2=0   G1   B
                % i%2=1   R    G2
                img.CDataGreen1 = img.CDataBayer(1:2:img.height,1:2:img.width);
                img.CDataBlue   = img.CDataBayer(1:2:img.height,2:2:img.width);
                img.CDataRed    = img.CDataBayer(2:2:img.height,1:2:img.width);
                img.CDataGreen2 = img.CDataBayer(2:2:img.height,2:2:img.width);
                
            case 'grbg'
                %       j%2=0 j%2=1
                % i%2=0   G1   R
                % i%2=1   B    G2
                img.CDataGreen1 = img.CDataBayer(1:2:img.height,1:2:img.width);
                img.CDataRed    = img.CDataBayer(1:2:img.height,2:2:img.width);
                img.CDataBlue   = img.CDataBayer(2:2:img.height,1:2:img.width);
                img.CDataGreen2 = img.CDataBayer(2:2:img.height,2:2:img.width);
                
            case 'bggr'
                %       j%2=0 j%2=1
                % i%2=0   B    G1
                % i%2=1   G2   R
                img.CDataBlue   = img.CDataBayer(1:2:img.height,1:2:img.width);
                img.CDataGreen1 = img.CDataBayer(1:2:img.height,2:2:img.width);
                img.CDataGreen2 = img.CDataBayer(2:2:img.height,1:2:img.width);
                img.CDataRed    = img.CDataBayer(2:2:img.height,2:2:img.width);
                
            case 'rggb'
                %       j%2=0 j%2=1
                % i%2=0   R    G1
                % i%2=1   G2   B
                img.CDataRed    = img.CDataBayer(1:2:img.height,1:2:img.width);
                img.CDataGreen1 = img.CDataBayer(1:2:img.height,2:2:img.width);
                img.CDataGreen2 = img.CDataBayer(2:2:img.height,1:2:img.width);
                img.CDataBlue   = img.CDataBayer(2:2:img.height,2:2:img.width);
                
            otherwise
                error('not supported image.sensorAlignment');
            end
        end
        
        function demosaic(img)
        % Demosaicing
        % from https://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
            img.CData = im2double(demosaic(im2uint16(img.CDataBayer),img.sensorAlignment));
            %image.colorSpace = 'cam';
        end
        
        function brigthnessCorrection(img)
        % Brigthness Correction
        % from https://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
            grayim = rgb2gray(img.CData);
            grayscale = 0.25 / mean(grayim(:));
            img.CData = min(1,img.CData*grayscale);
        end
        
        function whitebalanceBayer(img,wb_multipliers)
            %  White Balancing
            % from https://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
            mask = wbmask(...
                size(img.CDataBayer,1),...
                size(img.CDataBayer,2),...
                wb_multipliers,...
                img.sensorAlignment);
            
            img.CDataBayer = img.CDataBayer .* mask;
        end
        
        function whitebalance(img)
            % whitebalance test
            %mult = illuminantEstimator(image.CData, [true true], 20);
            %image.CData(1) = image.CData(1) / mult(1);
            %image.CData(2) = image.CData(2) / mult(2);
            %image.CData(3) = image.CData(3) / mult(3);
            
            img.CData = whitebalance(img.CData);
        end
        
        function transformRGBtoBayer(img, sensorAlignment)
            img.sensorAlignment = sensorAlignment;
            img.CDataBayer = zeros(size(img.CData,1),size(img.CData,2));
            
            switch sensorAlignment
            case 'rggb'
                img.CDataBayer(1:2:end,1:2:end) = img.CData(1:2:end,1:2:end,1); % r
                img.CDataBayer(2:2:end,2:2:end) = img.CData(2:2:end,2:2:end,2); % b
                img.CDataBayer(1:2:end,2:2:end) = img.CData(1:2:end,2:2:end,3); % g1
                img.CDataBayer(2:2:end,1:2:end) = img.CData(2:2:end,1:2:end,3); % g2
            case 'bggr'
                img.CDataBayer(2:2:end,2:2:end) = img.CData(2:2:end,2:2:end,1); % r
                img.CDataBayer(1:2:end,1:2:end) = img.CData(1:2:end,1:2:end,2); % b
                img.CDataBayer(1:2:end,2:2:end) = img.CData(1:2:end,2:2:end,3); % g1
                img.CDataBayer(2:2:end,1:2:end) = img.CData(2:2:end,1:2:end,3); % g2
            case 'grbg'
                img.CDataBayer(1:2:end,2:2:end) = img.CData(1:2:end,2:2:end,1); % r
                img.CDataBayer(2:2:end,1:2:end) = img.CData(2:2:end,1:2:end,2); % b
                img.CDataBayer(1:2:end,1:2:end) = img.CData(1:2:end,1:2:end,3); % g1
                img.CDataBayer(2:2:end,2:2:end) = img.CData(2:2:end,2:2:end,3); % g2
            case 'gbrg'
                img.CDataBayer(2:2:end,1:2:end) = img.CData(2:2:end,1:2:end,1); % r
                img.CDataBayer(1:2:end,2:2:end) = img.CData(1:2:end,2:2:end,2); % b
                img.CDataBayer(1:2:end,1:2:end) = img.CData(1:2:end,1:2:end,3); % g1
                img.CDataBayer(2:2:end,2:2:end) = img.CData(2:2:end,2:2:end,3); % g2
            end
            
            img.CData = [];
        end
        
        function value = get(img, prop, varargin)
        % image_data = get('CData', channel?: String, ('normalized', [Boolean=false])?)
        
            p = inputParser;
        
            validProp = {'CData'};
            checkProp = @(x) any(validatestring(x,validProp));
            addRequired(p,'prop',checkProp);
        
            switch prop
                case 'CData'
                    defaultChannel = '';
                    validChannels = {'','bayer','red','green','green1','green2','blue'};
                    checkChannel = @(x) any(validatestring(x,validChannels));
                    addParameter(p,'channel',defaultChannel,checkChannel);

                    defaultNormalized = false;
                    addParameter(p,'normalized',defaultNormalized,@islogical);
            end
            
            parse(p,prop,varargin{:});
            
%             disp(['prop: ', p.Results.prop])
% 
%             if ~isempty(fieldnames(p.Unmatched))
%                disp('Extra inputs:')
%                disp(p.Unmatched)
%             end
%             if ~isempty(p.UsingDefaults)
%                disp('Using defaults: ')
%                disp(p.UsingDefaults)
%             end

            switch p.Results.prop
                case 'CData'
                    if ~isempty(p.Results.channel)
                        Channel = [upper(p.Results.channel(1)),p.Results.channel(2:end)];
                        prop = [prop,Channel];
                    end
                    if p.Results.normalized
                        prop_normalized = [prop,'Normalized'];
                        if isempty(img.(prop_normalized))
                            image_data = img.(prop);
                            image_data_min = min(image_data(:));
                            image_data_max = max(image_data(:));
                            img.(prop_normalized) = (img.(prop) - image_data_min) / (image_data_max - image_data_min);
                        end
                        prop = prop_normalized;
                    end
                    value = img.(prop);
            end
        end
        
        function show(img, gamma_or_channel, gamma)
            figure;
            if nargin == 1
                imshow(img.get('CData'));
            else
                if nargin == 2
                    if ischar(gamma_or_channel)
                        channel = gamma_or_channel;
                        gamma = Image.defaultValue('gamma');
                    else
                        gamma = gamma_or_channel;
                        channel = '';
                    end
                end
                
                switch channel
                case {'red','green','green1','green2','blue','bayer'}
                    data = img.get('CData', 'channel', channel);
                    imshow(data .^ (1/gamma));
                    
                case 'histogram'
                    red = img.get('CData','channel', 'red');
                    h_red = histogram(red(:));
                    hold on
                    green1 = img.get('CData','channel', 'green1');
                    green2 = img.get('CData','channel', 'green2');
                    h_green = histogram([green1(:) green2(:)]);
                    hold on
                    blue = img.get('CData','channel', 'blue');
                    h_blue = histogram(blue(:));
                    
                    max_val = max(max([red(:) green1(:) green2(:) blue(:)]));
                    min_val = min(min([red(:) green1(:) green2(:) blue(:)]));
                    binWidth = (max_val - min_val) / 100;
                    h_red.FaceColor = [1 0 0];
                    h_red.EdgeColor = [1 0 0];
                    h_red.Normalization = 'probability';
                    h_red.BinWidth = binWidth;
                    h_green.FaceColor = [0 1 0];
                    h_green.EdgeColor = [0 1 0];
                    h_green.Normalization = 'probability';
                    h_green.BinWidth = binWidth;
                    h_blue.FaceColor = [0 0 1];
                    h_blue.EdgeColor = [0 0 1];
                    h_blue.Normalization = 'probability';
                    h_blue.BinWidth = binWidth;
                    
                otherwise
                    imshow(img.get('CData') .^ (1/gamma))
                end
            end
        end
    end
    
    methods(Static)
        
        function val = defaultValue(var)
            switch var
            case 'extractFeatures_method'
                val = 'MSER';
                
            case 'getTransformation_transformType'
                % TODO: change to projective as it returns better results
                % but fix for that the problem in SuperResolution
                val = 'affine'; % 'similarity' | 'affine' | 'projective'

            case 'gamma'
                val = 2.2;
                
            otherwise
                error('var is not supported');
            end
        end
        
        function img = read(filename, sensorAlignment, filename_ext, NaN_mask)
            if nargin <= 2
                filename_ext = ...
                    lower(regexp(filename,'[^\.]+$','match','lineanchors'));
                filename_ext = filename_ext{1};
            end
            
            switch filename_ext
                case {'raw','rw2','arw','dng'}
                    fprintf('Image.read(''%s'', ''%s'')\n', filename, filename_ext);
                    
                    global LIB_PATH
                    dcraw = [LIB_PATH,'/dcraw/dcraw.exe'];
                    
                    % load raw image data
                    % https://www.cybercom.net/~dcoffin/dcraw/dcraw.1.html
                    system([dcraw,' -4 -v -D -T ',pwd,'/',filename]);
                    filename_tiff = regexprep(filename,'[^\.]+$','tiff');
                    % t = Tiff(filename_tiff,'r');
                    data_bayer = imread(filename_tiff);
                    delete(filename_tiff);
                    
                    % get metadata
                    cmd = [dcraw,' -i -v ', filename];
                    [~,cmdout] = system(cmd);
                    metadata_lines = strsplit(strtrim(cmdout),'\n');
                    metadata_fields = regexprep(cellfun(@(line) strtrim(regexp(line,'^([^:]+)','match')),metadata_lines),'\s','_');
                    metadata_values = cellfun(@(line) strtrim(regexprep(line,'^[^:]+:','')),metadata_lines,'UniformOutput',false);
                    metadata = cell2struct(metadata_values.', metadata_fields);
                    
                    % correct metadata
                    metadata.Filter_pattern = regexprep(lower(metadata.Filter_pattern),'[^rgb]','');
                    %metadata.Camera_multipliers = str2double(strsplit(metadata.Camera_multipliers));
                    metadata.Daylight_multipliers = str2double(strsplit(metadata.Daylight_multipliers));
                    
                    % setup image object
                    img = Image();
                    img.CDataBayer = im2double(data_bayer);
                    if nargin >= 4
                        img.CDataBayer(NaN_mask) = NaN;
                    end
                    img.metaInfo = metadata;
                    img.sensorAlignment = img.metaInfo.Filter_pattern;
                    img.height = size(img.CDataBayer,1);
                    img.width  = size(img.CDataBayer,2);
                    
                    % whitebalance
                    %wb_multipliers = (metadata.Camera_multipliers(1:3) / metadata.Camera_multipliers(2)) .* metadata.Daylight_multipliers;
                    %image.whitebalanceBayer(wb_multipliers);
                    
                    % stretch image to have maximum on screen
                    % TODO: better check each channel when maximum of
                    % channel is reached
%                     image.CDataBayer = image.CDataBayer / max(image.CDataBayer(:));
                    %figure;
                    %imshow(image.CDataBayer);
                    
%                     image.demosaic();
%                     image.show();
                    
                    img.separateBayerChannels();
                    
                case {'jpg','png','tif'}
                    if nargin == 1
                        fprintf('Image.read(''%s'')\n', filename);
                    else
                        fprintf('Image.read(''%s'',''%s'')\n', filename, sensorAlignment);
                    end
                    
                    img = Image(im2double(imread(filename)));
                    
                    img.height = size(img.CData,1);
                    img.width  = size(img.CData,2);
                    
                    if nargin == 1
                        img.CDataRed = img.CData(:,:,1);
                        img.CDataGreen = img.CData(:,:,2);
                        img.CDataBlue = img.CData(:,:,3);
                        
                        if nargin >= 4
                            img.CDataRed = img.CDataRed .* mask;
                            img.CDataGreen = img.CDataGreen .* mask;
                            img.CDataBlue = img.CDataBlue .* mask;
                        end
                    else
                        % sensorAlignment set
                        img.transformRGBtoBayer(sensorAlignment);
                        img.separateBayerChannels();
                    end
                    
                otherwise
                    error('not supported file type');
            end
            
            %image.show('histogram');
            
        end
        
        function img = readDNG(filename, sensorAlignment)
        % https://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
            
            img = Image();
            
            if nargin == 1
                img.sensorAlignment = Image.defaultValue('sensorAlignment');
            else
                img.sensorAlignment = sensorAlignment;
            end

            % uncompressed DNG
            % from https://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
            warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
            t = Tiff(filename,'r');
            offsets = getTag(t,'SubIFD');
            setSubDirectory(t,offsets(1));
            % Create variable ’raw’, the Bayer CFA data
            img.CDataBayer = read(t);
            close(t);
            img.metaInfo = imfinfo(filename);
            warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');

            % Crop to only valid pixels
            % from https://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
            x_origin     = img.metaInfo.SubIFDs{1}.ActiveArea(2)+1; % +1 due to MATLAB indexing
            img.width  = uint16(img.metaInfo.SubIFDs{1}.DefaultCropSize(1));
            y_origin     = img.metaInfo.SubIFDs{1}.ActiveArea(1)+1;
            img.height = uint16(img.metaInfo.SubIFDs{1}.DefaultCropSize(2));
            img.CDataBayer = double(img.CDataBayer(...
                y_origin:y_origin+img.height-1,x_origin:x_origin+img.width-1));

            % Linearizing
            % from https://users.soe.ucsc.edu/~rcsumner/rawguide/RAWguide.pdf
            if isfield(img.metaInfo.SubIFDs{1},'LinearizationTable')
                ltab = img.metaInfo.SubIFDs{1}.LinearizationTable;
                img.CDataBayer  = ltab(img.CDataBayer+1);
            end
            black      = img.metaInfo.SubIFDs{1}.BlackLevel(1);
            saturation = img.metaInfo.SubIFDs{1}.WhiteLevel;
            img.CDataBayer = (img.CDataBayer-black)/(saturation-black);
            img.CDataBayer = max(0,min(img.CDataBayer,1));
            
            % white balance bayer
            %wb_multipliers = (image.metaInfo.AsShotNeutral).^-1;
            %wb_multipliers = wb_multipliers/wb_multipliers(2);
            
            img.separateBayerChannels();
            
            %image.demosaic();
            
            %image.brigthnessCorrection();
            
            % Gamma Correction
            %image.CData = image.CData .^ (1/2.2);
        end
        
        function transformation = getTransformation(img1, img2, transformType)
        % TODO: implement for all channels
        
            if nargin == 2
                transformType = Image.defaultValue('getTransformation_transformType');
            end
        
            map = containers.Map();
            if isempty(img1.CDataBayer) && isempty(img2.CDataBayer)
                map('red')    = {'featuresRed',    'validPtsRed'};
                map('green')  = {'featuresGreen',  'validPtsGreen'};
                map('blue')   = {'featuresBlue',   'validPtsBlue'};
            else
                map('red')    = {'featuresRed',    'validPtsRed'};
                map('green1') = {'featuresGreen1', 'validPtsGreen1'};
                map('green2') = {'featuresGreen2', 'validPtsGreen2'};
                map('blue')   = {'featuresBlue',   'validPtsBlue'};
            end
            
            matches_img1 = single([]);
            matches_img2 = single([]);
            
            tic
            for prop = map.keys
                map_value = map(prop{1});
                prop_features = map_value{1};
                prop_validPts = map_value{2};
                
                index_pairs = matchFeatures(...
                    img1.(prop_features),...
                    img2.(prop_features),...
                    'MaxRatio' , 0.7);

                matches_img1_channel = img1.(prop_validPts).Location(index_pairs(:,1),:);
                matches_img2_channel = img2.(prop_validPts).Location(index_pairs(:,2),:);
                
%                 figure;
%                 showMatchedFeatures(image1.CDataRed,image2.CDataRed,matches_img1_channel,matches_img1_channel);
                
                matches_img1 = cat(1, matches_img1, matches_img1_channel);
                matches_img2 = cat(1, matches_img2, matches_img2_channel);
            end
            toc
            
            % ~0.02 sec
            % remove similar matches (as channels are joint)
            matches_joined = sortrows(cat(2, matches_img1, matches_img2));
            tolerance = 1;
            % repeat several times to remove all close matches
            % TODO: for sure there are better methods for it, but it needs
            % also with many features max 0.005s
            for i = 1:4
                idx = [true; (sum(diff(matches_joined) > tolerance,2) > 0)];
                matches_joined = matches_joined(idx,:);
            end
            matches_img1 = matches_joined(:,1:2);
            matches_img2 = matches_joined(:,3:4);
            
            % 0.05 - 0.08 sec
            [transformation,inliers_img1,inliers_img2] = ...
                estimateGeometricTransform(...
                    matches_img2, matches_img1,...
                    transformType,...
                    'Confidence'   , 99.5,...
                    'MaxNumTrials' , 3000,...
                    'MaxDistance'  , 1.5);
            
%             figure;
%             showMatchedFeatures(...
%                 image1.get('CData','channel','red','normalized',true),...
%                 image1.get('CData','channel','red','normalized',true),...
%                 inliers_img1,inliers_img2);

            
            % scatteredInterpolant
        end
        
        function [CData1, CData2, correlation_curve] = adaptRelative(CData1, CData2)
            disp('adaptRelative(CData1, CData2)');
            
            % low pass pictures to find better correlation
%             H = fspecial('disk',5);
%             H = fspecial('gaussian',[13 13]);
%             CData1 = imfilter(CData1,H,'replicate');
%             CData2 = imfilter(CData2,H,'replicate');
            
            % filter out not finite data
            idx_finite = isfinite(CData1) & isfinite(CData2);
            
            samples = [CData1(idx_finite) CData2(idx_finite)];
            
            [overexposure_CData1, overexposure_CData2, p_1_to_2, correlation_curve] = ...
                Image.findOverexposure(samples(:,1),samples(:,2));
            
            if isfinite(overexposure_CData1)
                CData1(CData1 > overexposure_CData1) = NaN;
            end
            if isfinite(overexposure_CData2)
                CData2(CData2 > overexposure_CData2) = NaN;
            end
            
            CData2 = ( CData2 - p_1_to_2(2) ) / p_1_to_2(1);
            
%             if isfinite(overexp1)
%                 x = (samples(:,1) - start_overexposure);
%                 y = (samples(:,2) - mean_start_overexposure);
%                 a = -start_overexposure / mean_start_overexposure;
%                 selection_samples = y < a * x;
%                 selected_samples = samples(selection_samples,:);
%             else
%                 selected_samples = samples;
%             end

            global SAVE_FIGURES;
            global SHOW_FIGURES;
            if SAVE_FIGURES || SHOW_FIGURES
                if isfinite(overexposure_CData1) || isfinite(overexposure_CData2)
                    if SHOW_FIGURES
                        fig = figure;
                    else
                        fig = figure('Visible','Off');
                    end
                    
                    if isfinite(overexposure_CData1)
                        overexposed = CData1;
                    elseif isfinite(overexposure_CData2)
                        overexposed = CData2;
                    end
                    overexposed = image_mark_NaN(normalize_img(overexposed) ...
                        .^ (1/Image.defaultValue('gamma')));
                    image(overexposed);
                    
                    if SAVE_FIGURES
                        saveFigure(fig,'Image_adaptRelative_overexposed_areas');
                    end
                    if ~SHOW_FIGURES
                        close(fig);
                    end
                end
                
%                 % save compared values
%                 if SHOW_FIGURES
%                     fig = figure;
%                 else
%                     fig = figure('Visible','Off');
%                 end
%                 hold on;
%                 legends = {};
%                 
%                 % selected samples
%                 legends{end+1} = 'samples used for linear fit';
%                 plot(selected_samples(:,1), selected_samples(:,2), ...
%                     'LineStyle','none', ...
%                     'Marker','.', ...
%                     'MarkerSize',1, ...
%                     'MarkerEdgeColor',[0.2 0.2 1]);
%                 
%                 % not used samples
%                 if exist('selection_samples','var')
%                     legends{end+1} = 'not used samples';
%                     not_used_samples = samples(~selection_samples,:);
%                     plot(not_used_samples(:,1), not_used_samples(:,2), ...
%                         'LineStyle','none', ...
%                         'Marker','.', ...
%                         'MarkerSize',1, ...
%                         'MarkerEdgeColor',[0.7 0.7 1]);
%                 end
%                 
%                 % linear fit (selected samples)
%                 legends{end+1} = 'linear fit';
%                 plot([min_CData1 max_CData1], ...
%                     polyval(p_1_to_2,[min_CData1 max_CData1]), ...
%                     'Marker','none', ...
%                     'Color','black');
%                 
%                 if isfinite(overexposure_CData1)
%                     legends{end+1} = 'image 1 overexposure';
%                     plot([overexposure_CData1 overexposure_CData1],...
%                         [0.25*max_CData2 max_CData2], ...
%                         'Color','red', ...
%                         'Marker','none');
%                     legends{end+1} = 'start image 1 overexposure';
%                     plot([start_overexposure_CData1 start_overexposure_CData1],...
%                         [0.25*max_CData2 max_CData2], ...
%                         'LineStyle','--', ...
%                         'Color','red', ...
%                         'Marker','none');
%                     legends{end+1} = 'start image 1 overexposure mean';
%                     plot(start_overexposure_CData1,...
%                         mean_start_overexposure_CData1, ...
%                         'LineStyle','none', ...
%                         'Marker','.', ...
%                         'MarkerSize',10, ...
%                         'MarkerEdgeColor','red');
%                 end
%                 if isfinite(overexposure_CData2)
%                     legends{end+1} = 'image 2 overexposure';
%                     plot([0.25*max_CData1 max_CData1],...
%                         [overexposure_CData2 overexposure_CData2], ...
%                         'Color','red', ...
%                         'Marker','none');
%                     legends{end+1} = 'start image 2 overexposure';
%                     plot([0.25*max_CData1 max_CData1],...
%                         [start_overexposure_CData2 start_overexposure_CData2],...
%                         'LineStyle','--', ...
%                         'Color','red', ...
%                         'Marker','none');
%                     legends{end+1} = 'start image 2 overexposure mean';
%                     plot(mean_start_overexposure_CData2,...
%                         start_overexposure_CData2, ...
%                         'LineStyle','none', ...
%                         'Marker','.', ...
%                         'MarkerSize',10, ...
%                         'MarkerEdgeColor','red');
%                 end
%                 
%                 title('Brigthness of corresponding pixels');
%                 legend(legends,'Location','southeast');
%                 legend('boxoff');
%                 xlabel('brightness image 1');
%                 ylabel('brightness image 2');
%                 hold off;
%                 
%                 if SAVE_FIGURES
%                     saveFigure(fig,'Image_adaptRelative_usedValues',[2 1]);
%                 end
%                 if ~SHOW_FIGURES
%                     close(fig);
%                 end

                % TODO: set in CData1_filter all values that are inside CData2

                % save histogram
%                 if SHOW_FIGURES
%                     fig = figure;
%                 else
%                     fig = figure('Visible','Off');
%                 end
%                 histogram(CData1(:),edges);
%                 hold on;
%                 histogram(CData2(:),edges);
%                 hold on;
%                 histogram(CData2_adapted(:),edges);
%                 legend('reference values','values to adapt','adapted values');
%                 legend('boxoff');
%                 if SAVE_FIGURES
%                     saveFigure(fig,'Image_adaptRelative_histogram', [2 1]);
%                 end
%                 if ~SHOW_FIGURES
%                     close(fig);
%                 end
            end
        end
        
        function [overexp1,overexp2,p,diag_avg_samples] = findOverexposure(samples1,samples2)
            % find overexposure by comparing histogram of samples1 with values
            % of samples2
            % @param percentOverExposureArea - starting from which percentage
            % of the max value it's searched for the over exposure
            %
            % samples2 = p[2] * samples1 + p[1]
            
            disp('Image.findOverexposure(samples1,samples2)');
            
            NUM_BINS = 300;
            % tested by hand with which value the outlayers were removed
            % nice
            MIN_SAMPLES_PER_BIN = (length(samples1)/NUM_BINS)^0.5 * 0.15;
            
            % distance from beginning and inflection point
            DIST = NUM_BINS * 0.05;
            
            tic
            [hist_samples, bin_center] = hist3([samples1 samples2],[NUM_BINS NUM_BINS]);
            
            filtered_hist_samples = hist_samples;
            filtered_hist_samples( hist_samples < MIN_SAMPLES_PER_BIN ) = 0;
            
            % find maxima for [0 1; 1 0] diagonals => good results
            [spdiags_samples,d] = spdiags(flip(filtered_hist_samples));
            [~, diag_max_i] = max(spdiags_samples);
            diag_max_samples = [ ...
                bin_center{1}(diag_max_i)' ...
                bin_center{2}(d' - diag_max_i + NUM_BINS + 1)' ...
            ];
            diag_avg_i = ((1:size(spdiags_samples,1)) * spdiags_samples) ./ sum(spdiags_samples); 
            diag_avg_samples = [ ...
                interp1(1:NUM_BINS, bin_center{1}, diag_avg_i)' ...
                interp1(1:NUM_BINS, bin_center{2}, d' - diag_avg_i + NUM_BINS + 1)' ...
            ];
            
            % turn the found diagonal max values so the first and the last
            % value have the same height 
            diag_avg_samples_transp(:,2) = diag_avg_samples(:,2) - diag_avg_samples(1,2);
            diag_avg_samples_transp(:,1) = diag_avg_samples(:,1) - diag_avg_samples(1,1);
            dir = normr(diag_avg_samples_transp(end,:));
            rot = [dir(1) -dir(2); dir(2) dir(1)];
            diag_avg_samples_transp_rot = abs(diag_avg_samples_transp * rot);
            
            % as the maximum curve has after the overexposure nearly no
            % errors anymore:
            % search from the back for the first locale maximum
            inflection_idx = size(diag_avg_samples_transp_rot,1);
            found_max = 0;
            while found_max <= diag_avg_samples_transp_rot(inflection_idx, 2)
                found_max = diag_avg_samples_transp_rot(inflection_idx, 2);
                inflection_idx = inflection_idx - 1;
            end
            inflection_idx = inflection_idx+1;
            [found_max2, inflection_idx2] = max(diag_avg_samples_transp_rot(:,2));
            
            % sometimes when the overexposure is not linear the above
            % algorithm doesn't find the inflection - in that case use the
            % maximum
            if 1.1 * found_max < found_max2
                inflection_idx = inflection_idx2;
            end
            possible_inflection = diag_avg_samples(inflection_idx,:);
            
            % init
            angle_front_rear = 0;
            
            % for the tests and the calculations of the exact overexposure
            % point use values with a minimum distance to the found
            % inflection of the maximum line as the inflection has
            % sometimes some errors
            v_rear = diag_avg_samples(min([inflection_idx+DIST,...
                length(diag_avg_samples)]):end,:);
            if size(v_rear,1) >= 2
                v_front = diag_avg_samples(1+DIST:max([1,inflection_idx-DIST]),:);

                % calculate the dominant direction of each part of the curve
                % with the help of the covariance matrix
                cov_front_rear = [cov(v_front); cov(v_rear)];
                norm_cov_front_rear = sqrt(sum(cov_front_rear .^ 2,2));
                if norm_cov_front_rear(1) > norm_cov_front_rear(2)
                    dir_front_rear(1,1:2) = cov_front_rear(1,1:2);
                else
                    dir_front_rear(1,1:2) = cov_front_rear(2,1:2);
                end
                if norm_cov_front_rear(3) > norm_cov_front_rear(4)
                    dir_front_rear(2,1:2) = cov_front_rear(3,1:2);
                else
                    dir_front_rear(2,1:2) = cov_front_rear(4,1:2);
                end
                
                % scale values so angle can be computed also for samples
                % from very different ranges
                scaled_dir_front_rear(:,1) = dir_front_rear(:,1) / possible_inflection(1);
                scaled_dir_front_rear(:,2) = dir_front_rear(:,2) / possible_inflection(2);
                scaled_dir_front_rear = normr(scaled_dir_front_rear);

                % assume that an overexposure happened when the angle between
                % the front and the rear part is bigger then 5°
                angle_front_rear = acosd(scaled_dir_front_rear(1,:) * scaled_dir_front_rear(2,:)');
            end
            
            if angle_front_rear > 5 % degrees
                % calculate intersection of lines through v_front and
                % v_rear
                mean_v_front = mean(v_front);
                mean_v_rear = mean(v_rear);
                
                % solve:
                %   mean_v_front + t_1 * dir_front
                % = mean_v_rear + t_2 * dir_rear
                t = linsolve([dir_front_rear(1,:); -dir_front_rear(2,:)]', (mean_v_rear - mean_v_front)');
                
                inflection = t(1) * dir_front_rear(1,:) + mean_v_front;
                
                % cross product with [1,1,0]
                % > 0 => overexposure of samples1
                if dir_front_rear(2,1) - dir_front_rear(2,2) < 0
                    overexp1 = inflection(1);
                    overexp2 = NaN;
                % < 0 => overexposure of samples2
                else
                    overexp1 = NaN;
                    overexp2 = inflection(2);
                end
                % solve: y = p_1 * x + p_2 (x, y from mean_v_front)
                p_1 = dir_front_rear(1,2) / dir_front_rear(1,1);
                p_2 = mean_v_front(2) - p_1 * mean_v_front(1);
                p = [p_1, p_2];
            else
                overexp1 = NaN;
                overexp2 = NaN;
                % TODO: calculate out of cov & mean
                p = polyfit(diag_avg_samples(1+DIST:end-DIST,1),diag_avg_samples(1+DIST:end-DIST,2),1);
            end

            toc
            
            global SAVE_FIGURES;
            global SHOW_FIGURES;
            if SAVE_FIGURES || SHOW_FIGURES
                if SHOW_FIGURES
                    fig = figure;
                else
                    fig = figure('Visible','Off');
                end

                hold on
                legends = {};
                x = [bin_center{1}(1) bin_center{1}(end)];
                y = [bin_center{2}(1) bin_center{2}(end)];
                image(x, y, hist_samples.^(1/4),'CDataMapping','scaled',...
                    'AlphaData',hist_samples > 0)
                colormap(fig,flipud(colormap('bone')))
                
                plot(diag_max_samples(:,1),diag_max_samples(:,2),'r--')
                legends{end+1} = 'diagonal max';
                
                plot(diag_avg_samples(:,1),diag_avg_samples(:,2),'g--')
                legends{end+1} = 'diagonal average';
                
                plot(possible_inflection(1),possible_inflection(2),'bo')
                legends{end+1} = 'possible inflection of line';
                
                if isfinite(overexp1) || isfinite(overexp2)
                    plot(inflection(1),inflection(2),'bx')
                    legends{end+1} = 'overexposure point';
                    
                    x_front = [x(1); 1.02*inflection(1)];
                else
                    x_front = x;
                    
                end
                plot(x_front,polyval(p,x_front),...
                    'LineStyle','--','Color','black','Marker','none');
                legends{end+1} = 'linear fit';
                
                if isfinite(overexp1) || isfinite(overexp2)
                    % solve: y = p_1 * x + p_2 (x, y from mean_v_front)
                    p_1 = dir_front_rear(2,2) / dir_front_rear(2,1);
                    p_2 = mean_v_rear(2) - p_1 * mean_v_rear(1);
                    p_rear = [p_1, p_2];
                    
                    x_rear = [0.95*inflection(1); x(end)];
                    %plot(x_rear,polyval(p_rear,x_rear),...
                    %    'LineStyle','--','Color','black','Marker','none');
                end
                
                xlabel('values samples 1')
                ylabel('values samples 2')
                
                legend(legends,'Location','southeast');
                legend('boxoff');
                hold off
                
                if SAVE_FIGURES
                    saveFigure(fig,'Image_findOverexposure',[2 2]);
                end
                if ~SHOW_FIGURES
                    close(fig);
                end
            end
            
        end
        
        function [overexp1,overexp2,p] = findOverexposureOld3(samples1,samples2)
            % find overexposure by comparing histogram of samples1 with values
            % of samples2
            % @param percentOverExposureArea - starting from which percentage
            % of the max value it's searched for the over exposure
            %
            % samples2 = p[2] * samples1 + p[1]
            
            disp('Image.findOverexposure(samples1,samples2)');
            
            NUM_BINS = 300;
            MIN_SAMPLES_PER_BIN = max(20,(length(samples1)/NUM_BINS)^0.5 * 0.3);
            
            tic
            [hist_samples, bin_center] = hist3([samples1 samples2],[NUM_BINS NUM_BINS]);
            
            filtered_hist_samples = hist_samples;
            filtered_hist_samples( hist_samples < MIN_SAMPLES_PER_BIN ) = 0;
            
            % find maxima in x and y direction => not good result for
            % overexposured parts
%             [max12, idx12] = max(hist_samples);
%             idx_max12 = max12 > MIN_SAMPLES_PER_BIN;
%             [max21, idx21] = max(hist_samples,[],2);
%             idx_max21 = max21 > MIN_SAMPLES_PER_BIN;
%             max_samples12 = [bin_center{1}(idx_max12)' bin_center{2}(idx12(idx_max12))'];
%             max_samples21 = [bin_center{1}(idx21(idx_max21))' bin_center{2}(idx_max21)'];
            
            % find maxima for [0 1; 1 0] diagonals => good results
            [spdiags_samples,d] = spdiags(flip(filtered_hist_samples));
            [~, idx_max_diag] = max(spdiags_samples);
            diag_max_samples = [ ...
                bin_center{1}(idx_max_diag)' ...
                bin_center{2}(d' - idx_max_diag + NUM_BINS + 1)' ...
            ];
            
            function idx = find_inflection(data, levels)
                %[c1, l] = wavedec(data(:,2),levels,filter);
                [data_wavelet(:,1), l] = wavedec(data(:,1),levels,'haar');
                data_wavelet(:,2) = wavedec(data(:,2),levels,'haar');
                
                x_data = 1:size(data,1);
                
                l(1) = 1;
                
                fig = figure;
                axes = 1:levels;
                for level = 1:levels
                    signal = diff(data_wavelet(l(level+1)+1:l(level+2),:));
                    signal = sum(signal(1:end-1,:) .* signal(2:end,:),2);
                    
                    if level == 1
                        [~, pos_max] = max(signal);
                    else
                        % calculate max position by checking the area
                        % defined by the previous max position
                        pos_max = 2*pos_max;
                        [~, pos] = max(signal(pos_max-1:pos_max+1));
                        pos_max = pos_max -2 + pos;
                    end
                    
                    axes(level) = subplot(levels,1,level);
                    step = 2^(levels-level);
                    x_level = x_data(step:step:end);
                    x_level = x_level(1:size(signal,1));
                    plot(x_level,signal,'b-');
                    hold on
                    plot(x_level(pos_max),signal(pos_max),'ro')
                end
                linkaxes(axes,'x')
                
                idx = 2*pos_max;
            end

            center_inflection_idx = find_inflection(diag_max_samples, 5);
        
            start_inflection_idx = center_inflection_idx - 1;
            end_inflection_idx = center_inflection_idx + 1;
            
            inflection_range_idx = ...
                [start_inflection_idx,...
                 center_inflection_idx,...
                 end_inflection_idx];
            inflection_range = diag_max_samples(inflection_range_idx,:);
            
            v_front = diag_max_samples(1:start_inflection_idx,:);
            v_rear = diag_max_samples(end_inflection_idx:end,:);
            
            dir_front_rear = normr([v_front(end,:) - v_front(1,:); v_rear(end,:) - v_rear(1,:)]);
            ang_front_rear = acosd(dir_front_rear(1,:) * dir_front_rear(2,:)');
            
            if ang_front_rear > 5 % degrees
                % cross product with [1,1,0]
                % > 0 => overexposure of samples1
                if dir_front_rear(2,1) - dir_front_rear(2,2) < 0
                    overexp1 = inflection_range(1,1);
                    overexp2 = NaN;
                % < 0 => overexposure of samples2
                else
                    overexp1 = NaN;
                    overexp2 = inflection_range(1,2);
                end
                p = polyfit(v_front(:,1),v_front(:,2),1);
            else
                overexp1 = NaN;
                overexp2 = NaN;
                p = polyfit(diag_max_samples(:,1),diag_max_samples(:,2),1);
            end
            
            % TODO: compare p and S for front and rear and decide if
            % overflow occured
            
%             R = qr(v_front); 
%             x = R\(R'\(v_front'*b));
%             r = b - A*x;
%             err = R\(R'\(A'*r));
%             x = x + err;

            toc
            
            global SAVE_FIGURES;
            global SHOW_FIGURES;
            if SAVE_FIGURES || SHOW_FIGURES
                if SHOW_FIGURES
                    fig = figure;
                else
                    fig = figure('Visible','Off');
                end

                hold on
                legends = {};
                x = [bin_center{1}(1) bin_center{1}(end)];
                y = [bin_center{2}(1) bin_center{2}(end)];
                image(x, y, hist_samples.^(1/4),'CDataMapping','scaled','AlphaData',hist_samples > 0)
                colormap(fig,flipud(colormap('bone')))
                %plot(max_samples12(:,1),max_samples12(:,2))
                %plot(max_samples21(:,1),max_samples21(:,2))
                %plot(diag_max_samples(:,1),diag_max_samples(:,2),'r')
                
                plot(diag_max_samples(:,1),diag_max_samples(:,2),'r--')
                legends{end+1} = 'diagonal max';
                
                if isfinite(overexp1) || isfinite(overexp2)
                    plot(inflection_range(:,1),inflection_range(:,2),'bx')
                    legends{end+1} = 'found inflection area';
                    
                    x_fit = [x(1); 1.02*inflection_range(3,1)];
                else
                    x_fit = x;
                end
                plot(x_fit,polyval(p,x_fit),'LineStyle','--','Color','black','Marker','none');
                legends{end+1} = 'linear fit max histogram';
                
                if isfinite(overexp1)
                    % TODO: plot overexposure
                end
                
                xlabel('values samples 1')
                ylabel('values samples 2')
                
                legend(legends,'Location','southeast');
                legend('boxoff');
                hold off
                
                if SAVE_FIGURES
                    saveFigure(fig,'Image_findOverexposure',[2 2]);
                end
                if ~SHOW_FIGURES
                    close(fig);
                end
            end
        end

        function [overexp1,overexp2,p] = findOverexposureOld2(samples1,samples2)
            % find overexposure by comparing histogram of samples1 with values
            % of samples2
            % @param percentOverExposureArea - starting from which percentage
            % of the max value it's searched for the over exposure
            %
            % samples2 = p[2] * samples1 + p[1]
            
            disp('Image.findOverexposure(samples1,samples2)');
            
            NUM_BINS = 300;
            MIN_SAMPLES_PER_BIN = max(20,(length(samples1)/NUM_BINS)^0.5 * 0.3);
            % odd kernel size important to not move signal
            KERNEL_SIZE = round(NUM_BINS/100) * 2 + 1;
            
            tic
            [hist_samples, bin_center] = hist3([samples1 samples2],[NUM_BINS NUM_BINS]);
            
            filtered_hist_samples = hist_samples;
            filtered_hist_samples( hist_samples < MIN_SAMPLES_PER_BIN ) = 0;
            
            % find maxima in x and y direction => not good result for
            % overexposured parts
%             [max12, idx12] = max(hist_samples);
%             idx_max12 = max12 > MIN_SAMPLES_PER_BIN;
%             [max21, idx21] = max(hist_samples,[],2);
%             idx_max21 = max21 > MIN_SAMPLES_PER_BIN;
%             max_samples12 = [bin_center{1}(idx_max12)' bin_center{2}(idx12(idx_max12))'];
%             max_samples21 = [bin_center{1}(idx21(idx_max21))' bin_center{2}(idx_max21)'];
            
            % find maxima for [0 1; 1 0] diagonals => good results
            [spdiags_samples,d] = spdiags(flip(filtered_hist_samples));
            [~, idx_max_diag] = max(spdiags_samples);
            diag_max_samples = [ ...
                bin_center{1}(idx_max_diag)' ...
                bin_center{2}(d' - idx_max_diag + NUM_BINS + 1)' ...
            ];
        
            % smooth max samples to get rid of stairs structure
            kernel = fspecial('average',[KERNEL_SIZE 1]);
            conv_diag_max_samples = conv2(diag_max_samples,kernel,'valid');            
            
            diff_diag_max_samples = diff(conv_diag_max_samples);
            % cos(alpha) * 2^0.5/2 - 1 = [1 1] * diff_diag_max_samples(i,:) - 1
            % inflection_signal = sum(diff_diag_max_samples,2);
            
            % vector product between current and following difference
            inflection_signal = sum(diff_diag_max_samples(1:end-1,:) .* diff_diag_max_samples(2:end,:),2);
            
            dist_injection_samples = KERNEL_SIZE;
            injection_signal_samples = inflection_signal(1:dist_injection_samples:end);
            [~,center_inflection_idx] = max(diff(injection_signal_samples));
            center_inflection_idx = (center_inflection_idx - 1) * dist_injection_samples + 1 + KERNEL_SIZE;
            
            start_inflection_idx = center_inflection_idx - KERNEL_SIZE;
            end_inflection_idx = center_inflection_idx + KERNEL_SIZE;
            % in case the inflection is nearly at the end of the samples
            if end_inflection_idx > size(conv_diag_max_samples,1)
                end_inflection_idx = center_inflection_idx;
            end
            
            inflection_range_idx = ...
                [start_inflection_idx,...
                 center_inflection_idx,...
                 end_inflection_idx];
            inflection_range = conv_diag_max_samples(inflection_range_idx,:);
            
            % normalize x and y of diag_max_samples to [0 1] range
            norm_diag_max_samples(:,1) = mat2gray(diag_max_samples(:,1));
            norm_diag_max_samples(:,2) = mat2gray(diag_max_samples(:,2));
            
            v_front = norm_diag_max_samples(1:start_inflection_idx,:);
            v_rear = norm_diag_max_samples(end_inflection_idx:end,:);
            
            dir_front_rear = normr([v_front(end,:) - v_front(1,:); v_rear(end,:) - v_rear(1,:)]);
            ang_front_rear = acosd(dir_front_rear(1,:) * dir_front_rear(2,:)');
            
            if ang_front_rear > 5 % degrees
                % cross product with [1,1,0]
                % > 0 => overexposure of samples1
                if dir_front_rear(2,1) - dir_front_rear(2,2) < 0
                    overexp1 = inflection_range(1,1);
                    overexp2 = NaN;
                % < 0 => overexposure of samples2
                else
                    overexp1 = NaN;
                    overexp2 = inflection_range(1,2);
                end
                p = polyfit(v_front(:,1),v_front(:,2),1);
            else
                overexp1 = NaN;
                overexp2 = NaN;
                p = polyfit(diag_max_samples(:,1),diag_max_samples(:,2),1);
            end
            
            % TODO: compare p and S for front and rear and decide if
            % overflow occured
            
%             R = qr(v_front); 
%             x = R\(R'\(v_front'*b));
%             r = b - A*x;
%             err = R\(R'\(A'*r));
%             x = x + err;

            toc
            
            global SAVE_FIGURES;
            global SHOW_FIGURES;
            if SAVE_FIGURES || SHOW_FIGURES
                if SHOW_FIGURES
                    fig = figure;
                else
                    fig = figure('Visible','Off');
                end

                ax11 = subplot(2,2,1);
                hold on
                legends = {};
                x = [bin_center{1}(1) bin_center{1}(end)];
                y = [bin_center{2}(1) bin_center{2}(end)];
                image(x, y, hist_samples.^(1/4),'CDataMapping','scaled','AlphaData',hist_samples > 0)
                colormap(fig,flipud(colormap('bone')))
                %plot(max_samples12(:,1),max_samples12(:,2))
                %plot(max_samples21(:,1),max_samples21(:,2))
                %plot(diag_max_samples(:,1),diag_max_samples(:,2),'r')
                
                plot(diag_max_samples(:,1),diag_max_samples(:,2),'r--')
                legends{end+1} = 'diagonal max';
                
                plot(conv_diag_max_samples(:,1),conv_diag_max_samples(:,2),'r-')
                legends{end+1} = 'diagonal max (smoothed)';
                
                if isfinite(overexp1) || isfinite(overexp2)
                    plot(inflection_range(:,1),inflection_range(:,2),'bx')
                    legends{end+1} = 'found inflection area';
                    
                    x_fit = [x(1); 1.02*inflection_range(3,1)];
                else
                    x_fit = x;
                end
                plot(x_fit,polyval(p,x_fit),'LineStyle','--','Color','black','Marker','none');
                legends{end+1} = 'linear fit max histogram';
                
                if isfinite(overexp1)
                    % TODO: plot overexposure
                end
                
                xlabel('values samples 1')
                ylabel('values samples 2')
                
                legend(legends,'Location','southeast');
                legend('boxoff');
                hold off
                
                ax12 = subplot(2,2,2);
                hold on
                plot(inflection_signal,conv_diag_max_samples(1:end-2,2),'Color',[1 0.5 0.5]);
                plot(injection_signal_samples(1:end),conv_diag_max_samples(1:dist_injection_samples:end-2,2),'ro')
                plot(injection_signal_samples(round(inflection_range_idx/KERNEL_SIZE + 1)),conv_diag_max_samples(inflection_range_idx,2),'bx')
                xlabel('diff filtered values samples 1')
                ylabel('values samples 2')
                hold off

                ax21 = subplot(2,2,3);
                hold on
                plot(conv_diag_max_samples(1:end-2,1),inflection_signal,'Color',[1 0.5 0.5])
                plot(conv_diag_max_samples(1:dist_injection_samples:end-2,1),injection_signal_samples,'ro')
                plot(conv_diag_max_samples(inflection_range_idx,1),injection_signal_samples(round(inflection_range_idx/KERNEL_SIZE + 1)),'bx')
                xlabel('values samples 1')
                ylabel('diff filtered values samples 2')
                hold off

                linkaxes([ax11,ax21],'x')
                linkaxes([ax11,ax12],'y')
                
                if SAVE_FIGURES
                    saveFigure(fig,'Image_findOverexposure',[2 2]);
                end
                if ~SHOW_FIGURES
                    close(fig);
                end
            end
        end
        
        function [overexp, start_overexp, mean_start_overexp, std_start_overexp] = findOverexposureOld(samples1,samples2)
            % find overexposure by comparing histogram of samples1 with values
            % of samples2
            % @param percentOverExposureArea - starting from which percentage
            % of the max value it's searched for the over exposure
            
            disp('Image.findOverexposure(samples1,samples2)');
            tic
            
            NUM_BINS = 150;
            MIN_OVEREXPOSURE_PROP = 15;
            MIN_CHECKED_BIN = 0.7 * NUM_BINS;
            OVEREXP_START_FACTOR = 1.6;
            
            function [overexp, start_overexp, mean_start_overexp, bin_edges, bin_mean] = compare_samples(samples1, samples2)
                MIN_SAMPLES_PER_BIN = max(20,length(samples1)/NUM_BINS*0.001);

                % numBins = size(bin_edges,2) - 1;
                bin_edges = linspace(0, max(samples1(:)), NUM_BINS+1);

                % 0.1sec / 10mpix
                [h,~,whichBin] = histcounts(samples1, bin_edges);

                % NUM_BINS = 150: 1.8sec / 10mpix
                bin_mean(1:NUM_BINS) = NaN;
                for i = MIN_CHECKED_BIN:NUM_BINS
                    if h(i) >= MIN_SAMPLES_PER_BIN
                        bin_s2 = samples2((whichBin == i));
                        bin_mean(i) = mean(bin_s2);
                    end
                end

                % 0.005sec / 10mpix
                diff_bin_mean = diff(bin_mean);
                mean_diff_bin_mean = mean(diff_bin_mean,'omitnan');
                max_diff_bin_mean = nanmax(diff_bin_mean);

                overexp_prop = max_diff_bin_mean / mean_diff_bin_mean;
                if overexp_prop >= MIN_OVEREXPOSURE_PROP
                    [~,bin_overexp] = max(diff_bin_mean == max_diff_bin_mean);
                    overexp = mean(bin_edges(bin_overexp+1:bin_overexp+2));
                    bin_start_overexp = max([MIN_CHECKED_BIN round(NUM_BINS - OVEREXP_START_FACTOR * (NUM_BINS - bin_overexp))]);
                    start_overexp = mean(bin_edges(bin_start_overexp+1:bin_start_overexp+2));
                    mean_start_overexp = bin_mean(bin_start_overexp);
                    %std_start_overexp = std(samples2((whichBin == bin_start_overexp)));
                else
                    overexp = NaN;
                    start_overexp = NaN;
                    mean_start_overexp = NaN;
                    %std_start_overexp = NaN;
                end
            end
            
            [overexp_s1, start_overexp_s1, mean_start_overexp_s1, bin_edges_s1, bin_mean_s1] = ...
                compare_samples(samples1,samples2);
            [overexp_s2, start_overexp_s2, mean_start_overexp_s2, bin_edges_s2, bin_mean_s2] = ...
                compare_samples(samples2,samples1);
            
            if isfinite(overexp_s1)
                overexp = overexp_s1;
                start_overexp = start_overexp_s1;
                mean_start_overexp = mean_start_overexp_s1;
            elseif isfinite(overexp_s2)
                overexp = overexp_s2;
                start_overexp = mean_start_overexp_s2;
                mean_start_overexp = start_overexp_s2;
            end
            
            global SAVE_FIGURES;
            global SHOW_FIGURES;
            if SAVE_FIGURES || SHOW_FIGURES
                if SHOW_FIGURES
                    fig = figure;
                else
                    fig = figure('Visible','Off');
                end
                
                bin_center_s1 = (bin_edges_s1(1:end-1) + bin_edges_s1(2:end)) / 2;
                bin_center_s2 = (bin_edges_s2(1:end-1) + bin_edges_s2(2:end)) / 2;
                axis_s1_s2 = [min([bin_edges_s1(MIN_CHECKED_BIN) bin_mean_s2]) bin_edges_s1(end) min([bin_edges_s2(MIN_CHECKED_BIN) bin_mean_s1]) bin_edges_s2(end)];
                
                plot_style_s1 = 'b';
                plot_style_s2 = 'r';
                
                subplot(2,2,1)
                hold on
                plot(samples1, samples2, 'LineStyle', 'none', 'Marker', '.', 'MarkerEdgeColor', [0.7 0.7 0.7], 'MarkerSize', 1)
                plot(bin_center_s1, bin_mean_s1, plot_style_s1, ...
                    bin_mean_s2, bin_center_s2, plot_style_s2)
                hold off
                axis(axis_s1_s2)
                legend({'samples', 'mean samples 1','mean samples 2'},'Location','northwest')
                legend('boxoff');
                title('mean samples')
                xlabel('values samples 1')
                ylabel('values samples 2')
                
                diff_bin_mean_s1 = diff(bin_mean_s1);
                diff_bin_mean_s2 = diff(bin_mean_s2);
                axis_diff = [min([diff_bin_mean_s1 diff_bin_mean_s2]) ...
                    max([diff_bin_mean_s1 diff_bin_mean_s2])];
                
                subplot(2,2,3)
                plot(bin_center_s1(2:end), diff_bin_mean_s1,plot_style_s1)
                axis([axis_s1_s2(1:2) axis_diff])
                title('diff mean samples 2')
                xlabel('values samples 1')
                ylabel('diff mean samples 1')
                
                subplot(2,2,2);
                plot(diff_bin_mean_s2, bin_center_s2(2:end),plot_style_s2);
                axis([axis_diff axis_s1_s2(3:4)])
                title('diff mean samples 2')
                xlabel('values samples 1')
                ylabel('values samples 2')
                
                if SAVE_FIGURES
                    saveFigure(fig,'Image_findOverexposure',[2 2]);
                end
                if ~SHOW_FIGURES
                    close(fig);
                end
            end
            
            toc
        end
        
    end
    
end

