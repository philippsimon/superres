%mexOpenCV detectORBFeaturesOCV.cpp
%mexOpenCV extractORBFeaturesOCV.cpp

bayer_mode = 'gbrg';

if exist('images', 'var') == 0
    images = cell(2, 1);
    
    images{1} = readDNG('../Superresolution/DSC00036.dng', bayer_mode);
    % img1 = im2uint8(img1);
    
    images{2} = readDNG('../Superresolution/DSC00037.dng', bayer_mode);
    % img2 = im2uint8(img2);
end

%%
% You can experiment by varying the scale and rotation of the input image.
% However, note that there is a limit to the amount you can vary the scale
% and theta before the feature detector fails to find enough features.

%% Step 3: Find Matching Features Between Images
% Detect features in both images.
% ptsImg1 = detectORBFeaturesOCV(img1);
% ptsImg2 = detectORBFeaturesOCV(img2);
ptsImg1 = detectSURFFeatures(rgb2gray(images{1}));
ptsImg2 = detectSURFFeatures(rgb2gray(images{2}));

%%
% Extract feature descriptors.
%[featuresImg1_uint8, validPtsImg1] = extractORBFeaturesOCV(img1,  ptsImg1);
%[featuresImg2_uint8, validPtsImg2] = extractORBFeaturesOCV(img2, ptsImg2);
%featuresImg1 = binaryFeatures(featuresImg1_uint8);
%featuresImg2 = binaryFeatures(featuresImg2_uint8);

[featuresImg1, validPtsImg1] = extractFeatures(rgb2gray(images{1}), ptsImg1);
[featuresImg2, validPtsImg2] = extractFeatures(rgb2gray(images{2}), ptsImg2);

%%
% Match features by using their descriptors.
indexPairs = matchFeatures(featuresImg1, featuresImg2);

%%
% Retrieve locations of corresponding points for each image.
matchedImg1 = validPtsImg1.Location(indexPairs(:,1),:);
matchedImg2 = validPtsImg2.Location(indexPairs(:,2),:);

%% Step 4: Estimate Transformation
% Find a transformation corresponding to the matching point pairs using the
% statistically robust M-estimator SAmple Consensus (MSAC) algorithm, which
% is a variant of the RANSAC algorithm. It removes outliers while computing
% the transformation matrix. You may see varying results of the
% transformation computation because of the random sampling employed by the
% MSAC algorithm.
[tform, inlierImg1, inlierImg2] = estimateGeometricTransform(...
    matchedImg2, matchedImg1, 'similarity');

%%
% Display matching point pairs used in the computation of the
% transformation.
%figure; 
%showMatchedFeatures(images{1},images{2},inlierImg1,inlierImg2);
%title('Matching points (inliers only)');
%legend('ptsOriginal','ptsDistorted');

%% Step 5: Solve for Scale and Angle
% Use the geometric transform, tform, to recover the scale and angle.
% Since we computed the transformation from the distorted to the original
% image, we need to compute its inverse to recover the distortion.
%
%  Let sc = s*cos(theta)
%  Let ss = s*sin(theta)
%
%  Then, Tinv = [sc -ss  0;
%                ss  sc  0;
%                tx  ty  1]
%
%  where tx and ty are x and y translations, respectively.
%

%%
% Compute the inverse transformation matrix.
Tinv  = tform.invert.T;

ss = Tinv(2,1);
sc = Tinv(1,1);
scaleRecovered = sqrt(ss*ss + sc*sc);
thetaRecovered = atan2(ss,sc)*180/pi;

%%
% The recovered values should match your scale and angle values selected in
% *Step 2: Resize and Rotate the Image*.

%% Step 6: Recover the Original Image
% Recover the original image by transforming the distorted image.
outputView = imref2d(size(images{1}));
img2Translated  = imwarp(images{2},tform,'OutputView',outputView);

C = imfuse(images{1},img2Translated,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
imshow(C)

%%
% Compare |recovered| to |original| by looking at them side-by-side in a
% montage.
% figure, imshowpair(images{1},img2Translated,'montage')

%%
% The |recovered| (right) image quality does not match the |original|
% (left) image because of the distortion and recovery process. In
% particular, the image shrinking causes loss of information. The artifacts
% around the edges are due to the limited accuracy of the transformation.
% If you were to detect more points in *Step 3: Find Matching Features
% Between Images*, the transformation would be more accurate. For example,
% we could have used a corner detector, detectFASTFeatures, to complement
% the SURF feature detector which finds blobs. Image content and image size
% also impact the number of detected features.
