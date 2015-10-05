% test Image class with gbrg bayer image
CData =...
  [ 0.1 0.2 0.3 0.4;
    0.5 0.6 0.7 0.8;
    0.9 0.0 0.1 0.2;
    0.3 0.3 0.4 0.5 ];
test_image_bayer_gbrg = Image(CData, 'gbrg');
assert( isequal(test_image_bayer_gbrg.CDataBayer,CData),...
    'test_image_bayer_gbrg: CDataBayer wrong');
assert( isequal(test_image_bayer_gbrg.CDataRed,...
   [ 0.5 0.7;...
     0.3 0.4]), 'test_image_bayer_gbrg: CDataRed wrong');
assert( isequal(test_image_bayer_gbrg.CDataGreen1,...
   [ 0.1 0.3;
     0.9 0.1]), 'test_image_bayer_gbrg: CDataGreen1 wrong');
assert( isequal(test_image_bayer_gbrg.CDataGreen2,...
   [ 0.6 0.8;
     0.3 0.5]), 'test_image_bayer_gbrg: CDataGreen2 wrong');
assert( isequal(test_image_bayer_gbrg.CDataBlue,...
   [ 0.2 0.4;...
     0.0 0.2]), 'test_image_bayer_gbrg: CDataBlue wrong');
 
% test getLinearCorelation
CData1 = [  NaN  NaN  NaN  NaN  ;
            Inf  Inf  Inf  Inf  ;
           -Inf -Inf -Inf -Inf  ;
            1    1    1    1    ;
            1    2    3    4    ];
CData2 = [  NaN  Inf -Inf  1    ;
            NaN  Inf -Inf  1    ;
            NaN  Inf -Inf  1    ;
            NaN  Inf -Inf  3    ;
            3    5    7    9    ];
Image.getLinearCorrelation(CData1,CData2);
 
%test_image_dng = Image.read('test_data/test_text/P1080860.RW2');
%test_image_dng.whitebalanceBayer();
%test_image_dng.demosaic();
%test_image_dng.colorTransformation('sRGB');
%test_image_dng.whitebalance();
%test_image_dng.brigthnessCorrection();
%test_image_dng.show('red');
 
disp('all tests for class Image passed');