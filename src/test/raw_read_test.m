row = 2592;
col = 1944;
fin = fopen('test_data/test.raw','r');
I = fread(fin,row*col,'uint16=>uint16'); 
Z = reshape(I,row,col);
Z = Z';
k = imshow(normalize_img(double(Z)));