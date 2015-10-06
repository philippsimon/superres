function [img,img_min,img_max] = normalize_img( img, img_min, img_max )
%NORMALIZE_IMG Normalize image to [0 1] range

    if nargin == 1
        img_min = min(img(:));
        img_max = max(img(:));
    end
    
    img = ( img - img_min ) / (img_max - img_min);

end