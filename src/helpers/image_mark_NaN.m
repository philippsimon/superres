% Copyright (c) 2015, Philipp Simon Schmidt
% For more details see LICENSE.txt and AUTHORS.txt

function img = image_mark_NaN(img, color)
    if nargin == 1
        color = [1 0 0];
    end

    % find overexposed/unused pixels
    img_idx = isnan(img);
    
    % copy values to all layers to get gray image
    img_r = img;
    img_g = img_r;
    img_b = img_r;

    % mark overexposed/unused pixels in red 
    img_r(img_idx) = color(1);
    img_g(img_idx) = color(2);
    img_b(img_idx) = color(3);

    img = cat(3,img_r,img_g,img_b);
end