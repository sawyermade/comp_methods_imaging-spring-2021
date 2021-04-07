img_path_in = '~/Desktop/lenna_bw.png';
img_path_out = '~/Desktop/output_conv2d.png';
h = imread(img_path_in);
k = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1] / 320;

y = myconv2d(h, k);
imwrite(y, img_path_out);

function y = myconv2d(h, k)
    padding = (length(k) - 1) / 2;
    hdims = size(h);
    h = double(h);
    k = double(k);
    
    y = double(zeros(hdims));
    p = padding;
    for i = 1+p:hdims(1)-p
        for j = 1+p:hdims(2)-p
            hs = h(i-p:i+p, j-p:j+p);
            y(i,j) = sum(hs .* k, 'all');
        end
    end
    
    y = uint8(y);
end