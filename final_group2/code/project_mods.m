% Choose Measurements To Load

load data/data_s_u.mat;
% load data/data_exit_sign.mat;
% load data/data_mannequin.mat;
% load data/data_outdoor_s.mat;
reconstruction(rect_data, width)

function reconstruction(rect_data, width)

bin_res = 4e-12;
speed_light = 3e8;

backprop = 0;
diffuse = 0;
K = 2;
snr = 0.8;

z_trim = 600;


% load data/data_s_u.mat
z_offset = 800;

image_width = size(rect_data,1); %N
temporal_bins = size(rect_data,3); %M

hist_range = temporal_bins * speed_light * bin_res;


% 
for k = 1:K
    
    temporal_bins = temporal_bins / 2;
    
    bin_res = 2*bin_res;
    
    rect_data = rect_data(:,:,1:2:end) + rect_data(:,:,2:2:end);
    
    z_trim = round(z_trim/2);
    
    z_offset = round(z_offset/2);
    
end



rect_data(:,:,1:z_trim) = 0; 
    
psf = definePsf(image_width,temporal_bins,width./hist_range);

fpsf = fftn(psf);

invpsf = conj(fpsf) ./ (abs(fpsf).^2 + 1./snr);

[mtx, mtxi] = resampleOperator(temporal_bins);


data = permute(rect_data,[3 2 1]);

grid_z = repmat(linspace(0,1,temporal_bins)',[1 image_width image_width]);



if (diffuse)
    data = data.*(grid_z.^4);
else
    data = data.*(grid_z.^2);
end


tdata = zeros(2.*temporal_bins,2.*image_width,2.*image_width);
tdata(1:end./2,1:end./2,1:end./2)  = reshape(mtx*data(:,:),[temporal_bins image_width image_width]);

% tic;
% tdata = imgaussfilt(tdata,0.75,'filterDomain','frequency');
% toc;

tvol = ifftn(fftn(tdata).*invpsf);


tvol = tvol(1:end./2,1:end./2,1:end./2);

vol  = reshape(mtxi*tvol(:,:),[temporal_bins image_width image_width]);
vol  = max(real(vol),0);

tic_z = linspace(0,hist_range./2,size(vol,1));
tic_y = linspace(-width,width,size(vol,2));
tic_x = linspace(-width,width,size(vol,3));

ind = round(temporal_bins.*2.*width./(hist_range./2));
vol = vol(:,:,end:-1:1);
vol = vol((1:ind)+z_offset,:,:);


tic_z = tic_z((1:ind)+z_offset);

output_image = squeeze(max(vol,[],1));

figure;
imagesc(tic_x,tic_y, output_image);
title('Front Raw Output');
set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
xlabel('x (m)');
ylabel('y (m)');
colormap('gray');
axis square;

guassian_filtered_image = imgaussfilt(output_image,1);

figure;
imagesc(tic_x,tic_y,guassian_filtered_image );
title('Front Guassian Filter');
set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
xlabel('x (m)');
ylabel('y (m)');
colormap('gray');
axis square;

% 
% B = rescale(guassian_filtered_image)
% level = graythresh(B)
% BW = im2bw(B,.2);
% 
% figure;
% imagesc(tic_x,tic_y,BW );
% title('Front Guassian Filter');
% set(gca,'XTick',linspace(min(tic_x),max(tic_x),3));
% set(gca,'YTick',linspace(min(tic_y),max(tic_y),3));
% xlabel('x (m)');
% ylabel('y (m)');
% colormap('gray');
% axis square;
end


function psf = definePsf(w,b,slope)

x = linspace(-1,1,2*w);
y = linspace(-1,1,2*w);

z = linspace(0,2,2*b);

[grid_z,grid_y,grid_x] = ndgrid(z,y,x);

psf = abs(((4.*slope).^2).*(grid_x.^2 + grid_y.^2) - grid_z);
psf = double(psf == repmat(min(psf,[],1),[2.*b 1 1]));
psf = psf./sum(psf(:,w,w));
psf = psf./norm(psf(:));
psf = circshift(psf,[0 w w]);

end


function [mtx, mtxi] = resampleOperator(temporal_resolution)


mtx = sparse([],[],[],temporal_resolution^2,temporal_resolution,temporal_resolution^2);

x = 1:temporal_resolution^2;

mtx(sub2ind(size(mtx),x,ceil(sqrt(x)))) = 1;
mtx  = spdiags(1./sqrt(x)',0,temporal_resolution^2,temporal_resolution^2)*mtx;
mtxi = mtx';

I = log(temporal_resolution)./log(2);
for i = 1:round(I)
    mtx  = 0.5.*(mtx(1:2:end,:)  + mtx(2:2:end,:));
    mtxi = 0.5.*(mtxi(:,1:2:end) + mtxi(:,2:2:end));
end


end






