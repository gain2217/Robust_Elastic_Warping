%% APAP-train
% %{
imfolder = 'images\APAP-train';
im_n = 7;
imfile = cell(im_n,1);
for ii = 1:im_n
    imfile{ii} = sprintf('%s\\%1d.png', imfolder, ii - 1);
end

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2; 1,3; 2,3; 2,4; 3,4; 3,5; 4,5; 4,6; 4,7; 5,6; 5,7; 6,7];

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        im{ii} = imresize(im{ii}, scale);

        imsize(ii,:) = size(im{ii});
    end
end

mosaic = REW_mosaic( im, edge_list, 0, 'equi', 0.05, imfolder );
%}