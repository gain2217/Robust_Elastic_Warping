%% CAVE-seaport (not good)
% %{
imfolder = 'images\CAVE-seaport';
im_n = 15;
imfile = cell(im_n,1);
for ii = 1:im_n
    imfile{ii} = sprintf('%s\\%02d.jpg', imfolder, ii-1);
end

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2; 1,4; 1,14; 1,15; 2,3; 2,5; 3,6; 3,9; 3,13; 4,5; 4,6; 4,8; 5,6; 5,7; 6,7; 6,9; 7,8; 7,9;...
    9,11; 10,11; 10,12; 11,12; 11,13; 12,13; 13,14; 14,15];

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        im{ii} = imresize(im{ii}, scale);

        imsize(ii,:) = size(im{ii});
    end
end

mosaic = REW_mosaic( im, edge_list, 0, 'equi', 0.02, imfolder );
%}