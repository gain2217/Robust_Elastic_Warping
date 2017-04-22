%% SPHP-building
% %{
imfolder = 'images\SPHP-building';
im_n = 3;
imfile = cell(im_n,1);
imfile{1} = [imfolder '\' 'building_01.jpg'];
imfile{2} = [imfolder '\' 'building_02.jpg'];
imfile{3} = [imfolder '\' 'building_03.jpg'];

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2; 2,3];

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        im{ii} = imresize(im{ii}, scale);

        imsize(ii,:) = size(im{ii});
    end
end

mosaic = REW_mosaic( im, edge_list, 2, 'persp', 0, imfolder );
%}