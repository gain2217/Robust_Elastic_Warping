%% APAP-conssite
% %{
imfolder = 'images\APAP-conssite';
im_n = 5;
imfile = cell(im_n,1);
imfile{1} = [imfolder '\' 'DSC03002.jpg'];
imfile{2} = [imfolder '\' 'DSC03003.jpg'];
imfile{3} = [imfolder '\' 'DSC03004.jpg'];
imfile{4} = [imfolder '\' 'DSC03005.jpg'];
imfile{5} = [imfolder '\' 'DSC03006.jpg'];

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2; 1,3; 1,4; 2,3; 2,4; 2,5; 3,4; 3,5; 4,5];

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        im{ii} = imresize(im{ii}, scale);

        imsize(ii,:) = size(im{ii});
    end
end

mosaic = REW_mosaic( im, edge_list, 1, 'equi', 0, imfolder );
%}