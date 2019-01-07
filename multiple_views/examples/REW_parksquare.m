%% REW_songyalake_square2
% %{
imfolder = 'images\REW_parksquare';
im_n = 2;
imfile = cell(im_n,1);
imfile{1} = [imfolder '\' 'parksquare_01.jpg'];
imfile{2} = [imfolder '\' 'parksquare_02.jpg'];

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2];

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        im{ii} = imresize(im{ii}, scale);

        imsize(ii,:) = size(im{ii});
    end
end

mosaic = REW_mosaic( im, edge_list, 1, 'persp', 0, imfolder );
%}
