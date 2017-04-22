%% ANAP-skyline
% %{
imfolder = 'images\ANAP-skyline';
im_n = 4;
imfile = cell(im_n,1);
imfile{1} = [imfolder '\' 'skyline0.png'];
imfile{2} = [imfolder '\' 'skyline1.png'];
imfile{3} = [imfolder '\' 'skyline2.png'];
imfile{4} = [imfolder '\' 'skyline3.png'];

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2; 1,3; 1,4; 2,3; 2,4; 3,4];

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        im{ii} = imresize(im{ii}, scale);

        imsize(ii,:) = size(im{ii});
    end
end

mosaic = REW_mosaic( im, edge_list, 4, [], 0, imfolder );
%}