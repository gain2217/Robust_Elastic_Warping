%% NISwGSP-cise
% %{
imfolder = 'images\NISwGSP-cise';
im_n = 11;
imfile = cell(im_n,1);
imfile{1} = [imfolder '\' '00-IMG_8711.JPG'];
imfile{2} = [imfolder '\' '01-IMG_8709.JPG'];
imfile{3} = [imfolder '\' '02-IMG_8712.JPG'];
imfile{4} = [imfolder '\' '03-IMG_8710.JPG'];
imfile{5} = [imfolder '\' '04-IMG_8713.JPG'];
imfile{6} = [imfolder '\' '05-IMG_8714.JPG'];
imfile{7} = [imfolder '\' '06-IMG_8715.JPG'];
imfile{8} = [imfolder '\' '07-IMG_8716.JPG'];
imfile{9} = [imfolder '\' '08-IMG_8717.JPG'];
imfile{10} = [imfolder '\' '09-IMG_8718.JPG'];
imfile{11} = [imfolder '\' '10-IMG_8719.JPG'];

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = zeros(im_n-1,2);
ei = 0;
for ii = 1:im_n-1
    ei = ei + 1;
    edge_list(ei,:) = [ii,ii+1];
end

imsize = zeros(im_n,3);

for ii = 1:im_n
    imsize(ii,:) = size(im{ii});
    if imsize(ii,1) > 720
        scale = 720/size(im{ii}, 1);
        im{ii} = imresize(im{ii}, scale);

        imsize(ii,:) = size(im{ii});
    end
end

mosaic = REW_mosaic( im, edge_list, 0, 'equi', 0.04, imfolder );
%}