%% NISwGSP-13_PalazzoPubblico2
% %{
imfolder = 'images\NISwGSP-13_PalazzoPubblico2';
im_n = 11;
imfile = cell(im_n,1);
for ii = 1:im_n
    imfile{ii} = sprintf('%s\\image%02d.jpg', imfolder, ii);
end

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
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

mosaic = REW_mosaic( im, [], 0, 'equi', 0.01, imfolder );
%}