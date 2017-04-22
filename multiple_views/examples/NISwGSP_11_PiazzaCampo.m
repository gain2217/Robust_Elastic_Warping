%% NISwGSP-11_PiazzaCampo
% %{
imfolder = 'images\NISwGSP-11_PiazzaCampo';
im_n = 18;
imfile = cell(im_n,1);
for ii = 1:im_n
    imfile{ii} = sprintf('%s\\image%02d.jpg', imfolder, ii);
end

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = zeros(im_n-3,2);
ei = 0;
for ii = 3:im_n-1
    ei = ei + 1;
    edge_list(ei,:) = [ii,ii+1];
end
edge_list = cat(1, edge_list, [3,18; 1,2; 1,18; 1,17; 1,3; 2,17; 2,16; 2,18]);

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