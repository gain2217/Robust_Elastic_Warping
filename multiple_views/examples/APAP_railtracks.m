%% APAP-railtracks
% %{
imfolder = 'images\APAP-railtracks';
im_n = 2;
imfile = cell(im_n,1);
imfile{1} = [imfolder '\' 'P1010517.JPG'];
imfile{2} = [imfolder '\' 'P1010520.JPG'];

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2];

imsize = zeros(im_n,3);

% scale the large images for faster computing, if needed
% for ii = 1:im_n
%     imsize(ii,:) = size(im{ii});
%     if imsize(ii,1) > 720
%         scale = 720/size(im{ii}, 1);
%         im{ii} = imresize(im{ii}, scale);
% 
%         imsize(ii,:) = size(im{ii});
%     end
% end

refi = 1; % the index of the reference image, 0 for automatic straithtening
mosaic = REW_mosaic( im, edge_list, refi, 'persp', 0, imfolder );
%}