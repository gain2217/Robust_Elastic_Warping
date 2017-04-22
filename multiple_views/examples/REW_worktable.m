%% REW_worktable
% %{
imfolder = 'images\REW_worktable'; % The folder containing the input images,
                                      % and saving the resulted mosaic, if needed.
im_n = 2; % The number of input images
imfile = cell(im_n,1); % The file paths of the input images
imfile{1} = [imfolder '\' 'worktable_01.jpg'];
imfile{2} = [imfolder '\' 'worktable_02.jpg'];

im = cell(im_n,1);
for ii = 1:im_n
    im{ii} = imread(imfile{ii});
end

edge_list = [1,2]; % Each row represents an image pair to be aligned.
                   % [] for automatic detection. 

% imsize = zeros(im_n,3);
% 
% for ii = 1:im_n
%     imsize(ii,:) = size(im{ii});
%     if imsize(ii,1) > 720
%         scale = 720/size(im{ii}, 1);
%         im{ii} = imresize(im{ii}, scale);
% 
%         imsize(ii,:) = size(im{ii});
%     end
% end

refi = 1; % The index of the reference image, 0 for automatic straithtening
projection_type = 'persp'; % the projection type of the mosaic, 
                           % 'persp' (default) for perspective projection
                           % and 'equi' for equirectangular projection.
                           % The equirectangular projection is recommended
                           % for mosaic with large field-of-view.
ransac_threshold = 0.1;% threshold of global ransac.
                       % 0 for default value which is set to 0.1.
  
mosaic = REW_mosaic( im, edge_list, refi, projection_type, ransac_threshold, imfolder );
%}
