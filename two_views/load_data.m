close all; clear all;
run('../vlfeat-0.9.20/toolbox/vl_setup');
addpath('transforms');

% exp_path = 'images/APAP-railtracks';
% imfile1 = [exp_path '/' 'railtracks_01.jpg'];
% imfile2 = [exp_path '/' 'railtracks_02.jpg'];
% exp_path = 'images/DHW-temple';
% imfile1 = [exp_path '/' '4.jpg'];
% imfile2 = [exp_path '/' '5.jpg'];
% exp_path = 'images/ANAP-intersection';
% imfile1 = [exp_path '/' 'intersection_01.jpg'];
% imfile2 = [exp_path '/' 'intersection_02.jpg'];
% exp_path = 'images/REW_workingtable';
% imfile1 = [exp_path '/' 'workingtable_01.jpg'];
% imfile2 = [exp_path '/' 'workingtable_02.jpg'];

exp_path = 'images/REW_theater';
imfile1 = [exp_path '/' 'theater_01.jpg'];
imfile2 = [exp_path '/' 'theater_02.jpg'];

im1 = imread(imfile1);
im2 = imread(imfile2);

imsize1 = size(im1);
imsize2 = size(im2);

% scale the input images if a smaller size is needed
% if size(im1, 1) > 960
%     scale = 960/size(im1, 1);
%     im1 = imresize(im1, scale);
%     
%     imsize1 = size(im1);
% end
% if size(im2, 1) > 960
%     scale = 960/size(im2, 1);
%     im2 = imresize(im2, scale);
%     
%     imsize2 = size(im2);
% end

im_ch = size(im1,3);
if im_ch > 1
    gray1 = im2single(rgb2gray(im1));
    gray2 = im2single(rgb2gray(im2));
elseif im_ch == 1
    gray1 = im2single(im1);
    gray2 = im2single(im2);
end

if exist('vl_sift', 'file')
    [ kp1,ds1 ] = vl_sift(single(gray1),'PeakThresh', 0,'edgethresh',500);
    [ kp2,ds2 ] = vl_sift(single(gray2),'PeakThresh', 0,'edgethresh',500);
    matches = vl_ubcmatch(ds1,ds2);
    
    X1 = [ kp1(1:2,matches(1,:)) ; ones(1,size(matches,2)) ];
    X2 = [ kp2(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
else
    metric_threshold = 200;
    num_octaves = 3;
    num_scale_levels = 4;
    points1 = detectSURFFeatures(gray1, 'MetricThreshold', metric_threshold, 'NumOctaves', num_octaves, 'NumScaleLevels', num_scale_levels);
    points2 = detectSURFFeatures(gray2, 'MetricThreshold', metric_threshold, 'NumOctaves', num_octaves, 'NumScaleLevels', num_scale_levels);
    
    [features1, valid_points1] = extractFeatures(gray1, points1);
    [features2, valid_points2] = extractFeatures(gray2, points2);
    
    [indexPairs,matchmetric] = matchFeatures(features1, features2);
    
    matched_points1 = valid_points1(indexPairs(:, 1), :);
    matched_points2 = valid_points2(indexPairs(:, 2), :);
    
    X1 = matched_points1.Location';
    X2 = matched_points2.Location';
    X1 = [X1;ones(1,size(X1, 2))];
    X2 = [X2;ones(1,size(X2, 2))];
end