% Parameters
lambda = 0.001 * imsize1(1)*imsize1(2); % weighting parameter to balance the fitting term and the smoothing term
intv_mesh = 10; % interval in pixels for the computing of deformation functions
K_smooth = 5; % the smooth transition width in the non-overlapping region is set to K_smooth times of the maximum bias.

% Mosaic
box1 = [1:size(im1,2)        1:size(im1,2)                    ones(1,size(im1,1))  size(im1,2)*ones(1,size(im1,1)) ;
        ones(1,size(im1,2))  size(im1,1)*ones(1,size(im1,2))  1:size(im1,1)        1:size(im1,1) ;
        ones(1,size(im1,2))  ones(1,size(im1,2))              ones(1,size(im1,1))  ones(1,size(im1,1)) ] ;
box2 = [1:size(im2,2)        1:size(im2,2)                    ones(1,size(im2,1))  size(im2,2)*ones(1,size(im2,1)) ;
        ones(1,size(im2,2))  size(im2,1)*ones(1,size(im2,2))  1:size(im2,1)        1:size(im2,1) ;
        ones(1,size(im2,2))  ones(1,size(im2,2))              ones(1,size(im2,1))  ones(1,size(im2,1)) ] ;

fe = max(M1(1,1), M2(1,1));
box1_ = zeros(size(box1));
box2_ = zeros(size(box2));
[box1_(1,:), box1_(2,:)] =  trans_persp2equi(box1(1,:), box1(2,:), eye(3), M1, D1, fe);
[box2_(1,:), box2_(2,:)] =  trans_persp2equi(box2(1,:), box2(2,:), R', M2, D2, fe);

u0 = min([box1_(1,:), box2_(1,:)]);
u1 = max([box1_(1,:), box2_(1,:)]);
ur = u0:u1;
v0 = min([box1_(2,:), box2_(2,:)]);
v1 = max([box1_(2,:), box2_(2,:)]);
vr = v0:v1;
mosaicw = size(ur, 2);
mosaich = size(vr, 2);

% align the sub coordinates with the mosaic coordinates
margin = 0.2 * min(imsize1(1),imsize1(2)); % additional margin of the reprojected image region
u0_im_ = max(min(box2_(1,:)) - margin, u0);
u1_im_ = min(max(box2_(1,:)) + margin, u1);
v0_im_ = max(min(box2_(2,:)) - margin, v0);
v1_im_ = min(max(box2_(2,:)) + margin, v1);
offset_u0_ = ceil(u0_im_ - u0 + 1);
offset_u1_ = floor(u1_im_ - u0 + 1);
offset_v0_ = ceil(v0_im_ - v0 + 1);
offset_v1_ = floor(v1_im_ - v0 + 1);
imw_ = floor(offset_u1_ - offset_u0_ + 1);
imh_ = floor(offset_v1_ - offset_v0_ + 1);

% boundaries of the overlapping region in the image coordiantes of image 2
box1_2 = zeros(size(box1));
[box1_2(1,:), box1_2(2,:)] =  trans_persp2persp(box1(1,:), box1(2,:), R, M1, D1, M2, D2);
sub_u0_ = max([1, min(box1_2(1,:))]);
sub_u1_ = min([imsize2(2), max(box1_2(1,:))]);
sub_v0_ = max([1, min(box1_2(2,:))]);
sub_v1_ = min([imsize2(1), max(box1_2(2,:))]);

% TPS
% merge the coincided points
ok_nd1 = false(size(X1_ok,2),1);
[~, idx1] = unique(round(X1_ok'), 'rows', 'stable');
ok_nd1(idx1) = true;
ok_nd2 = false(size(X2_ok,2),1);
[~, idx2] = unique(round(X2_ok'), 'rows', 'stable');
ok_nd2(idx2) = true;
ok_nd = ok_nd1 & ok_nd2;
X1_nd = X1_ok(:,ok_nd);
X2_nd = X2_ok(:,ok_nd);

% form the linear system
x1 = X1_nd(1,:);
y1 = X1_nd(2,:);
x2 = X2_nd(1,:);
y2 = X2_nd(2,:);
% metrix = metrix_nd;

[x1_, y1_] = trans_persp2persp(x1, y1, R, M1, D1, M2, D2);
gxn = x1_ - x2;
hyn = y1_ - y2;

n = size(x1_, 2);
xx = x1_(ones(1,n),:);
yy = y1_(ones(1,n),:);
dist2 = (xx - xx').^2 + (yy - yy').^2;
dist2(1:n+1:n*n) = ones(1,n);
K = 0.5 * dist2 .* log(dist2);
% K(1:n+1:n*n) = lambda * 8*pi * ones(1,n) ./ metrix_ok;
K(1:n+1:n*n) = lambda * 8*pi * ones(1,n);
K_ = zeros(n+3,n+3);
K_(1:n,1:n) = K;
K_(n+1,1:n) = x1_;
K_(n+2,1:n) = y1_;
K_(n+3,1:n) = ones(1,n);
K_(1:n,n+1) = x1_';
K_(1:n,n+2) = y1_';
K_(1:n,n+3) = ones(n,1);
G_ = zeros(n+3,2);
G_(1:n,1) = gxn';
G_(1:n,2) = hyn';

% solve the linear system
W_ = K_\G_;
wx = W_(1:n,1);
wy = W_(1:n,2);
a = W_(n+1:n+3,1);
b = W_(n+1:n+3,2);

% remove outliers based on the distribution of weights
outlier = abs(wx)>3*std(wx) | abs(wy)>3*std(wy);

inlier_idx = 1:size(x1, 2);
for kiter = 1:10
%     if ~any(outlier)
    if sum(outlier) < 0.0027*n
        break;
    end
    ok = ~outlier;
    inlier_idx = inlier_idx(ok);
    K_ = K_([ok;true(3,1)],[ok;true(3,1)]);
    G_ = G_([ok;true(3,1)],:);
    W_ = K_\G_;
    n = size(inlier_idx,2);
    wx = W_(1:n,1);
    wy = W_(1:n,2);
    a = W_(n+1:n+3,1);
    b = W_(n+1:n+3,2);
    outlier = abs(wx)>3*std(wx) | abs(wy)>3*std(wy);  
end
ok = false(size(x1, 2),1);
ok(inlier_idx) = true;
x1 = x1(ok);
y1 = y1(ok);
x2 = x2(ok);
y2 = y2(ok);
x1_ = x1_(ok);
y1_ = y1_(ok);
% metrix = metrix(ok);

figure; clf;
marker_size = 15;
imshow([im1,im2], 'border', 'tight'); hold on;
plot(X1_nd(1,~ok), X1_nd(2,~ok), 'r.', 'MarkerSize', marker_size);
plot(X1_nd(1,ok), X1_nd(2,ok), 'b.', 'MarkerSize', marker_size);

plot(X2_nd(1,~ok)+imsize1(2), X2_nd(2,~ok), 'r.', 'MarkerSize', marker_size);
plot(X2_nd(1,ok)+imsize1(2), X2_nd(2,ok), 'b.', 'MarkerSize', marker_size);

% deform image
gx = zeros(mosaich, mosaicw);
hy = zeros(mosaich, mosaicw);
[u,v] = meshgrid(ur,vr) ;

[u_, v_] = trans_equi2persp(u, v, eye(3), M1, D1, fe);
if exist('vl_imwbackward','file')
    im1_p = vl_imwbackward(im2double(im1),u_,v_) ;
else
    im1_p = zeros(mosaich,mosaicw,im_ch);
    for kc = 1:size(im1,3)
        im1_p(:,:,kc) = interp2(im2double(im1(:,:,kc)),u_,v_);
    end
end

[u_, v_] = trans_equi2persp(u, v, R, M2, 0, fe);
u_im_ = u_(offset_v0_:intv_mesh:offset_v1_,offset_u0_:intv_mesh:offset_u1_);
v_im_ = v_(offset_v0_:intv_mesh:offset_v1_,offset_u0_:intv_mesh:offset_u1_);
gx_sub = zeros(ceil(imh_/intv_mesh), ceil(imw_/intv_mesh));
hy_sub = zeros(ceil(imh_/intv_mesh), ceil(imw_/intv_mesh));
for kf = 1:n
    dist2 = (u_im_ - x1_(kf)).^2 + (v_im_ - y1_(kf)).^2;
    rbf = 0.5 * dist2 .* log(dist2);
    gx_sub = gx_sub + wx(kf)*rbf;
    hy_sub = hy_sub + wy(kf)*rbf;
end
gx_sub = gx_sub + a(1).*u_im_+a(2).*v_im_+a(3);
hy_sub = hy_sub + b(1).*u_im_+b(2).*v_im_+b(3);
gx_sub = imresize(gx_sub, [imh_,imw_]);
hy_sub = imresize(hy_sub, [imh_,imw_]);
gx(offset_v0_:offset_v1_,offset_u0_:offset_u1_) = gx_sub;
hy(offset_v0_:offset_v1_,offset_u0_:offset_u1_) = hy_sub;

%smooth tansition to global transform
eta_d0 = 0; % lower boundary for smooth transition area
eta_d1 = K_smooth * max(abs([gxn, hyn])); % higher boundary for smooth transition area
dist_horizontal = max(sub_u0_-u_, u_-sub_u1_);
dist_vertical = max(sub_v0_-v_, v_-sub_v1_);
dist_sub = max(dist_horizontal, dist_vertical);
dist_sub = max(0, dist_sub);
eta = (eta_d1 - dist_sub) ./ (eta_d1 - eta_d0);
eta(dist_sub < eta_d0) = 1;
eta(dist_sub > eta_d1) = 0;
gx = gx .* eta;
hy = hy .* eta;

u_ = u_ - gx;
v_ = v_ - hy;

if exist('vl_imwbackward','file')
    im2_p = vl_imwbackward(im2double(im2),u_,v_) ;
else
    im2_p = zeros(mosaich,mosaicw,im_ch);
    for kc = 1:size(im2,3)
        im2_p(:,:,kc) = interp2(im2double(im2(:,:,kc)),u_,v_);
    end
end

mask1 = ~isnan(im1_p);
mask2 = ~isnan(im2_p);
mass = mask1 + mask2 ;
im1_p(~mask1) = 0 ;
im2_p(~mask2) = 0 ;
mosaic = (im1_p + im2_p) ./ mass ;
mosaic(mass==0) = 1;% white background

figure;
imshow(mosaic, 'border', 'tight') ;
imwrite(mosaic,  [exp_path, '/', 'mosaic_ours_equi.jpg']) ;