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
[box1_(1,:), box1_(2,:)] =  trans_persp2equi(box1(1,:), box1(2,:), eye(3), M1, 0, fe);
[box2_(1,:), box2_(2,:)] =  trans_persp2equi(box2(1,:), box2(2,:), R', M2, 0, fe);

u0 = min([box1_(1,:), box2_(1,:)]);
u1 = max([box1_(1,:), box2_(1,:)]);
ur = u0:u1;
v0 = min([box1_(2,:), box2_(2,:)]);
v1 = max([box1_(2,:), box2_(2,:)]);
vr = v0:v1;
mosaicw = size(ur, 2);
mosaich = size(vr, 2);
[u,v] = meshgrid(ur,vr) ;

[u_, v_] = trans_equi2persp(u, v, eye(3), M1, 0, fe);
if exist('vl_imwbackward','file')
    im1_p = vl_imwbackward(im2double(im1),u_,v_) ;
else
    im1_p = zeros(mosaich,mosaicw,im_ch);
    for kc = 1:size(im1,3)
        im1_p(:,:,kc) = interp2(im2double(im1(:,:,kc)),u_,v_);
    end
end

[u_, v_] = trans_equi2persp(u, v, R, M2, 0, fe);
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
imwrite(mosaic,  [exp_path, '/', 'mosaic_global_equi.jpg']) ;