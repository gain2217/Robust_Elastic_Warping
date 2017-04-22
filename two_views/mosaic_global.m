% Mosaic
box2 = [1  size(im2,2) size(im2,2)  1 ;
        1  1           size(im2,1)  size(im2,1) ;
        1  1           1            1 ] ;
box2_ = H \ box2 ;
box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
u0 = min([1 box2_(1,:)]);
u1 = max([size(im1,2) box2_(1,:)]);
ur = u0:u1;
v0 = min([1 box2_(2,:)]);
v1 = max([size(im1,1) box2_(2,:)]) ;
% v0 = max(-500, v0); v1 = min(size(im1,1)+500, v1); %% restrict the size for fast computation
vr = v0:v1;
mosaicw = size(ur, 2);
mosaich = size(vr, 2);

[u,v] = meshgrid(ur,vr) ;
if exist('vl_imwbackward','file')
    im1_p = vl_imwbackward(im2double(im1),u,v) ;
else
    im1_p = zeros(mosaich,mosaicw,im_ch);
    for kc = 1:size(im1,3)
        im1_p(:,:,kc) = interp2(im2double(im1(:,:,kc)),u,v);
    end
end

z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
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
im1_p(isnan(im1_p)) = 0 ;
im2_p(isnan(im2_p)) = 0 ;
mosaic = (im1_p + im2_p) ./ mass ;
mosaic(mass==0) = 1;% white background

figure;
imshow(mosaic, 'border', 'tight') ;

imwrite(mosaic, [exp_path, '/', 'mosaic_global.jpg']) ;