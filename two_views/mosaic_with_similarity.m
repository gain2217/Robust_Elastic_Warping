for no_rotation = [false, true]
    
    if no_rotation
        A_S = zeros(2*n,3);
        b_S = zeros(2*n,1);
        
        A_S(1:2:2*n,1) = x1;
        A_S(1:2:2*n,2) = 1;
        A_S(2:2:2*n,1) = y1;
        A_S(2:2:2*n,3) = 1;
        b_S(1:2:2*n) = x2;
        b_S(2:2:2*n) = y2;
        
        s = A_S \ b_S;
        
        S = [s(1) 0 s(2)
            0 s(1) s(3)];
    else
        A_S = zeros(2*n,4);
        b_S = zeros(2*n,1);
        
        A_S(1:2:2*n,1) = x1;
        A_S(1:2:2*n,2) = -y1;
        A_S(1:2:2*n,3) = 1;
        A_S(2:2:2*n,1) = y1;
        A_S(2:2:2*n,2) = x1;
        A_S(2:2:2*n,4) = 1;
        b_S(1:2:2*n) = x2;
        b_S(2:2:2*n) = y2;
        
        s = A_S \ b_S;
        
        S = [s(1) -s(2) s(3)
            s(2) s(1) s(4)];
    end
    
    H_S = [S;0 0 1];
    
    % Mosaic
    box2 = [1  size(im2,2) size(im2,2)  1 ;
        1  1           size(im2,1)  size(im2,1) ;
        1  1           1            1 ] ;
    box2_ = H_S \ box2 ;
    box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
    box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
    u0 = min([1 box2_(1,:)]);
    u1 = max([size(im1,2) box2_(1,:)]);
    ur = u0:u1;
    v0 = min([1 box2_(2,:)]);
    v1 = max([size(im1,1) box2_(2,:)]) ;
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
    box1_2 = H * box1 ;
    box1_2(1,:) = box1_2(1,:) ./ box1_2(3,:) ;
    box1_2(2,:) = box1_2(2,:) ./ box1_2(3,:) ;
    
    sub_u0_ = max([1, min(box1_2(1,:))]);
    sub_u1_ = min([imsize2(2), max(box1_2(1,:))]);
    sub_v0_ = max([1, min(box1_2(2,:))]);
    sub_v1_ = min([imsize2(1), max(box1_2(2,:))]);
    
    % deform image
    gx = zeros(mosaich, mosaicw);
    hy = zeros(mosaich, mosaicw);
    [u,v] = meshgrid(ur,vr) ;
    
    zH_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
    uH_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ zH_ ;
    vH_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ zH_ ;
    uS_ = S(1,1) * u + S(1,2) * v + S(1,3) ;
    vS_ = S(2,1) * u + S(2,2) * v + S(2,3) ;
    
    mu_S = uS_./imsize2(2);
    mu_S(mu_S < 0) = 0;
    mu_S(mu_S > 1) = 1;
    mu_H = 1 - mu_S;
    u_ = mu_H .* uH_ + mu_S .* uS_;
    v_ = mu_H .* vH_ + mu_S .* vS_;
    
    % compute the transformation for the reference image
    H_inv = inv(H);
    z = H_inv(3,1) * u_ + H_inv(3,2) * v_ + H_inv(3,3) ;
    u = (H_inv(1,1) * u_ + H_inv(1,2) * v_ + H_inv(1,3)) ./ z ;
    v = (H_inv(2,1) * u_ + H_inv(2,2) * v_ + H_inv(2,3)) ./ z ;
    
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
    
    % smooth tansition to global transformationsub_u0_ = sub_u0_ + min(gxn);
    sub_u0_ = sub_u0_ + min(gxn);
    sub_u1_ = sub_u1_ + max(gxn);
    sub_v0_ = sub_v0_ + min(hyn);
    sub_v1_ = sub_v1_ + max(hyn);
    dist_horizontal = max(sub_u0_-uH_, uH_-sub_u1_);
    dist_vertical = max(sub_v0_-vH_, vH_-sub_v1_);
    dist_sub = max(dist_horizontal, dist_vertical);
    dist_sub = max(0, dist_sub);
    eta = (eta_d1 - dist_sub) ./ (eta_d1 - eta_d0);
    eta(dist_sub < eta_d0) = 1;
    eta(dist_sub > eta_d1) = 0;
    gx = gx .* eta;
    hy = hy .* eta;
    
    u_ = u_ - gx;
    v_ = v_ - hy;
    
    % mosaiking
    if exist('vl_imwbackward','file')
        im1_p = vl_imwbackward(im2double(im1),u,v) ;
    else
        im1_p = zeros(mosaich,mosaicw,im_ch);
        for kc = 1:size(im1,3)
            im1_p(:,:,kc) = interp2(im2double(im1(:,:,kc)),u,v);
        end
    end
    
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
    im1_p(~mask1) = 1 ;% white background
    im2_p(~mask2) = 1 ;% white background
    mosaic(mass==0) = 1;% white background
    
    figure;
    imshow(mosaic, 'border', 'tight') ;
    
    if no_rotation
        imwrite(im1_p,  [exp_path, '/', 'im1_p_0.jpg']) ;
        imwrite(im2_p,  [exp_path, '/', 'im2_p_0.jpg']) ;
        imwrite(mosaic,  [exp_path, '/', 'mosaic_similarity_0.jpg']) ;
    else
        imwrite(im1_p,  [exp_path, '/', 'im1_p_1.jpg']) ;
        imwrite(im2_p,  [exp_path, '/', 'im2_p_1.jpg']) ;
        imwrite(mosaic,  [exp_path, '/', 'mosaic_similarity_1.jpg']) ;
    end
    
end