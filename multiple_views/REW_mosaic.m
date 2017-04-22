function [ mosaic ] = REW_mosaic( im, edge_list, refi, projection_type, ransac_threshold, imfolder )
% Robust Elastic Warping for Parallax-Tolerant Image Stitcing
%   im : the input images
%   edge_list: the overllapped edge list, each row represents an image pair to be aligned
%   refi: the index of the reference image, 0 for automatic straithtening
%   projection_type: the projection type of the mosaic, 
%                    'persp' (default) for perspective projection 
%                    and 'equi' for equirectangular projection.
%   ransac_threshold£º threshold of global ransac
%                     0 for default value which is set to 0.1;
%   imfolder£º The folder containing the input images,
%              and saving the resulted mosaic, if needed.


if nargin <= 5
    imfolder = 'temp';
    if ~exist(imfolder, 'dir')
        mkdir(imfolder) ;
    end
end
if nargin <= 4
    ransac_threshold = 0.1;
elseif ransac_threshold == 0
    ransac_threshold = 0.1;
end
if nargin <= 3
    projection_type = 'persp';
elseif isempty(projection_type)
    projection_type = 'persp';
end
if ~strcmp(projection_type,'persp') && ~strcmp(projection_type,'equi')
    error('Unsupported projection type!');
end

save_results = true;
recomp_global_paras = false;
compute_global_results = true;
show_intermediate_results = false;
blend_output = true;
bgcolor = 1; % 0 for black, 1 for white

im_n = size(im, 1);
edge_n = size(edge_list, 1);

imsize = zeros(im_n,3);
for i = 1:im_n
    imsize(i,:) = size(im{i});
end

%% Parameters
lambda = 0.001 * imsize(1,1)*imsize(1,2); % weighting parameter to balance the fitting term and the smoothing term
intv_mesh = 10; % interval in pixels for the computing of deformation functions
K_smooth = 5; % the smooth transition width in the non-overlapping region is set to K_smooth times of the maximum bias.

%% feature detection and matching
im_gray = cell(im_n, 1);
points = cell(im_n, 1);
features = cell(im_n, 1);
valid_points = cell(im_n, 1);

% detection
for i = 1:im_n
    im_ch = size(im{i},3);
    if im_ch > 1
        im_gray{i} = im2single(rgb2gray(im{i}));
    elseif im_ch == 1
        im_gray{i} = im2single(im{i});
    end
    if exist('vl_sift', 'file')
        [ points{i},features{i} ] = vl_sift(single(im_gray{i}),'PeakThresh', 0,'edgethresh',500);
    else
        points{i} = detectSURFFeatures(im_gray{i},...
            'MetricThreshold', 0, 'NumOctaves', 3, 'NumScaleLevels', 4);
        [features{i}, valid_points{i}] = extractFeatures(im_gray{i}, points{i});
    end
end

% matching
if edge_n == 0
    % image match verfication referring to AutoStitch
    X = cell(im_n*(im_n-1)/2, 2);
    for i = 1:im_n-1
        for j = i+1:im_n
            if exist('vl_sift', 'file')
                matches = vl_ubcmatch(features{i}, features{j});
                X_1 = [ points{i}(1:2,matches(1,:)) ; ones(1,size(matches,2)) ];
                X_2 = [ points{j}(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
            else
                indexPairs = matchFeatures(features{i}, features{j});
                matched_points_1 = valid_points{i}(indexPairs(:, 1), :);
                matched_points_2 = valid_points{j}(indexPairs(:, 2), :);
                X_1 = [matched_points_1.Location'; ones(1,size(matched_points_1, 1))];
                X_2 = [matched_points_2.Location'; ones(1,size(matched_points_2, 1))];
            end
            nf = size(X_1, 2);
            if nf <= 14 % minimum allown matches number
                continue;
            end
            [ ~, ok, score ] = HM_ransac(X_1, X_2, 200, ransac_threshold);
            if score > 8 + 0.4 * nf
                edge_n = edge_n + 1;
                edge_list = cat(1,edge_list,[i,j]);
                X{edge_n,1} = X_1(:,ok);
                X{edge_n,2} = X_2(:,ok);
            end
        end
    end
    X = X(1:edge_n,:);
else
    X = cell(edge_n, 2);
    for ei = 1 : edge_n
        i = edge_list(ei, 1);
        j = edge_list(ei, 2);
        
        if exist('vl_sift', 'file')
            matches = vl_ubcmatch(features{i}, features{j});
            X_1 = [ points{i}(1:2,matches(1,:)) ; ones(1,size(matches,2)) ];
            X_2 = [ points{j}(1:2,matches(2,:)) ; ones(1,size(matches,2)) ];
        else
            indexPairs = matchFeatures(features{i}, features{j});
            matched_points_1 = valid_points{i}(indexPairs(:, 1), :);
            matched_points_2 = valid_points{j}(indexPairs(:, 2), :);
            X_1 = [matched_points_1.Location'; ones(1,size(matched_points_1, 1))];
            X_2 = [matched_points_2.Location'; ones(1,size(matched_points_2, 1))];
        end
        [ ~, ok, ~ ] = HM_ransac(X_1, X_2, 200, ransac_threshold);
        X{ei,1} = X_1(:,ok);
        X{ei,2} = X_2(:,ok);
    end
end

if show_intermediate_results
    figure(1); clf;
    marker_size = 15;
    en_root = ceil(sqrt(2*edge_n)/2);
    em_root = ceil(edge_n/en_root);
    for ei = 1 : edge_n
        i = edge_list(ei, 1);
        j = edge_list(ei, 2);
        
        subplot(em_root,en_root,ei);
        imshow([im{i},im{j}], 'border', 'tight'); hold on;
        plot(X{ei,1}(1,:), X{ei,1}(2,:), 'b.', 'MarkerSize', marker_size);
        plot(X{ei,2}(1,:)+imsize(i,2), X{ei,2}(2,:), 'b.', 'MarkerSize', marker_size);
        
        line([X{ei,1}(1,:) ; X{ei,2}(1,:)+imsize(i,2)], [X{ei,1}(2,:) ; X{ei,2}(2,:)], 'LineWidth', 1) ;
    end
end

%% global transforms estimation
if exist([imfolder,'\global_paras.mat'],'file') && ~recomp_global_paras
    load([imfolder,'\global_paras.mat']);
else
    options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'Display','final',...
        'MaxFunEvals',1000*im_n, 'MaxIter',1e3, 'TolFun',1e-6, 'TolX',1e-6, 'Jacobian','off');
    paras_init = [1000*ones(1,im_n),zeros(1,3*(im_n-1))];
    for sigma = [1000, 100, 10]
        [paras, ~ ,~ ,exitflag] = lsqnonlin(...
            @(p)residual_all_robust(X, imsize, edge_list, p, sigma), paras_init,...
            [],[],options);
        if exitflag > 0
            paras_init = paras;
        end
    end
    save([imfolder,'\global_paras.mat'],'paras');
end

M = cell(im_n,1);
D = cell(im_n,1);
R_pair = cell(im_n, im_n); % pairwise rotation matrics
for i = 1 : im_n
    ki = paras(i);
    M{i} = [ki, 0, imsize(i,2)/2;
             0, ki, imsize(i,1)/2;
             0,  0, 1];
    D{i} = 0;
end
for i = 1 : im_n
    R_pair{i, i} = eye(3);
end
for i = 2:im_n
    theta = paras(im_n+3*(i-2)+1:im_n+3*(i-2)+3);
    theta_m = [0         -theta(3) theta(2)
               theta(3)  0         -theta(1)
               -theta(2) theta(1)  0];
    R_pair{1,i} = expm(theta_m);
    R_pair{i,1} = R_pair{1,i}';
end
for i = 2:im_n-1
    for j = i+1:im_n
        R_pair{i,j} = R_pair{1,j}*R_pair{i,1};
        R_pair{j,i} = R_pair{i,j}';
    end
end

if refi == 0
    % automatic straithtening
    xz_vecs = zeros(im_n*2,3);
    for i = 1:im_n
        xz_vecs(2*i-1:2*i,:) = R_pair{1,i}([1,3],:);
    end
    [~,~,V] = svd(xz_vecs,'econ');
    y_vec_glb = V(:,3);
    if y_vec_glb(2) < 0
        y_vec_glb = -y_vec_glb;
    end
    x_temp = zeros(3,1);
    for i = 1:im_n
        x_temp = x_temp + R_pair{1,i}(1,:)';
    end
    x_temp = x_temp ./ im_n;
    if norm(x_temp) < 0.1
        x_temp = [1;0;0];
    end
    z_vec_glb = cross(x_temp,y_vec_glb);
    z_vec_glb = z_vec_glb/norm(z_vec_glb);
    x_vec_glb = cross(y_vec_glb,z_vec_glb);
    % global rotation matrix for the reference image which is set to the first image here
    R_ref = [x_vec_glb,y_vec_glb,z_vec_glb]'; % global rotation matrix
    
    refi = 1;
else
    R_ref = eye(3);
end
R = cell(im_n, 1); % global rotation matrics
for i = 1 : im_n
    R{i} = R_pair{refi,i}*R_ref';
end

%% computing mosaic parameters
if strcmp(projection_type,'equi')
    fe = max(M{refi}(1,1),M{refi}(2,2));
else
    Mp = M{refi};
    Dp = D{refi};
end
ubox = cell(im_n,1);
vbox = cell(im_n,1);
ubox_ = cell(im_n,1);
vbox_ = cell(im_n,1);
ubox_all_ = [];
vbox_all_ = [];
for i = 1 : im_n
    ubox{i} = [1:imsize(i,2)        1:imsize(i,2)                     ones(1,imsize(i,1))  imsize(i,2)*ones(1,imsize(i,1))] ;
    vbox{i} = [ones(1,imsize(i,2))  imsize(i,1)*ones(1,imsize(i,2))  1:imsize(i,1)        1:imsize(i,1) ];
    if strcmp(projection_type,'equi')
        [ubox_{i}, vbox_{i}] =  trans_persp2equi(ubox{i}, vbox{i}, R{i}', M{i}, D{i}, fe);
    else
        [ubox_{i}, vbox_{i}] =  trans_persp2persp(ubox{i}, vbox{i}, R{i}', M{i}, D{i}, Mp, Dp);
    end
    ubox_all_ = cat(2,ubox_all_,ubox_{i});
    vbox_all_ = cat(2,vbox_all_,vbox_{i});
end
u0 = min(ubox_all_);
u1 = max(ubox_all_);
ur = u0:u1;
v0 = min(vbox_all_);
v1 = max(vbox_all_);
vr = v0:v1;
mosaicw = size(ur, 2);
mosaich = size(vr, 2);

m_u0_ = zeros(im_n,1);
m_u1_ = zeros(im_n,1);
m_v0_ = zeros(im_n,1);
m_v1_ = zeros(im_n,1);
imw_ = zeros(im_n,1);
imh_ = zeros(im_n,1);
for i = 1 : im_n
    % align the sub coordinates with the mosaic coordinates
    margin = 0.2 * min(imsize(1,1),imsize(1,2)); % additional margin of the reprojected image region considering the possilbe deformation
    u0_im_ = max(min(ubox_{i}) - margin, u0);
    u1_im_ = min(max(ubox_{i}) + margin, u1);
    v0_im_ = max(min(vbox_{i}) - margin, v0);
    v1_im_ = min(max(vbox_{i}) + margin, v1);
    m_u0_(i) = ceil(u0_im_ - u0 + 1);
    m_u1_(i) = floor(u1_im_ - u0 + 1);
    m_v0_(i) = ceil(v0_im_ - v0 + 1);
    m_v1_(i) = floor(v1_im_ - v0 + 1);
    imw_(i) = floor(m_u1_(i) - m_u0_(i) + 1);
    imh_(i) = floor(m_v1_(i) - m_v0_(i) + 1);
end

%% global mosaic
if compute_global_results
    [u,v] = meshgrid(ur,vr) ;
    
    im_p = cell(im_n,1);
    mask = cell(im_n,1);
    mass = zeros(mosaich, mosaicw, im_ch);
    mosaic = zeros(mosaich, mosaicw, im_ch);
    for i = 1 : im_n
        u_im = u(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i));
        v_im = v(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i));
        if strcmp(projection_type,'equi')
            [u_im_, v_im_] = trans_equi2persp(u_im, v_im, R{i}, M{i}, D{i}, fe);
        else
            [u_im_, v_im_] = trans_persp2persp(u_im, v_im, R{i}, Mp, Dp, M{i}, D{i});
        end
        im_p{i} = zeros(imh_(i),imw_(i),imsize(i,3));
        for kc = 1:imsize(i,3)
            im_p{i}(:,:,kc) = interp2(im2double(im{i}(:,:,kc)),u_im_,v_im_);
        end
        mask{i} = ~isnan(im_p{i});
        mass(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:)...
            = mass(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:) + mask{i};
        im_p{i}(isnan(im_p{i})) = 0;
        mosaic(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:)...
            = mosaic(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:) + im_p{i};
    end
    mosaic = mosaic ./ mass;
    mosaic(isnan(mosaic)) = bgcolor;
    
    figure(2) ;
    imshow(mosaic, 'border', 'tight') ;
    drawnow;
    if save_results
        imwrite(mosaic, [imfolder, '\mosaic_global.jpg']);
    end
end

%% local mosaic

% elastic local alignment and mosaiking
Adj = zeros(im_n,im_n); % adjacent matrix describing the topological 
                        % relationship of the input images,
                        % the elements of Adj indicate the edge indexes
                        % between the two overlapping images, the 0
                        % elements indicate not overlapped
for ei = 1:edge_n
    i = edge_list(ei, 1);
    j = edge_list(ei, 2);
    Adj(i,j) = ei;
    Adj(j,i) = ei;
end

[u,v] = meshgrid(ur,vr) ;

imi_ = cell(im_n,1); % only for intermediate results
im_p = cell(im_n,1);
mask = cell(im_n,1);
mass = zeros(mosaich, mosaicw, im_ch);
mosaic = zeros(mosaich, mosaicw, im_ch);
for ki = 1:im_n
    i = mod(ki + refi - 2, im_n) + 1; % start from the reference image
    % align image i
    u_im = u(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i));
    v_im = v(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i));
    if strcmp(projection_type,'equi')
        [u_im_, v_im_] = trans_equi2persp(u_im, v_im, R{i}, M{i}, D{i}, fe);
    else
        [u_im_, v_im_] = trans_persp2persp(u_im, v_im, R{i}, Mp, Dp, M{i}, D{i});
    end
    
    need_deform = false;
    sub_u0_ = [];
    sub_u1_ = [];
    sub_v0_ = [];
    sub_v1_ = [];
    Pi = [];
    Pi_ = [];
    for kj = 1:ki-1 % for every image that has already been aligned
        j = mod(kj + refi - 2, im_n) + 1;
        if Adj(i,j) > 0
            need_deform = true;
            
            [ubox_ji, vbox_ji] =  trans_persp2persp(ubox{j}, vbox{j}, R_pair{j,i}, M{j}, D{j}, M{i}, D{i});
            sub_u0_ = cat(1,sub_u0_,max([1, min(ubox_ji)]) );
            sub_u1_ = cat(1,sub_u1_,min([imsize(i,2), max(ubox_ji)]) );
            sub_v0_ = cat(1,sub_v0_,max([1, min(vbox_ji)]) );
            sub_v1_ = cat(1,sub_v1_,min([imsize(i,1), max(vbox_ji)]) );
            
            ei = Adj(i,j);
            if i == edge_list(ei, 1) && j == edge_list(ei, 2)
                Xi = X{ei,1};
                Xj = X{ei,2};
            else
                Xi = X{ei,2};
                Xj = X{ei,1};
            end
            [xj_i, yj_i] = trans_persp2persp(Xj(1,:), Xj(2,:), R_pair{j,i}, M{j}, D{j}, M{i}, D{i});
            Pi = cat(2, Pi, Xi(1:2,:));
            Pi_ = cat(2, Pi_, [xj_i;yj_i]);
        end
    end
    if need_deform
        sub_u0_ = min(sub_u0_);
        sub_u1_ = max(sub_u1_);
        sub_v0_ = min(sub_v0_);
        sub_v1_ = max(sub_v1_);
        
        % merge the coincided points
        ok_Pi = false(size(Pi,2),1);
        [~, idx_Pi] = unique(round(Pi'), 'rows', 'stable');
        ok_Pi(idx_Pi) = true;
        ok_Pi_ = false(size(Pi_,2),1);
        [~, idx_Pi_] = unique(round(Pi_'), 'rows', 'stable');
        ok_Pi_(idx_Pi_) = true;
        ok_nd = ok_Pi & ok_Pi_;
        Pi_nd = Pi(:,ok_nd);
        Pi_nd_ = Pi_(:,ok_nd);

        % form the linear system
        xi = Pi_nd(1,:);
        yi = Pi_nd(2,:);
        xj_ = Pi_nd_(1,:);
        yj_ = Pi_nd_(2,:);
        gxn = xj_ - xi;
        hyn = yj_ - yi;
        
        n = size(xj_, 2);
        xx = xj_(ones(1,n),:);
        yy = yj_(ones(1,n),:);
        dist2 = (xx - xx').^2 + (yy - yy').^2;
        dist2(1:n+1:n*n) = ones(1,n);
        K = 0.5 * dist2 .* log(dist2);
        K(1:n+1:n*n) = lambda * 8*pi * ones(1,n);
        K_ = zeros(n+3,n+3);
        K_(1:n,1:n) = K;
        K_(n+1,1:n) = xj_;
        K_(n+2,1:n) = yj_;
        K_(n+3,1:n) = ones(1,n);
        K_(1:n,n+1) = xj_';
        K_(1:n,n+2) = yj_';
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
        
        inlier_idx = 1:size(xi, 2);
        for kiter = 1:10
%             if ~any(outlier)
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
        ok = false(size(xj_, 2),1);
        ok(inlier_idx) = true;
        xj_ = xj_(ok);
        yj_ = yj_(ok);
        gxn = gxn(ok);
        hyn = hyn(ok);
        
        eta_d0 = 0; % lower boundary for smooth transition area
        eta_d1 = K_smooth * max(abs([gxn, hyn])); % higher boundary for smooth transition area
%         eta_d1 = 0.3*max(imsize(1,1:2));
        sub_u0_ = sub_u0_ + min(gxn);
        sub_u1_ = sub_u1_ + max(gxn);
        sub_v0_ = sub_v0_ + min(hyn);
        sub_v1_ = sub_v1_ + max(hyn);
        
        if show_intermediate_results
%             margin = df_max;
            [u_im,v_im] = meshgrid(1:intv_mesh:imsize(i,2),1:intv_mesh:imsize(i,1));
            gx_im = zeros(ceil(imsize(i,1)/intv_mesh),ceil(imsize(i,2)/intv_mesh));
            hy_im = zeros(ceil(imsize(i,1)/intv_mesh),ceil(imsize(i,2)/intv_mesh));
            for kf = 1:n
                dist2 = (u_im - xj_(kf)).^2 + (v_im - yj_(kf)).^2;
                rbf = 0.5 * dist2 .* log(dist2);
                gx_im = gx_im + wx(kf)*rbf;
                hy_im = hy_im + wy(kf)*rbf;
            end
            gx_im = gx_im + a(1).*u_im+a(2).*v_im+a(3);
            hy_im = hy_im + b(1).*u_im+b(2).*v_im+b(3);
            gx_im = imresize(gx_im, [imsize(i,1),imsize(i,2)]);
            hy_im = imresize(hy_im, [imsize(i,1),imsize(i,2)]);
            
            [u_im,v_im] = meshgrid(1:imsize(i,2),1:imsize(i,1)) ;
            dist_horizontal_im = max(sub_u0_-u_im, u_im-sub_u1_);
            dist_vertical_im = max(sub_v0_-v_im, v_im-sub_v1_);
            dist_sub_im = max(dist_horizontal_im, dist_vertical_im);
            dist_sub_im = max(0, dist_sub_im);
            eta_im = (eta_d1 - dist_sub_im) ./ (eta_d1 - eta_d0);
            eta_im(dist_sub_im < eta_d0) = 1;
            eta_im(dist_sub_im > eta_d1) = 0;
            gx_im = gx_im .* eta_im;
            hy_im = hy_im .* eta_im;
            
            u_im = u_im - gx_im;
            v_im = v_im - hy_im;
            imi_{i} = zeros(imsize(i,:));
            for kc = 1:imsize(i,3)
                imi_{i}(:,:,kc) = interp2(im2double(im{i}(:,:,kc)),u_im,v_im);
            end
            
            figure(4);
            intv_showmesh = 30;
            subplot(im_n-1, 2, 2*(ki-2)+1);
            mesh(gx_im(1:intv_showmesh:imsize(i,1), 1:intv_showmesh:imsize(i,2)));
            xlabel('x');
            ylabel('y');
            subplot(im_n-1, 2, 2*(ki-2)+2);
            mesh(hy_im(1:intv_showmesh:imsize(i,1), 1:intv_showmesh:imsize(i,2)));
            xlabel('x');
            ylabel('y');
            
            figure(5);
            subplot(im_n-1, 2, 2*(ki-2)+1);
            imshow(im2double(im{i}), 'border', 'tight');
            subplot(im_n-1, 2, 2*(ki-2)+2);
            imshow(imi_{i}, 'border', 'tight');
        end
        
        u_mesh_ = u_im_(1:intv_mesh:imh_(i),1:intv_mesh:imw_(i));
        v_mesh_ = v_im_(1:intv_mesh:imh_(i),1:intv_mesh:imw_(i));
        gx_mesh_ = zeros(ceil(imh_(i)/intv_mesh), ceil(imw_(i)/intv_mesh));
        hy_mesh_ = zeros(ceil(imh_(i)/intv_mesh), ceil(imw_(i)/intv_mesh));
        for kf = 1:n
            dist2 = (u_mesh_ - xj_(kf)).^2 + (v_mesh_ - yj_(kf)).^2;
            rbf = 0.5 * dist2 .* log(dist2);
            gx_mesh_ = gx_mesh_ + wx(kf)*rbf;
            hy_mesh_ = hy_mesh_ + wy(kf)*rbf;
        end
        gx_mesh_ = gx_mesh_ + a(1).*u_mesh_+a(2).*v_mesh_+a(3);
        hy_mesh_ = hy_mesh_ + b(1).*u_mesh_+b(2).*v_mesh_+b(3);
        gx_im_ = imresize(gx_mesh_, [imh_(i),imw_(i)]);
        hy_im_ = imresize(hy_mesh_, [imh_(i),imw_(i)]);
        
        %smooth tansition to global transform
        dist_horizontal = max(sub_u0_-u_im_, u_im_-sub_u1_);
        dist_vertical = max(sub_v0_-v_im_, v_im_-sub_v1_);
        dist_sub = max(dist_horizontal, dist_vertical);
        dist_sub = max(0, dist_sub);
        eta = (eta_d1 - dist_sub) ./ (eta_d1 - eta_d0);
        eta(dist_sub < eta_d0) = 1;
        eta(dist_sub > eta_d1) = 0;
        gx_im_ = gx_im_ .* eta;
        hy_im_ = hy_im_ .* eta;
        
        u_im_ = u_im_ - gx_im_;
        v_im_ = v_im_ - hy_im_;
        
        % update the feature locations
        for kj = ki+1:im_n % for every image that has not been aligned
            j = mod(kj + refi - 2, im_n) + 1;
            if Adj(i,j) > 0
                ei = Adj(i,j);
                if i == edge_list(ei, 1) && j == edge_list(ei, 2)
                    Xi = X{ei,1};
                else
                    Xi = X{ei,2};
                end
                newXi = Xi;
                % Iterative solver
                for kiter = 1:20
                    u_f = newXi(1,:);
                    v_f = newXi(2,:);
                    gx_f = zeros(1,size(newXi,2));
                    hy_f = zeros(1,size(newXi,2));
                    for kf = 1:n
                        dist2 = (u_f - xj_(kf)).^2 + (v_f - yj_(kf)).^2;
                        rbf = 0.5 * dist2 .* log(dist2);
                        gx_f = gx_f + wx(kf)*rbf;
                        hy_f = hy_f + wy(kf)*rbf;
                    end
                    gx_f = gx_f + a(1).*u_f+a(2).*v_f+a(3);
                    hy_f = hy_f + b(1).*u_f+b(2).*v_f+b(3);
                    dist_horizontal_f = max(sub_u0_-u_f, u_f-sub_u1_);
                    dist_vertical_f = max(sub_v0_-v_f, v_f-sub_v1_);
                    dist_sub_f = max(dist_horizontal_f, dist_vertical_f);
                    dist_sub_f = max(0, dist_sub_f);
                    eta_f = (eta_d1 - dist_sub_f) ./ (eta_d1 - eta_d0);
                    eta_f(dist_sub_f < eta_d0) = 1;
                    eta_f(dist_sub_f > eta_d1) = 0;
                    gx_f = gx_f .* eta_f;
                    hy_f = hy_f .* eta_f;
%                     disp([sum(abs(Xi(1,:) + gx_f - newXi(1,:))),sum(abs(Xi(2,:) + hy_f - newXi(2,:)))]);
                    newXi(1,:) = Xi(1,:) + gx_f;
                    newXi(2,:) = Xi(2,:) + hy_f;
                end
                if i == edge_list(ei, 1) && j == edge_list(ei, 2)
                    X{ei,1} = newXi;
                else
                    X{ei,2} = newXi;
                end
            end
        end
    end
    im_p{i} = zeros(imh_(i),imw_(i),imsize(i,3));
    for kc = 1:imsize(i,3)
        im_p{i}(:,:,kc) = interp2(im2double(im{i}(:,:,kc)),u_im_,v_im_);
    end
    mask{i} = ~isnan(im_p{i});
    mass(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:)...
        = mass(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:) + mask{i};
    im_p{i}(isnan(im_p{i})) = 0;
    mosaic(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:)...
        = mosaic(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:) + im_p{i};
end
mosaic = mosaic ./ mass;
mosaic(isnan(mosaic)) = bgcolor;

figure(6) ;
imshow(mosaic, 'border', 'tight') ;
drawnow;
if save_results
    imwrite(mosaic, [imfolder, '\mosaic_ours.jpg']);
end

%% blending the output mosaic
if blend_output
    for i = 1:im_n
        im_p_save = zeros(mosaich,mosaicw,im_ch);
        im_p_save(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i),:) = im_p{i};
        mask_save = zeros(mosaich,mosaicw);
        mask_save(m_v0_(i):m_v1_(i),m_u0_(i):m_u1_(i)) = mask{i}(:,:,1);
        
        imp_path = sprintf('%s\\p_%02d.png', imfolder, i);
        imwrite(im_p_save, imp_path, 'png', 'Alpha', mask_save) ;
    end
    if bgcolor == 0
        sysorder = sprintf('enblend --output=%s\\mosaic_blend.jpg', imfolder) ;
    else
        sysorder = sprintf('enblend --output=%s\\mosaic_blend.png', imfolder) ;
    end
    if strcmp(projection_type,'equi')
        if abs(mosaicw / fe - pi) < 0.2
            sysorder = sprintf('%s --wrap=horizontal', sysorder) ;
        end
    end
    for i = 1:im_n
        sysorder = sprintf('%s %s\\p_%02d.png',sysorder, imfolder, i) ;
    end
    sysorder = sprintf('%s --gpu', sysorder) ;
    system(sysorder) ;
    
    for i = 1:im_n
        imp_path = sprintf('%s\\p_%02d.png', imfolder, i);
        delete(imp_path);
    end
end

end

