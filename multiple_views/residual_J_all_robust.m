function [ err, J] = residual_J_all_robust( X, imsize, edge_list, paras, sigma )

% parameretes sigma indicates the distance scope for inliers with homopraphy matrix, in pixels

cal_Jacobian = false;

im_n = size(imsize, 1);
edge_n = size(X, 1);

M = cell(im_n, 1);
Minv = cell(im_n, 1);
R = cell(im_n, 1);
if cal_Jacobian
    DM_Dk = [1 0 0; 0 1 0; 0 0 0]; % for Jacobian
    DMinv_Dk = cell(im_n, 1); % for Jacobian
    DR_Dthetai1 = cell(im_n, 1); % for Jacobian
    DR_Dthetai2 = cell(im_n, 1); % for Jacobian
    DR_Dthetai3 = cell(im_n, 1); % for Jacobian
    DRinv_Dthetai1 = cell(im_n, 1); % for Jacobian
    DRinv_Dthetai2 = cell(im_n, 1); % for Jacobian
    DRinv_Dthetai3 = cell(im_n, 1); % for Jacobian
end

for i = 1 : im_n
    ki = paras(i);
    cx = imsize(i,2)/2;
    cy = imsize(i,1)/2;
    M{i} = [ki, 0, cx;
            0,  ki,cy;
            0,  0, 1];
    Minv{i} = [1/ki, 0,   -cx/ki;
               0     1/ki -cy/ki;
               0     0    1];
    
    if cal_Jacobian
        %% for Jacobian
        DMinv_Dk{i} = [-1/ki^2 0       cx/ki^2;
                       0       -1/ki^2 cy/ki^2;
                       0       0       0];
    end
end
R{1} = eye(3);
for i = 2:im_n
    theta = paras(im_n+3*(i-2)+1:im_n+3*(i-2)+3);
    theta_m = [0         -theta(3) theta(2)
               theta(3)  0         -theta(1)
               -theta(2) theta(1)  0];
    R{i} = expm(theta_m);
    
    if cal_Jacobian
        %% for Jacobian
        DR_Dthetai1{i} = R{i} * [0 0 0; 0 0 -1; 0 1 0];
        DR_Dthetai2{i} = R{i} * [0 0 1; 0 0 0; -1 0 0];
        DR_Dthetai3{i} = R{i} * [0 -1 0; 1 0 0; 0 0 0];
        DRinv_Dthetai1{i} = -R{i}' * [0 0 0; 0 0 -1; 0 1 0];
        DRinv_Dthetai2{i} = -R{i}' * [0 0 1; 0 0 0; -1 0 0];
        DRinv_Dthetai3{i} = -R{i}' * [0 -1 0; 1 0 0; 0 0 0];
    end
end

err = [];
J = []; % Jacobian
p_n = size(edge_n, 1);
for ei = 1 : edge_n
    p_n(ei) = size(X{ei,1},2);
    for si = 1:2
        if si == 1
            i = edge_list(ei, 1);
            j = edge_list(ei, 2);
            xi = double(X{ei,1}(1,:));
            yi = double(X{ei,1}(2,:));
            xj = double(X{ei,2}(1,:));
            yj = double(X{ei,2}(2,:));
        else
            i = edge_list(ei, 2);
            j = edge_list(ei, 1);
            xi = double(X{ei,2}(1,:));
            yi = double(X{ei,2}(2,:));
            xj = double(X{ei,1}(1,:));
            yj = double(X{ei,1}(2,:));
        end
        
        H = M{j}*R{j}*R{i}'*Minv{i};
        
        z = H(3,1)*xi + H(3,2)*yi + H(3,3);
        xi_j = (H(1,1)*xi + H(1,2)*yi + H(1,3)) ./ z;
        yi_j = (H(2,1)*xi + H(2,2)*yi + H(2,3)) ./ z;
        err_ij = [xi_j - xj; yi_j - yj];
        
        err = cat(1,err,err_ij(1,:)',err_ij(2,:)');
        
        if cal_Jacobian
           %% for Jacobian
            J_ei = zeros(2*p_n(ei), 4*im_n-3);
            
            Xi = [xi;yi;zeros(1,p_n(ei))];
            
            % about instrinsic parameters
            Dph_Dkj = DM_Dk * R{j} * R{i}' * Minv{i} * Xi;
            Dph_Dki = M{j} * R{j} * R{i}' * DMinv_Dk{i} * Xi;
            
            Dx_Dph = [(1./z)', zeros(p_n(ei),1), -(xi_j./z.^2)'];
            Dy_Dph = [zeros(p_n(ei),1), (1./z)', -(yi_j./z.^2)'];
            
            Dx_Dkj = sum(Dx_Dph .* Dph_Dkj',2);
            Dx_Dki = sum(Dx_Dph .* Dph_Dki',2);
            Dy_Dkj = sum(Dy_Dph .* Dph_Dkj',2);
            Dy_Dki = sum(Dy_Dph .* Dph_Dki',2);
            
            J_ei(        1:  p_n(ei),j) = Dx_Dkj;
            J_ei(p_n(ei)+1:2*p_n(ei),j) = Dy_Dkj;
            J_ei(        1:  p_n(ei),i) = Dx_Dki;
            J_ei(p_n(ei)+1:2*p_n(ei),i) = Dy_Dki;
            
            % about extrinsic parameters
            if j > 1
                Dph_Dtheta1_j = M{j} * DR_Dthetai1{j} * R{i}' * Minv{i} * Xi;
                Dph_Dtheta2_j = M{j} * DR_Dthetai2{j} * R{i}' * Minv{i} * Xi;
                Dph_Dtheta3_j = M{j} * DR_Dthetai3{j} * R{i}' * Minv{i} * Xi;
                
                Dx_Dtheta1_j = sum(Dx_Dph .* Dph_Dtheta1_j',2);
                Dx_Dtheta2_j = sum(Dx_Dph .* Dph_Dtheta2_j',2);
                Dx_Dtheta3_j = sum(Dx_Dph .* Dph_Dtheta3_j',2);
                
                Dy_Dtheta1_j = sum(Dy_Dph .* Dph_Dtheta1_j',2);
                Dy_Dtheta2_j = sum(Dy_Dph .* Dph_Dtheta2_j',2);
                Dy_Dtheta3_j = sum(Dy_Dph .* Dph_Dtheta3_j',2);
                
                J_ei(        1:  p_n(ei),im_n+3*(j-2)+1) = Dx_Dtheta1_j;
                J_ei(p_n(ei)+1:2*p_n(ei),im_n+3*(j-2)+1) = Dy_Dtheta1_j;
                J_ei(        1:  p_n(ei),im_n+3*(j-2)+2) = Dx_Dtheta2_j;
                J_ei(p_n(ei)+1:2*p_n(ei),im_n+3*(j-2)+2) = Dy_Dtheta2_j;
                J_ei(        1:  p_n(ei),im_n+3*(j-2)+3) = Dx_Dtheta3_j;
                J_ei(p_n(ei)+1:2*p_n(ei),im_n+3*(j-2)+3) = Dy_Dtheta3_j;
                
            end
            if i > 1
                Dph_Dtheta1_i = M{j} * R{j} * DRinv_Dthetai1{i} * Minv{i} * Xi;
                Dph_Dtheta2_i = M{j} * R{j} * DRinv_Dthetai2{i} * Minv{i} * Xi;
                Dph_Dtheta3_i = M{j} * R{j} * DRinv_Dthetai3{i} * Minv{i} * Xi;
                
                Dx_Dtheta1_i = sum(Dx_Dph .* Dph_Dtheta1_i',2);
                Dx_Dtheta2_i = sum(Dx_Dph .* Dph_Dtheta2_i',2);
                Dx_Dtheta3_i = sum(Dx_Dph .* Dph_Dtheta3_i',2);
                
                Dy_Dtheta1_i = sum(Dy_Dph .* Dph_Dtheta1_i',2);
                Dy_Dtheta2_i = sum(Dy_Dph .* Dph_Dtheta2_i',2);
                Dy_Dtheta3_i = sum(Dy_Dph .* Dph_Dtheta3_i',2);
                
                J_ei(        1:  p_n(ei),im_n+3*(i-2)+1) = Dx_Dtheta1_i;
                J_ei(p_n(ei)+1:2*p_n(ei),im_n+3*(i-2)+1) = Dy_Dtheta1_i;
                J_ei(        1:  p_n(ei),im_n+3*(i-2)+2) = Dx_Dtheta2_i;
                J_ei(p_n(ei)+1:2*p_n(ei),im_n+3*(i-2)+2) = Dy_Dtheta2_i;
                J_ei(        1:  p_n(ei),im_n+3*(i-2)+3) = Dx_Dtheta3_i;
                J_ei(p_n(ei)+1:2*p_n(ei),im_n+3*(i-2)+3) = Dy_Dtheta3_i;
            end
            
            J = cat(1,J,J_ei);
        end
    end
    
end

outlier = (abs(err) > sigma);
if any(outlier)
    if cal_Jacobian
        J(outlier,:) = J(outlier,:) .* (sigma./err(outlier)*ones(1,4*im_n-3));
    end
    err(outlier) = sign(err(outlier)) .* (sigma + sigma * log(abs(err(outlier))/sigma));
    % err(outlier) = sign(err(outlier)) .* sqrt(2*sigma*abs(err(outlier)) - sigma*sigma);
end

end

