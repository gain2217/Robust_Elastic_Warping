function [ err ] = residual_all_robust( X, imsize, edge_list, paras, sigma )

% parameretes sigma indicates the distance scope for inliers with homopraphy matrix, in pixels

im_n = size(imsize, 1);
edge_n = size(X, 1);

R_pair = cell(im_n, im_n);
for i = 1 : im_n
    R_pair{i, i} = eye(3);
end
for ii = 2:im_n
    theta = paras(im_n+3*(ii-2)+1:im_n+3*(ii-2)+3);
    theta_m = [0         -theta(3) theta(2)
               theta(3)  0         -theta(1)
               -theta(2) theta(1)  0];
    R_pair{1,ii} = expm(theta_m);
    R_pair{ii,1} = R_pair{1,ii}';
end
for i = 2:im_n-1
    for j = i+1:im_n
        R_pair{i,j} = R_pair{1,j}*R_pair{i,1};
        R_pair{j,i} = R_pair{i,j}';
    end
end

err = [];
for ei = 1 : edge_n
    i = edge_list(ei, 1);
    j = edge_list(ei, 2);
    
    ki = paras(i);
    kj = paras(j);
    Ki = [ki, 0, imsize(i,2)/2;
          0, ki, imsize(i,1)/2;
          0,  0, 1];
    Kj = [kj, 0, imsize(i,2)/2;
          0, kj, imsize(j,1)/2;
          0,  0, 1];
    
    err_ij = residual_KR( double(X{ei,1}), double(X{ei,2}), Ki, Kj, R_pair{i,j});
    err_ji = residual_KR( double(X{ei,2}), double(X{ei,1}), Kj, Ki, R_pair{j,i});
    err = cat(1,err,err_ij,err_ji);
end

outlier = (abs(err) > sigma);
err(outlier) = sign(err(outlier)) .* (sigma + sigma * log(abs(err(outlier))/sigma));
% err(outlier) = sign(err(outlier)) .* sqrt(2*sigma*abs(err(outlier)) - sigma*sigma);

end

