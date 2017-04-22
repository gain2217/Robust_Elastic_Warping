function [ err ] = residual_KR_robust( X1, X2, imsize1, imsize2, paras, sigma )

% parameretes sigma indicates the distance scope for inliers with homopraphy matrix, in pixels

k1 = paras(1);
k2 = paras(2);
theta = paras(3:5);
% yaw = paras(3);
% pitch = paras(4);
% roll = paras(5);

K1 = [k1, 0, imsize1(2)/2;
     0, k1, imsize1(1)/2;
     0,  0, 1];
K2 = [k2, 0, imsize2(2)/2;
      0, k2, imsize2(1)/2;
      0,  0, 1];
theta_m = [0         -theta(3) theta(2)
           theta(3)  0         -theta(1)
           -theta(2) theta(1)  0];
R = expm(theta_m);
% Ry = [cos(yaw),     0,              -sin(yaw);
%       0,            1,              0;
%       sin(yaw),     0,              cos(yaw)] ;
% Rp = [1,            0               0;
%       0,            cos(pitch),     sin(pitch);
%       0,            -sin(pitch),    cos(pitch)] ;
% Rr = [cos(roll),    sin(roll),      0;
%       -sin(roll),   cos(roll),      0;
%       0,            0,              1] ;
% R = Rr * Rp * Ry;

err = residual_H(X1, X2, K1, K2, R);

outlier = (abs(err) > sigma);
err(outlier) = sign(err(outlier)) .* (sigma + sigma * log(abs(err(outlier))/sigma));
% err(outlier) = sign(err(outlier)) .* sqrt(2*sigma*abs(err(outlier)) - sigma*sigma);

end

