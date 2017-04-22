function [ H, ok, score ] = HM_ransac(X1, X2, Nr, min_dis)

N = size(X1,2);

u = X1(1,:)';
v = X1(2,:)';
u_ = X2(1,:)';
v_ = X2(2,:)';

% normalization (isotropic scaling)
% Reference: In Defense of the Eight-Point Algorithm
scale = 1/mean([u;u_;v;v_]);
u = u * scale;
v = v * scale;
u_ = u_ * scale;
v_ = v_ * scale;

% normalization (non-isotropic scaling)
% Reference: In Defense of the Eight-Point Algorithm
% Us = [sum(u.^2), sum(u.*v), sum(u);
%     sum(v.*u), sum(v.^2), sum(v);
%     sum(u), sum(v), N;];
% K = chol(Us/N);
% X1n = K' \ X1;
% 
% Us_ = [sum(u_.^2), sum(u_.*v_), sum(u_);
%     sum(v_.*u_), sum(v_.^2), sum(v_);
%     sum(u_), sum(v_), N;];
% K_ = chol(Us_/N);
% X2n = K_' \ X2;

% u = X1n(1,:)';
% v = X1n(2,:)';
% w = X1n(3,:)';
% u_ = X2n(1,:)';
% v_ = X2n(2,:)';
% w_ = X2n(3,:)';

% coefficient matrix
% A1 = [zeros(N,3),          -w_.*u, -w_.*v, -w_.*w, v_.*u, v_.*v, v_.*w];
% A2 = [w_.*u, w_.*v, w_.*w, zeros(N,3),             -u_.*u, -u_.*v, -u_.*w];
% % A3 = [-v_.*u,-v_.*v ,-v_.*w , u_*u, u_*v, u_*w, zeros(N,3)];

% % coefficient matrix
A1 = [zeros(N,3),      -u, -v, -ones(N,1), v_.*u, v_.*v, v_];
A2 = [u, v, ones(N,1), zeros(N,3),         -u_.*u, -u_.*v, -u_];

if (min_dis > 0)
    H = cell(Nr,1);
    ok = cell(Nr,1);
    score = zeros(Nr,1);
    for t = 1:Nr
        % estimate foundamental matrix
%         subset = vl_colsubset(1:size(X1,2), 4) ;
        subset = randsample(N,4);
        A = [A1(subset,:);A2(subset,:)];
        [U,S,V] = svd(A);
        h = V(:,9);
        H{t} = reshape(h,3,3)';

        % score foundamental matrix
        dis2 = (A1 * h).^2 + (A2 * h).^2;
        ok{t} = dis2 < min_dis * min_dis;
        score(t) = sum(ok{t}) ;
        
    end
    [score, best] = max(score) ;
    ok = ok{best} ;
    A = [A1(ok,:);A2(ok,:)];
    [U,S,V] = svd(A,'econ');
    h = V(:,9);
    H = reshape(h,3,3)';
else
    A = [A1;A2];
    [U,S,V] = svd(A,'econ');
    h = V(:,9);
    H = reshape(h,3,3)';
end

% denormalization (isotropic scaling)
H = [1/scale, 0, 0; 0, 1/scale, 0; 0, 0, 1] * H * [scale, 0, 0; 0, scale, 0; 0, 0, 1];

% denormalization (non-isotropic scaling)
% H = K_' * H / K';

end
