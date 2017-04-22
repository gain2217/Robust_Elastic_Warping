function [ err ] = residual_KR( X1, X2, K1, K2, R)

% homopraphy matrix
H = K2*R/K1;
n = size(X1,2);

x1 = X1(1,:); y1 = X1(2,:);
x2 = X2(1,:); y2 = X2(2,:);

w = H(3,1)*x1 + H(3,2)*y1 + H(3,3);
x1_2 = (H(1,1)*x1 + H(1,2)*y1 + H(1,3)) ./ w;
y1_2 = (H(2,1)*x1 + H(2,2)*y1 + H(2,3)) ./ w;
err1_2 = [x1_2 - x2; y1_2 - y2];

err = reshape(err1_2,[2*n,1]);

end

