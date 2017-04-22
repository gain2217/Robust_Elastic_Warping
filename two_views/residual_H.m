function [ errH ] = residual_H( X1, X2, K1, K2, R)

% homopraphy matrix
H = K2*R/K1;

X1_p2 = H * X1;
X1_p2(1,:) = X1_p2(1,:) ./ X1_p2(3,:) ;
X1_p2(2,:) = X1_p2(2,:) ./ X1_p2(3,:) ;
errH1_2 = X2 - X1_p2;
X2_p1 = H \ X2;
X2_p1(1,:) = X2_p1(1,:) ./ X2_p1(3,:) ;
X2_p1(2,:) = X2_p1(2,:) ./ X2_p1(3,:) ;
errH2_1 = X1 - X2_p1;

errH = [errH1_2(1,:)', errH1_2(2,:)', errH2_1(1,:)', errH2_1(2,:)'];

end

