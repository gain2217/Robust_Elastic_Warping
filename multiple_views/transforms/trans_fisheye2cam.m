function [xcam, ycam, zcam] = trans_fisheye2cam(x, y, M, D)
fx = M(1,1); fy = M(2,2);
cx = M(1,3); cy = M(2,3);
if D ~= 0
    k1 = D(1); k2 = D(2); k3 = D(3); k4 = D(4);
end

alpha_d_x = (x - cx) ./ fx;
alpha_d_y = (y - cy) ./ fy;
alpha_d = sqrt(alpha_d_x .^ 2 + alpha_d_y .^ 2) ;

if D ~= 0
    % Iterative solver
    alpha = alpha_d;
    for kk = 1:20
        alpha_2 = alpha .^ 2;
        k_radial =  1 + k1 * alpha_2 + k2 * alpha_2.^2 + k3 * alpha_2.^3 + k4 * alpha_2.^4;
        alpha = alpha_d ./ k_radial;
    end
else
    alpha = alpha_d;
end

zcam = cos(alpha) ;
xcam = sin(alpha) .* alpha_d_x ./ alpha_d ;
ycam = sin(alpha) .* alpha_d_y ./ alpha_d ;

end