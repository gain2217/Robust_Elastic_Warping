function [x, y] = trans_cam2fisheye(xcam, ycam, zcam, M, D)
fx = M(1,1); fy = M(2,2);
cx = M(1,3); cy = M(2,3);
if D ~= 0
    k1 = D(1); k2 = D(2); k3 = D(3); k4 = D(4);
end

rcam = sqrt(xcam .^ 2 + ycam .^ 2 + zcam .^ 2) ;
alpha = acos(zcam ./ rcam) ;

if D ~= 0
    alpha_2 = alpha .^ 2;
    k_radial =  1 + k1 * alpha_2 + k2 * alpha_2.^2 + k3 * alpha_2.^3 + k4 * alpha_2.^4;
    alpha_d = alpha .* k_radial;
else
    alpha_d = alpha;
end

rcam_xy = sqrt(xcam.^2 + ycam.^2);
alpha_d_x = xcam ./ rcam_xy .* alpha_d;
alpha_d_y = ycam ./ rcam_xy .* alpha_d;

x = alpha_d_x .* fx + cx ;
y = alpha_d_y .* fy + cy ;

end