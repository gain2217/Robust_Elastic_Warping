function [u, v] = trans_cam2persp(xcam, ycam, zcam, M, D)
fx = M(1,1); fy = M(2,2);
cx = M(1,3); cy = M(2,3);
if D ~= 0
    k1 = D(1); k2 = D(2); p1 = D(3); p2 = D(4);
    k3 = D(5); k4 = D(6); k5 = D(7); k6 = D(8);
end

x = xcam./zcam;
y = ycam./zcam;
% zcam should > 0, otherwise some unexpected projection errors may occur
x(zcam<=0) = 100*cx;
y(zcam<=0) = 100*cy;

if D ~= 0
    r2 = x.^2 + y.^2;
    k_radial = (1+k1*r2+k2*r2.^2+k3*r2.^3)./(1+k4*r2+k5*r2.^2+k6*r2.^3);
    delta_x = 2*p1*x.*y+p2*(r2+2*x.^2);
    delta_y = p1*(r2+2*y.^2)+2*p2*x.*y;
    x_d = x.*k_radial + delta_x;
    y_d = y.*k_radial + delta_y;
else
    x_d = x;
    y_d = y;
end

u = fx * x_d + cx;
v = fy * y_d + cy;

end

