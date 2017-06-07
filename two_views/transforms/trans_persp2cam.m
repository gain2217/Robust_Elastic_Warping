function [xcam, ycam, zcam] = trans_persp2cam(u, v, M, D)
fx = M(1,1); fy = M(2,2);
cx = M(1,3); cy = M(2,3);
if D ~= 0
    k1 = D(1); k2 = D(2); p1 = D(3); p2 = D(4);
    k3 = D(5); k4 = D(6); k5 = D(7); k6 = D(8);
end

x_d = (u - cx) / fx;
y_d = (v - cy) / fy;

if D ~= 0
    % Iterative solver
    x = x_d;
    y = y_d;
    for kk = 1:20
        r2 = x.^2 + y.^2;
        k_radial = (1+k1*r2+k2*r2.^2+k3*r2.^3)./(1+k4*r2+k5*r2.^2+k6*r2.^3);
        delta_x = 2*p1*x.*y+p2*(r2+2*x.^2);
        delta_y = p1*(r2+2*y.^2)+2*p2*x.*y;
        x = (x_d - delta_x) ./ k_radial;
        y = (x_d - delta_y) ./ k_radial;
    end
else
    x = x_d;
    y = y_d;
end

xcam = x;
ycam = y;
zcam = ones(size(xcam));

end

