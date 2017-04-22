function [u, v] = trans_cam2equi(xcam, ycam, zcam, fe)

rr = sqrt(xcam .^ 2 + ycam .^ 2 + zcam .^ 2) ;

theta = acos(zcam ./ sqrt(xcam .^ 2 + zcam .^ 2)) .* (2 * (xcam > 0) - 1) ;
phi = asin(ycam ./ rr);

u = fe * theta ;
v = fe * phi ;

end

