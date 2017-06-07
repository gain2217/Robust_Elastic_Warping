function [x, y] = trans_equi2fisheye(u, v, R, M, D, fe)

[xcam, ycam, zcam] = trans_equi2cam(u, v, fe);

xr = R(1, 1) * xcam + R(1, 2) * ycam + R(1, 3) * zcam ;
yr = R(2, 1) * xcam + R(2, 2) * ycam + R(2, 3) * zcam ;
zr = R(3, 1) * xcam + R(3, 2) * ycam + R(3, 3) * zcam ;

[x, y] = trans_cam2fisheye(xr, yr, zr, M, D);

end