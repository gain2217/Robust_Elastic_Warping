function [u, v] = trans_persp2equi(x, y, R, M, D, fe)

[xcam, ycam, zcam] = trans_persp2cam(x, y, M, D) ;

xr = R(1, 1) * xcam + R(1, 2) * ycam + R(1, 3) * zcam ;
yr = R(2, 1) * xcam + R(2, 2) * ycam + R(2, 3) * zcam ;
zr = R(3, 1) * xcam + R(3, 2) * ycam + R(3, 3) * zcam ;

[u, v] = trans_cam2equi(xr, yr, zr, fe);

end