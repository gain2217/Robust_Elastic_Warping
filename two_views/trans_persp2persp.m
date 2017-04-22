function [x, y] = trans_persp2persp(u, v, R, M1, D1, M2, D2)

[xcam, ycam, zcam] = trans_persp2cam(u, v, M1, D1);

xr = R(1, 1) * xcam + R(1, 2) * ycam + R(1, 3) * zcam ;
yr = R(2, 1) * xcam + R(2, 2) * ycam + R(2, 3) * zcam ;
zr = R(3, 1) * xcam + R(3, 2) * ycam + R(3, 3) * zcam ;

[x, y] = trans_cam2persp(xr, yr, zr, M2, D2);

end