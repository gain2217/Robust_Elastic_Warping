function [xcam, ycam, zcam] = trans_equi2cam(u, v, fe)
theta = u ./ fe ;
phi = v ./ fe;

ycam = sin(phi) ;
zcam = cos(phi) .* cos(theta) ;
xcam = cos(phi) .* sin(theta) ;

end