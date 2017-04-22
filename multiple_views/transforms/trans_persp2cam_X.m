function Xcam = trans_persp2cam_X(X, M, D)

u = X(1, :) ./ X(3, :) ;
v = X(2, :) ./ X(3, :) ;

[xcam, ycam, zcam] = trans_persp2cam(u, v, M, D);

Xcam = [xcam; ycam; zcam] ;

end

