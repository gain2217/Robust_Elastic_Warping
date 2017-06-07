function Xcam = trans_fisheye2cam_X(X, M, D)

x = X(1, :) ./ X(3, :) ;
y = X(2, :) ./ X(3, :) ;

[xcam, ycam, zcam] = trans_fisheye2cam(x, y, M, D);

Xcam = [xcam; ycam; zcam] ;

end