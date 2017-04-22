[ H, ok, score ] = HM_ransac(X1, X2, 500, 0.1);

figure; clf;
marker_size = 15;
imshow([im1,im2], 'border', 'tight'); hold on;
plot(X1(1,ok), X1(2,ok), 'b.', 'MarkerSize', marker_size);
plot(X2(1,ok)+imsize1(2), X2(2,ok), 'b.', 'MarkerSize', marker_size);
plot(X1(1,~ok), X1(2,~ok), 'r.', 'MarkerSize', marker_size);
plot(X2(1,~ok)+imsize1(2), X2(2,~ok), 'r.', 'MarkerSize', marker_size);

% h = line([X1(1,ok) ; X2(1,ok)+imsize1(2)], [X1(2,ok) ; X2(2,ok)], 'LineWidth', 2) ;

X1_ok = X1(:,ok);
X2_ok = X2(:,ok);

% bundle adjustment on the normallized matching data
options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'Display','final',...
    'MaxFunEvals',2000, 'MaxIter',1e3, 'TolFun',1e-6, 'TolX',1e-6, 'Jacobian','off');
[paras_init,resnorm_init,residual_init,exitflag_init,output_init] = lsqnonlin(...
    @(p)residual_KR_robust(double(X1_ok), double(X2_ok), imsize1, imsize2, p, 1000), [1000 1000 0 0 0],...
    [],[],options);
[paras,resnorm,residual,exitflag,output] = lsqnonlin(...
    @(p)residual_KR_robust(double(X1_ok), double(X2_ok), imsize1, imsize2, p, 10), paras_init,...
    [],[],options);
k1 = paras(1);
k2 = paras(2);
theta = paras(3:5);

K1 = [k1, 0, imsize1(2)/2;
     0, k1, imsize1(1)/2;
     0,  0, 1];
K2 = [k2, 0, imsize2(2)/2;
      0, k2, imsize2(1)/2;
      0,  0, 1];
theta_m = [0         -theta(3) theta(2)
           theta(3)  0         -theta(1)
           -theta(2) theta(1)  0];
R = expm(theta_m);

H = K2*R/K1;
M1 = K1; M2 = K2; D1 = 0; D2 = 0;