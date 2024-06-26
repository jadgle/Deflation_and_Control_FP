function result = R_complete(m1,m2,beta,alpha,kappa)
    
    V_hkb = @(x) alpha*cos(2*x);
    Z = @(x,m1,m2) exp(-beta*V_hkb(x)+(beta*kappa*(m1*cos(x)+m2*sin(x))));
    rho = @(x,m1,m2,Z) (1/Z)*exp(-beta*V_hkb(x)+(beta*kappa*(m1*cos(x)+m2*sin(x))));
    
    n_gauss = 100;
    [xi,wi]=Gauss_quad(n_gauss,0,2*pi);

    Z_m = sum(Z(xi,m1,m2).*wi);
    R_m1 = sum(rho(xi,m1,m2,Z_m).*cos(xi).*wi);
    R_m2 = sum(rho(xi,m1,m2,Z_m).*sin(xi).*wi);
    result = [R_m1,R_m2];   
end