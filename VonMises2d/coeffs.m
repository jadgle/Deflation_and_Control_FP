function u = coeffs(X,Y,wx,wy,L,rho0)

a = zeros(1 + 4*L + 4*L^2,1);
k = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 0
a(k)   = 2*sum(rho0);

k = 2;
for m = 1:L
        sin_m = sin(m*Y');
        cos_m = cos(m*Y');

        a(k) = 2*sin_m*(rho0.*wx.*wy);
        a(k+1) = 2*cos_m*(rho0.*wx.*wy);
        k = k+2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:L
    cos_n = cos(n*X');
    sin_n = sin(n*X');
    a(k)   = 2*sin_n*(rho0.*wx.*wy);
    a(k+1) = 2*cos_n*(rho0.*wx.*wy);
    k = k+2;
    for m = 1:L
        cos_m = cos(m*Y');
        sin_m = sin(m*Y');
        a(k)   = 4*sin_n.*sin_m*(rho0.*wx.*wy);
        a(k+1) = 4*sin_n.*cos_m*(rho0.*wx.*wy);
        a(k+2) = 4*cos_n.*sin_m*(rho0.*wx.*wy);
        a(k+3) = 4*cos_n.*cos_m*(rho0.*wx.*wy);
        k = k+4;
    end
end

u = 1/(4*pi^2)*a;