% Coefficients
theta    =   1;
beta_m1  =   1e-1;

%% Discretization

% prepping for Gauss quadratures
n_gauss = 200;
[xi,wi]=Gauss_quad(n_gauss,-pi,pi);
[yj,zj]=Gauss_quad(n_gauss,-pi,pi);
[Xi, Yj] = meshgrid(xi,yj);
[wi, zj] = meshgrid(wi,zj);
xi = reshape(Xi,(n_gauss)^2,1);
yj = reshape(Yj,(n_gauss)^2,1);
wi = reshape(wi,(n_gauss)^2,1);
zj = reshape(zj,(n_gauss)^2,1);

% potentials
I02 = (besselj(0, 1i))^2;

h = @(x,y) 0; % confining force
W = @(x,y) -exp(theta*(cos(x)+cos(y)))/I02; % interaction potential

l = 10; %number of modes
% 12. 15. 20



%% Discretization
% basis functions multi index (for n and m)
multi_index = zeros(1 + 4*l + 4*l^2,2);
psi = cell(1 + 4*l + 4*l^2,1);
dpsi = cell(1 + 4*l + 4*l^2,1);

psi{1}  = @(x,y) ones(size(x));
dpsi{1} = @(x,y) [zeros(size(x)),zeros(size(x))];

i = 2;
for m = 1:l
        multi_index(i:i+1,:)   =  repmat([0,m],2,1);
        psi{i}   = @(x,y) sin(m*y);
        psi{i+1} = @(x,y) cos(m*y);
        dpsi{i}   = @(x,y) [zeros(size(x)),m*cos(m*y)];
        dpsi{i+1} = @(x,y) [zeros(size(x)),-m*sin(m*y)];
        i = i+2;
end
for n = 1:l
    multi_index(i:i+1,:)   =  repmat([n,0],2,1);

    psi{i}   = @(x,y) sin(n*x);
    psi{i+1} = @(x,y) cos(n*x);
    dpsi{i}   = @(x,y) [n*cos(n*x),zeros(size(x))];
    dpsi{i+1} = @(x,y) [-n*sin(n*x),zeros(size(x))];

    i = i+2;
    for m = 1:l
        multi_index(i:i+3,:)   =  repmat([n,m],4,1);

        psi{i}    = @(x,y) sin(n*x).*sin(m*y);
        dpsi{i}   = @(x,y) [n*cos(n*x).*sin(m*y),m*cos(m*y).*sin(n*x)];
        psi{i+1}  = @(x,y) sin(n*x).*cos(m*y);
        dpsi{i+1} = @(x,y) [n*cos(n*x).*cos(m*y),-m*sin(m*y).*sin(n*x)];
        psi{i+2}  = @(x,y) cos(n*x).*sin(m*y);
        dpsi{i+2} = @(x,y) [-n*sin(n*x).*sin(m*y),m*cos(m*y).*cos(n*x)];
        psi{i+3}  = @(x,y) cos(n*x).*cos(m*y);
        dpsi{i+3} = @(x,y) [-n*sin(n*x).*cos(m*y),-m*sin(m*y).*cos(n*x)];

        i = i+4;
    end
end

%%
% Generating the system matrices
L = 1 + 4*l + 4*l^2; 
D1 = zeros(L,L,L);
D2 = D1;
parfor k = 1:L
    f_k = psi{k}(xi,yj);
    for m = 1:L
        df_m = dpsi{m}(xi,yj); 
        for n = 1:L
            f_n  = psi{n}(xi,yj); 
            D1(n,k,m) = sum(f_k.*df_m(:,1).*f_n.*wi.*zj);
            D2(n,k,m) = sum(f_k.*df_m(:,2).*f_n.*wi.*zj);
        end
    end
end



writematrix(D1, 'D1_vonmises')
writematrix(D2, 'D2_vonmises')