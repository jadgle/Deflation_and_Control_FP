%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                    %
%%                     Solving Von Mises 2d steady states with Newton (and deflation) method                        %% 
%                                                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global L D1 D2 M1 B nt a_infty lambda

% Coefficients
theta    =   1;
kappa    =   1;
beta_m1  =   0.3701; %phase transition around 0.3701 
beta = 1/beta_m1;

%% Discretization


% potential
W = @(x,y) -exp(theta*(cos(x)+cos(y)))/(besselj(0, 1i))^2; % interaction potential with kappa=1

l = 10; %number of modes




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

%% Generating/importing system matrices

L = 1 + 4*l + 4*l^2; 
A1 = zeros(L,1);
M = zeros(L,1);
C1 = zeros(L);

% prepping for Gauss quadratures
n_gauss = 100;
[xi,wi]=Gauss_quad(n_gauss,-pi,pi);
[yj,zj]=Gauss_quad(n_gauss,-pi,pi);
[Xi, Yj] = meshgrid(xi,yj);
[wi, zj] = meshgrid(wi,zj);
xi = reshape(Xi,(n_gauss)^2,1);
yj = reshape(Yj,(n_gauss)^2,1);
wi = reshape(wi,(n_gauss)^2,1);
zj = reshape(zj,(n_gauss)^2,1);


for i = 1:L
    df_i = dpsi{i};
    f_i  = psi{i};
    A1(i)       =  sum(sum(df_i(xi,yj).^2,2).*wi.*zj);
    M(i)        =  sum(f_i(xi,yj).^2.*wi.*zj);
end

A2 = diag(beta_m1*A1);
M1 = diag(M.^-1);

tic
B1 = readmatrix("B_vonmises.txt");
D1 = readmatrix("D1_vonmises.txt");
D2 = readmatrix("D2_vonmises.txt");
toc

B1 = reshape(B1,L,L,L); %bilinear form
D1 = reshape(D1,L,L,L); %control form1
D2 = reshape(D2,L,L,L); %control form2

B = kappa*B1;


%% Simulation for the initial guess - iterative scheme Vanden-Eijnden 2022 
rng('default')  % For reproducibility

basis = zeros((n_gauss)^2,L);
for i = 1:L
    f_i = psi{i};
    basis(:,i) = f_i(xi,yj);
end
beta_m1 = 0.3701;
beta = 1/beta_m1;
guess0 =rand(size(xi));
guess0 = guess0/(sum(guess0.*wi.*zj));

convegence_flag = 1;
tol = 1e-2;
max_iteration = 5;
iteration = 1;
f_old = guess0;
f_new = guess0*0;
temp = f_new;
%%
iteration=1;
while iteration < max_iteration
    for i = 1:n_gauss^2
            temp = W(xi(i)-xi,yj(i)-yj).*f_old;
            f_new(i) = exp(-beta*sum(temp.*wi.*zj));
    end
    f_new = f_new/(sum(f_new.*wi.*zj));
    diff = norm(f_old-f_new);
    if diff<tol
        break
    end
    f_old = f_new;
    iteration = iteration + 1;
end

%%
Z = reshape(f_new,n_gauss,n_gauss);

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'DefaultLegendFontSize',24)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',26)
set(0,'DefaultLineLineWidth',1.2);
set(gcf,'color','w');
figure(2)
surf(Xi,Yj,Z)
title('Simulated guess for Newton')

guess = coeffs(xi,yj,wi,zj,l,f_new);


%% Newton with deflation

A = [A2;zeros(1,L)]; % integral=1 condition
C = [C1;zeros(1,L)]; % integral=1 condition
for i = 1:L
    f_i = psi{i};
    A(end,i) = sum(f_i(xi,yj).*wi.*zj);
end
update = 100;        % turning off the convergence flag
tolerance = 1e-10;    % tolerance for convergence
max_iteration = 50; % number of before assuming failed convergence

n_deflation = 0; % counter for the number of solutions
iteration = 0;   % counter for the number of iterations

solutions = [];  % list of solutions

% parameters and start for deflation
p = 1;       % deflation power
shift = .01; % deflation shift
eta0 = @(u) 1; 
new_eta = @(u) 1;
deta0 = @(u) zeros(size(u));
new_deta = @(u) zeros(size(u));

trying = 1;
while trying < 2 % this is trying>2 is useful when we are changing initial condition
    u = guess;

    while (update>tolerance && iteration < max_iteration)
        iteration = iteration + 1;
        
        % undeflated system at u
        f  = F(u,A,B,C);
        df = dF(u,A,B,C);
        
        % updating the deflation operator
        eta_now  = new_eta(u);
        deta_now = new_deta(u);
        eta_pre  = eta0(u);
        deta_pre = deta0(u);
        
        % deflated system at u
        G  = f./(eta_pre*eta_now)+shift*f;
        dG = df./(eta_pre*eta_now) - (f./(eta_pre.^2*eta_now.^2))*(deta_pre*eta_now + eta_pre*deta_now)' + shift*df;
        
        % update for u 
        u_new  = u - pinv(dG)*G;
        update = norm(u_new-u);
        u      = u_new;
    end
    
    if (iteration == max_iteration || sum(isnan(u))>1) % if the previous loop has finished without convergence, we start
        trying = trying + 1;                         % again from a different initial guess
        guess = zeros(size(u));   % exploring different initial guess
        %guess = u_unif; % checking if the uniform distribution is a solution    

    else   % if the previous loop has finished with convergence, we deflate w.r.t. the reached root
        % updating the deflation eta
        deta0 =  @(r) new_deta(r)*eta0(r)+new_eta(r)*deta0(r);
        eta0  =  @(r) new_eta(r)*eta0(r);
        new_eta  =  @(r) norm(r-u_new).^p;
        new_deta =  @(r) (p*norm(r-u_new)^(p-2))*(r-u_new); % problem with dimension in dG, see handwritten notes 
        
        sprintf('Solution for the %d-ith deflation found after %d iterations',n_deflation,iteration)
        n_deflation = n_deflation + 1;
        solutions = [solutions, u];
        
        update    = 100;  % reset of the convergence flag 
        iteration = 0;    % reset of the divergence  flag
        guess = zeros(size(u)); % checking if the uniform distribution is a solution   
    end
     
end


%% Plotting the results
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex')
set(0,'defaultAxesXGrid','on')
set(0,'defaultAxesYGrid','on')
set(0,'DefaultLegendFontSize',24)
set(0,'DefaultTextFontSize',20)
set(0,'DefaultAxesFontSize',26)
set(0,'DefaultLineLineWidth',1.2);
set(gcf,'color','w');


k = 0;
tol = -1e-1;

for n = 1:n_deflation
    u = solutions(:,n);
    f = sum(u.*basis',1);
    if min(f)>tol % we filter out all the solutions with negative values, as f is to be intended as a density
        norm(F(u,A,B,C))
        k = k+1;
        Z = reshape(f,n_gauss,n_gauss);
        figure(k+1)
        set(gcf,'color','w');
        surf(Xi,Yj,Z);
    end
end
return

%% Stabilizing control

a_infty = solutions(:,2);
rho_infty = sum(a_infty'.*basis,2); % we want to stabilize around an unstable steady state
gamma = 1e-3; %5 

% temporal discretization
tf=1;
dt=0.01;
nt=tf/dt+1;

% two dimensional control u = (u1,u2)
guess = coeffs(xi,yj,wi,zj,l,f_new);
u00 = rand(L,nt,2)*0;

diff=2;
tol=1;
niter=1;

Mass = diag(M);
forward   = @(a,u) odefw_c(a,u,A2,C1);
dforward  = @(a,u) odedfw_c(a,u,A2,C1);
backward  = @(a,p,u) odebw_c(a,p,u,A2,C1,Mass,gamma);
dbackward = @(a,p,u) odedbw_c(a,p,u,A2,C1,Mass);


%% Receding horizon control sinthesys

X_MPC = zeros(L,nt);
u_MPC = zeros(L,nt-1,2);
X_MPC(:,1) = x0;



for t = 1:nt 

    [X00]=IRK(x0,u00,dt,forward,dforward);   %forward pass

    % first iteration
    [P0,dJ0]=IRKadj(X00,u00,dt,backward,dbackward,gamma,Mass);  %backward pass
    u0 = u00*0;
    for k=1:nt-1
        u0(:,k,:) = u00(:,k,:) - .001*dJ0(:,k,:);
    end
    [X0]=IRK(x0,u0,dt,forward,dforward);                       %forward pass  
    disp('first iteration')
    
    
    % second iteration
    [P1,dJ1]=IRKadj(X0,u0,dt,backward,dbackward,gamma,Mass);  %backward pass
    u1 = u0*0;
    for k=1:nt-1
        u1(:,k,:) = u0(:,k,:) - 0.001*dJ1(:,k,:);
    end
    [X]=IRK(x0,u1,dt,forward,dforward);  
    disp('second iteration')
    diff = 2;
    while(diff>tol)
        [P,dJ_new]=IRKadj(X,u1,dt,backward,dbackward,gamma,Mass);  %backward pass
        p = squeeze(P(:,1:nt-1));            %squeeze adjoint
        u_new = u1*0;
        for k=1:nt-1
            alpha_k = (u1(:,k)-u0(:,k))'*(dJ1(:,k)-dJ0(:,k))/norm(dJ1(:,k)-dJ0(:,k))^2; % BB step
            u_new(:,k) = u1(:,k) - alpha_k*dJ_new(:,k); %gradient update
        end

        [X]=IRK(x0,u_new,dt,forward,dforward);       %forward pass

        niter = niter+1
        diff  = norm(u_new(:)-u1(:))/(norm(u1(:)))   %iteration relative error


        %update variables...
        u0 = u1;
        u1 = u_new;
        dJ0 = dJ1;
        dJ1 = dJ_new;
    end
    toc
    x0 = X(:,2);
    X_MPC(:,t) = x0;
    u_MPC(:,t-1,:) = u1(:,1,:);
end

%% uncontrolled simulation

[X_uncontrolled]=IRK(X_MPC(:,1),u00*0,dt,forward,dforward); 

%% comparing the controlled and uncontrolled simulation

figure(2)
norma = zeros(1,nt);
norma_uncontrolled = norma;

k = 1;
for t = 1:nt
    norma(k) = (X_MPC(:,t)-a_infty)'*Mass*(X_MPC(:,t)-a_infty);
    k = k+1;
end

k=1;
for t = 1:nt
    norma_uncontrolled(k) = (X(:,t)-a_infty)'*Mass*(X(:,t)-a_infty);
    k = k+1;
end

plot(1:nt,norma,'ro')
hold on
plot(1:nt,norma_uncontrolled(1:51),'b*')
    


%%
figure(3)

norma = zeros(1,nt);
for n = 1:nt
    a = X_uncontrolled(:,n);
    f = sum(a'.*basis,2);
    norma(n) = norm(f);
    Z = reshape(f,n_gauss,n_gauss);
    surf(Xi,Yj,Z);
    zlim([-0.1,2])
    drawnow
    pause(0.01)
    
    %hold on
end


