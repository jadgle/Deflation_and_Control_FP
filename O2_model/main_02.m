%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                    %
%%                       Solving O2 model steady states with Newton (and deflation) method                          %% 
%                                                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc




beta_m1 = 0.415; %0.25, 0.4

alpha = 0.05;
kappa = 1;

%%%%%%
% comparing with reference paper, to see if we can parametrize the solutions using bessel functions
% formulation
W = @(x)-cos(x); 
h = @(x) -alpha*sin(x);
V = @(x) -alpha*cos(x);

%% Discretization

l = 10; %number of modes

% basis functions
psi  = cell(2*l+1,1);
dpsi = cell(2*l+1,1);

psi{1}  = @(y) ones(size(y));
dpsi{1} = @(y) zeros(size(y));

for i = 1:l
        psi{i+1}   = @(y) sin(i*y);
        psi{l+i+1} = @(y) cos(i*y);
        dpsi{i+1}   = @(y) i*cos(i*y);
        dpsi{l+i+1} = @(y) -i*sin(i*y);
end



%% Generating the system matrices

% prepping for Gauss quadratures
n_gauss = 100;
[xi,wi]=Gauss_quad(n_gauss,0,2*pi);

% linear components of the equation 
L = size(psi,1); 
C1 = zeros(L);
A1 = zeros(L,1);
M = zeros(L,1);

for i = 1:L
    df_i = dpsi{i};
    f_i  = psi{i};
    A1(i)       =  sum(df_i(xi).^2.*wi);
    M(i)        =  sum(f_i(xi).^2.*wi);
    for j = 1:L
        df_j = dpsi{j};
        C1(j,i)        =  sum(-h(xi).*f_i(xi).*df_j(xi).*wi);
    end
end
A1 = beta_m1*diag(A1);
M1 = diag(M.^-1);


% bilinear term
B = zeros(L,L,L); 

convo_integrand = @(s,x_bar,df_k) (df_k(x_bar-s).*(W(s)));
integrand1 = @(s,f_n,df_m,conv) f_n(s).*conv.*df_m(s);

tic
for k = 1:L
    convolution = zeros(n_gauss,1);
    df_k = dpsi{k};
    for i = 1:n_gauss
        x_bar = xi(i);
        convolution(i) = sum(convo_integrand(xi,x_bar,df_k).*wi);
    end

    for m = 1:L
        df_m = dpsi{m}; 
        for n = 1:L
            f_n  = psi{n}; 
            B(n,k,m) = sum(integrand1(xi,f_n,df_m,convolution).*wi);
        end
    end
end
toc
B = kappa*B;

%% Simulation for the initial guess

basis = zeros(n_gauss,L);
for i = 1:L
    f_i = psi{i};
    basis(:,i) = f_i(xi);
end

rng('default')  % For reproducibility

guess = zeros(L,1);
guess(1) = 0.152;

%% Newton with deflation

A = [A1;zeros(1,L)]; % integral=1 condition
C = [C1;zeros(1,L)]; % integral=1 condition

for i = 1:L
    f_i = psi{i};
    A(end,i) = sum(f_i(xi).*wi);
end

update = 100;         % turning off the convergence flag
tolerance = 1e-15;    % tolerance for convergence
max_iteration = 1e3;  % number of before assuming failed convergence

n_deflation = 0; % counter for the number of solutions
iteration = 0;   % counter for the number of iterations

solutions = [];  % list of solutions

% parameters and start for deflation
p = 2;       % deflation power
shift = .05; % deflation shift
eta0 = @(u) 1; 
new_eta = @(u) 1;
deta0 = @(u) zeros(size(u));
new_deta = @(u) zeros(size(u));

trying = 1;
while trying < 2 % this is trying>2 is useful when we are changing initial condition
    u = guess;
    while (update>tolerance && iteration < max_iteration)% && norm(F(u,A,B,C))>1e-10) 
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
        guess =rand(size(u));   % exploring different initial guess
        %guess = u_uniform; % checking if the uniform distribution is a solution    

    else   % if the previous loop has finished with convergence, we deflate w.r.t. the reached root
        % updating the deflation eta
        deta0 =  @(r) new_deta(r)*eta0(r)+new_eta(r)*deta0(r);
        eta0  =  @(r) new_eta(r)*eta0(r);
        new_eta  =  @(r) norm(r-u_new).^p;
        new_deta =  @(r) (p*norm(r-u_new)^(p-2))*(r-u_new); % problem with dimension in dG, see handwritten notes 
        
        disp(['Solution for the ',int2str(n_deflation),'-th deflation found after ',int2str(iteration),' iterations'])
        n_deflation = n_deflation + 1;
        solutions = [solutions, u];
        
        update    = 100;  % reset of the convergence flag 
        iteration = 0;    % reset of the divergence  flag
    end
     
end
%% Plotting and verification of solutions
figure(2)
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
for n = 1:n_deflation
    u = solutions(:,n_deflation-n+1);
    f = sum(u'.*basis,2);
    
    if min(f)>-0.1 && norm(F(u,A,B,C))<1e-3
        plot(xi,f','linewidth',3)
        hold on
    end
end
title(['$\beta^{-1} = $' num2str(beta_m1) '$,\;\eta= $' num2str(alpha) '$,\;\kappa= $' num2str(kappa)])

%%


S = [ones(n_gauss,1)*-1,cos(2*pi*xi)];
c = [];
labels = [];
k=0;
disp('------------- |alpha|  -------------')
tol = 3e-2;
for n = 1:n_deflation
    u = solutions(:,n);
    f = sum(u'.*basis,2);
    
    if min(f)>-0.1 && norm(F(u,A,B,C))<1e-3
        Y = log(f);
        c(:,n) = pinv(S'*S)*S'*Y;
        disp([num2str(k),'-th solution: ', '   ',mat2str(real(c(2,n)))])       
        k = k+1;
    end
end
