%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                    %
%%                          Solving HKB steady states with Newton (and deflation) method                            %% 
%                                                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
global L nt a_infty D M1 B


alpha = -1;
beta_m1 = 1;
beta= 1/beta_m1;

kappa = 5;

%%%%%%
% 
W = @(x) -cos(x); 
h = @(x) 2*alpha*sin(2*x);
V = @(x) alpha*cos(2*x);

%% Discretization

l = 50; %number of modes

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

%% Basis functions and initial guess

basis = zeros(n_gauss,L);
for i = 1:L
    f_i = psi{i};
    basis(:,i) = f_i(xi);
end

rng('default')  % For reproducibility

guess0 =rand(size(xi));
guess0 = guess0/(sum(guess0.*wi));
guess = coeffs(xi,wi,l,guess0);

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
p = 4;       % deflation power
shift = 1; % deflation shift
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

N = 100;
m1 = linspace(-2,2,N);
m2 = m1;
[m1,m2] = meshgrid(m1,m2);

sol = zeros(N,N,2);
V_hkb = @(x) alpha*cos(2*x);

% fixed points of the Euler-Lagrange equations for the free energy functional 
Z = @(m1,m2) sum(exp(-beta*V_hkb(xi)+(beta*kappa*(m1*cos(xi)+m2*sin(xi)))).*wi);
rho = @(x,m1,m2) (1/Z(m1,m2))*exp(-beta*V_hkb(x)+(beta*kappa*(m1*cos(x)+m2*sin(x))));

% free energy functional, used for determining the stability of the steady states
freeE = zeros(size(m1));
Free_energy = @(x,w,m1,m2,rho,Z) sum(V(x).*rho(x,m1,m2).*w) + ...
                                 sum(beta_m1*rho(x,m1,m2).*log(rho(x,m1,m2)).*w)-...
                                 (kappa/2)*convol_complete(x,w,m1,m2,rho,Z);

% determining the self consistency system (m1 = R1(m1,m2) and m2 = R2(m1,m2))
for j = 1:N
    for i=1:N
        m_2 = m2(i,j);
        m_1 = m1(i,j);
        R = R_complete(m_1,m_2,beta,alpha,kappa);
        sol(i,j,:) = R;
        freeE(i,j) = Free_energy(xi,wi,m_1,m_2,rho,Z);
    end
    
end

figure(2)
set(gcf,'color','w');
p1 = pcolor(m1,m2,freeE);
p1.EdgeAlpha = 0.1;
colormap('summer')
caxis([-3  -1.75])
hold on
contourf(m1,m2,sol(:,:,1)-m1,[0,0],'Fill','off','LineColor',[0.8350 0.0780 0.1840],'LineWidth',3)
hold on
contourf(m1,m2,sol(:,:,2)-m2,[0,0],'Fill','off','LineColor',[0.3010 0.7450 0.9330],'LineWidth',3)

%% LSS estimation of the coefficients (m1,m2) associated with the identified steady states

S = [ones(n_gauss,1),beta*kappa*cos(xi),beta*kappa*sin(xi)];
c = [];

k=0;
disp('----- Norm of the Residual and |e^b-1/Z|  -----')
tol = 3e-2;
for n = 1:n_deflation
    u = solutions(:,n);
    f = sum(u'.*basis,2);
    
    if sum(f<0)<100 && norm(F(u,A,B,C))<1e-3
        Y = log(f)+beta*V(xi);
        c(:,n) = pinv(S'*S)*S'*Y; %LSS
        disp([num2str(k),'-th solution: ', '   ',num2str(norm(F(u,A,B,C))), '      ',mat2str(c(:,n)')])       
        k = k+1;
        plot(c(2,n),c(3,n),'x','MarkerSize', 25,'LineWidth',4)
        hold on
    end
end
legend('','$R_1=m_1$','$R_2=m_2$')
xlabel('$m_1$')
ylabel('$m_2$')


hold off

%% Plotting the identified steady states
figure(3)
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
    u = solutions(:,n);
    f = sum(u'.*basis,2);
    
    if sum(f<0)<100 && norm(F(u,A,B,C))<1e-3
        plot(xi,f','linewidth',3) 
        hold on
    end
end
title(['$\beta^{-1} = $' num2str(beta_m1) '$,\;\alpha= $' num2str(alpha) '$,\;\kappa= $' num2str(kappa)])

