%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                    %
%%                          Solving HKB steady states with Newton (and deflation) method                            %% 
%                                                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc
global L nt a_infty D M1 B


alpha = 1;
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

%% Simulation for the initial guess

basis = zeros(n_gauss,L);
for i = 1:L
    f_i = psi{i};
    basis(:,i) = f_i(xi);
end

rng('default')  % For reproducibility

guess0 =rand(size(xi));
guess0 = guess0/(sum(guess0.*wi));
guess = coeffs(xi,wi,l,guess0);

% set(0,'DefaultTextInterpreter','latex')
% set(0,'DefaultLegendInterpreter','latex')
% set(0,'DefaultAxesTickLabelInterpreter','latex')
% set(0,'defaultAxesXGrid','on')
% set(0,'defaultAxesYGrid','on')
% set(0,'DefaultLegendFontSize',24)
% set(0,'DefaultTextFontSize',20)
% set(0,'DefaultAxesFontSize',26)
% set(0,'DefaultLineLineWidth',1.2);
% set(gcf,'color','w');
% t_interval = [0,1];
% [t,u] = ode45(@(t,u) odefcn(t,u,A1,B,C1,M1) , t_interval , guess);
% 
% 
% f = sum(u(end,:).*basis,2);
% figure(1)
% plot(xi,f)
% guess = u(end,:)';

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
% set(0,'defaultAxesXGrid','on')
% set(0,'defaultAxesYGrid','on')
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
Z = @(m1,m2) sum(exp(-beta*V_hkb(xi)+(beta*kappa*(m1*cos(xi)+m2*sin(xi)))).*wi);
rho = @(x,m1,m2) (1/Z(m1,m2))*exp(-beta*V_hkb(x)+(beta*kappa*(m1*cos(x)+m2*sin(x))));
freeE = zeros(size(m1));
Free_energy = @(x,w,m1,m2,rho,Z) sum(V(x).*rho(x,m1,m2).*w) + ...
                                 sum(beta_m1*rho(x,m1,m2).*log(rho(x,m1,m2)).*w)-...
                                 (kappa/2)*convol_complete(x,w,m1,m2,rho,Z);

for j = 1:N
    for i=1:N
        m_2 = m2(i,j);
        m_1 = m1(i,j);
        R = R_complete(m_1,m_2,beta,alpha,kappa);
        sol(i,j,:) = R;
        freeE(i,j) = Free_energy(xi,wi,m_1,m_2,rho,Z);
    end
    
end
%%
figure(2)
set(gcf,'color','w');
pcolor(m1,m2,freeE)
colormap('jet')
caxis([-3  -1.75])
hold on
contourf(m1,m2,sol(:,:,1)-m1,[0,0],'Fill','off','LineColor',[0.8350 0.0780 0.1840],'LineWidth',3)
hold on
contourf(m1,m2,sol(:,:,2)-m2,[0,0],'Fill','off','LineColor',[0.3010 0.7450 0.9330],'LineWidth',3)


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
        c(:,n) = pinv(S'*S)*S'*Y;
        disp([num2str(k),'-th solution: ', '   ',num2str(norm(F(u,A,B,C))), '      ',mat2str(c(:,n)')])       
        k = k+1;
        plot(c(2,n),c(3,n),'x','MarkerSize', 25,'LineWidth',4)
        %plot(x,f','linewidth',3) 
        %xlim([-0,a])
        %ylim([0.13,0.18])
        hold on
    end
end
legend('','$R_1=m_1$','$R_2=m_2$')
xlabel('$m_1$')
ylabel('$m_2$')



hold off
%%
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
        %xlim([-0,a])
        %ylim([0.13,0.2])
        hold on
    end
end
title(['$\beta^{-1} = $' num2str(beta_m1) '$,\;\alpha= $' num2str(alpha) '$,\;\kappa= $' num2str(kappa)])

%% Control matrix D
D = zeros(L,L,L);
for k = 1:L
    f_k = psi{k}(xi);
    for m = 1:L
        df_m = dpsi{m}(xi); 
        for n = 1:L
            f_n  = psi{n}(xi); 
            D(n,k,m) = sum(f_k.*df_m.*f_n.*wi);
        end
    end
end

%% Stabilizing control 
a_infty = solutions(:,1);
rho_infty = sum(a_infty'.*basis,2); % we want to stabilize around an unstable steady state
gamma = 5; %5 
% temporal discretization
tf=5; %3
dt=0.01; %0.015
nt=tf/dt+1;
t=0;

u00=[ones(1,nt-1)*0.1;zeros(L-1,nt-1)];
X=zeros(L,nt);
guess0 =rand(size(xi));
guess0 = guess0/(sum(guess0.*wi));
a0 = coeffs(xi,wi,l,guess0);
x0 = a0-a_infty;

diff=1;
tol=0.01;
niter=1;

forward   = @(a,u) odefw_c(a,u,A1,C1);
dforward  = @(a,u) odedfw_c(a,u,A1,C1);
backward  = @(a,p,u) odebw_c(a,u,p,A1,C1);
dbackward = @(a,p,u) odedbw_c(a,u,p,A1,C1);

% first iteration
[X]=IRK(x0,u00,dt,forward,dforward);                      %forward pass
[P0,dJ0]=IRKadj(X,x0*0,u00,dt,backward,dbackward,gamma);  %backward pass
u0 = u00*0;
for k=1:nt-1
    u0(:,k) = u00(:,k) - 0.01*dJ0(:,k);
end
[X0]=IRK(x0,u0,dt,forward,dforward);                       %forward pass  
disp('first iteration')
% second iteration
[P1,dJ1]=IRKadj(X0,x0*0,u0,dt,backward,dbackward,gamma);  %backward pass
u1 = u0*0;
for k=1:nt-1
    u1(:,k) = u0(:,k) - 0.001*dJ1(:,k);
end
[X]=IRK(x0,u1,dt,forward,dforward);  
disp('second iteration')
% now we have the info for using the BB algorithm
while(diff>tol)
    [P,dJ_new]=IRKadj(X,x0*0,u1,dt,backward,dbackward,gamma);  %backward pass
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

%%
guess0 =rand(size(xi));
guess0 = guess0/(sum(guess0.*wi));
a0 = coeffs(xi,wi,l,guess0);
x0 = a0;


[X_uncontrolled]=IRK(x0,u00,dt,forward,dforward);
norma = zeros(1,101);
norma_uncontrolled = norma;
k = 1;
for t = 1:5:nt
    norma(k) = norm(X(:,k))^2;
    norma_uncontrolled(k) = norm(X_uncontrolled(:,k))^2;
    k = k+1;
end
figure(1)
plot(1:101,norma_uncontrolled,'b*',1:101,norma,'ro')

%%
figure(2)
for n = 1:5:nt-1
    a = X(:,n);
    f = sum(a'.*basis,2);
    plot3(xi,ones(size(xi))*dt*n,f)
    hold on
end
hold off
drawnow
figure(3)
for n = 1:5:nt-1
    a = u0(:,n);
    u = sum(a'.*basis,2);
    plot3(xi,ones(size(xi))*dt*n,u)
    hold on
end
hold off
drawnow


%% Plotting the trajectory
guess0 =rand(size(xi));
guess0 = guess0/(sum(guess0.*wi));
a0 = coeffs(xi,wi,l,guess0);
x0 = a0-a_infty;

forward_sdre   = @(a) odefw_sdre(a,A1,C1,Q,R);
dforward_sdre  = @(a) odedfw_sdre(a,A1,C1,Q,R);
[X,u]=IRK_sdre(x0,dt,forward_sdre,dforward_sdre);

for n = 1:nt
    a = X(:,n);
    f = sum(a'.*basis,2);
    plot3(xi,ones(size(xi))*dt*n,f)
    hold on
end
hold off
figure(2)
for n = 1:nt-1
    a = u(:,n);
    un = sum(a'.*basis,2);
    plot3(xi,ones(size(xi))*dt*n,un)
    hold on
end
hold off


%% SDRE approach 
a_infty = solutions(:,3);
rho_infty = sum(a_infty'.*basis,2); % we want to stabilize around an unstable steady state
gamma = .1;
% temporal discretization
tf=2;
dt=0.01;
nt=tf/dt+1;
Q = eye(L)*1/2;
R = eye(L)*gamma/2;

X=zeros(L,nt+1);
a0 = rand(size(a_infty));
x0 = a0-a_infty;
X(:,1)=x0;
u = zeros(L,nt+1);


%%

for n = 1:nt
    a = X0(:,n);
    f = sum(a'.*basis,2);
    plot3(xi,ones(size(xi))*dt*n,f)
    hold on
end
hold off
%%
figure(2)
for n = 1:8%nt-1
    a = u(:,n)
    ut = sum(a'.*basis,2);
    plot3(xi,ones(size(xi))*dt*n,ut)
    hold on
end
hold off

%%
X=zeros(L,nt+1);
guess0 =rand(size(xi));
guess0 = guess0/(sum(guess0.*wi));
a0 = coeffs(xi,wi,l,guess0);
x0 = a0-a_infty;

X(:,1)=x0;
for t=1:nt
    [b_a,g_a,~,] = semilin(x0);
    a = M1*(-A1-C1-b_a);
    b = M1*(-g_a);
    [K,~,~] = lqr(a,b,Q,R);
    u0     = -K*x0;
    f = a*x0+b*u0;
    x1 = x0+dt*f;
    X(:,t+1) = x1;
    x0=x1;
end

