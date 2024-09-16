function [p,dJ]=IRKadj(a,u,h,f,df,gamma,M)
global L nt a_infty
lambda = 1000;
max_iterations = 2e3;
tol = 1e-3;

a11 = 1/4;
a12 = 1/4 - sqrt(3)/6;
a21 = 1/4 + sqrt(3)/6;
a22 = 1/4;

GL_iterate = @(k1,k2,w_n,p_n,u_n) [k1 - f(w_n,p_n + h*(k1*a11 + k2*a12),u_n);
                                   k2 - f(w_n,p_n + h*(k1*a21 + k2*a22),u_n)];
                       
er = @ (k1,k2,w_n,p_n,u_n) norm(GL_iterate(k1,k2,w_n,p_n,u_n))^2; 

adj=zeros(L,nt);
aux=zeros(L,nt);
u0 = zeros(L,nt-1,2);
p0 = 2*lambda*M*(a(:,end)-a_infty);
adj(:,1)=p0;
for k=1:nt
    aux(:,k)=a(:,nt-k+1);
    u0(:,k,:)=u(:,nt-k+1,:);
end 

for k=1:nt-1
    
    p_n = adj(:,k);
    w_n = aux(:,k);
    u_n = u0(:,k,:);
    
    %f_n = f(w_n,p_n,u_n);
    f_n = f(w_n,p_n,u_n);
    p1_guess = p_n + (1/4 - sqrt(3)/6)*h*f_n;
    p2_guess = p_n + (1/4 + sqrt(3)/6)*h*f_n;
    k1 = f(w_n,p1_guess,u_n);
    k2 = f(w_n,p2_guess,u_n);
    iteration = 0;
    while (er(k1,k2,w_n,p_n,u_n) > tol && iteration < max_iterations)
        iteration = iteration + 1;
        df_n1 = df(w_n,p_n+a11*h*k1+a12*h*k2,u_n);
        df_n2 = df(w_n,p_n+a21*h*k1+a22*h*k2,u_n);
        d_GL = [eye(L)-df_n1*a11*h, -df_n1*a12*h; 
                -df_n2*a21*h, eye(L)-df_n2*a22*h];
        GL = GL_iterate(k1,k2,w_n,p_n,u_n);
        k_next = [k1;k2] - linsolve(d_GL, GL);
        k1 = k_next(1:L);
        k2 = k_next(L+(1:L));
    end   
    if er(k1,k2,w_n,p_n,u_n) > tol
        disp(['Newton did not converge at timestep ',num2str(k)])
        return
        
    end
    new_adj = p_n + h / 2 * (k1 + k2);
    if any(isnan(new_adj))
        error('NaN at timestep %d',k)
    else
        adj(:,k+1) = new_adj;
    end
    
end

p  = aux;
dJ = u;
for k=1:nt
    p_n=adj(:,nt-k+1);
    w_n = a(:,k);
    [ga1,ga2] = g(w_n);
    
    dJ(:,k,1) = (2*gamma*u(:,k,1)'*M + p_n'*ga1);
    dJ(:,k,2) = (2*gamma*u(:,k,2)'*M + p_n'*ga2);
    p(:,k) = p_n;
end

