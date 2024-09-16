function [x]=IRK(x0,u,h,f,df)
global L nt

max_iterations = 2e3;
tol = 1e-3;
x   = zeros(L,nt);
x(:,1)=x0;

a11 = 1/4;
a12 = 1/4 - sqrt(3)/6;
a21 = 1/4 + sqrt(3)/6;
a22 = 1/4;
GL_iterate = @(k1,k2,w_n,u_n) [k1 - f(w_n + h*(k1*a11 + k2*a12),u_n);
                               k2 - f(w_n + h*(k1*a21 + k2*a22),u_n)];
                       
err = @ (k1,k2,w_n,u_n) norm(GL_iterate(k1,k2,w_n,u_n));

for k=1:nt-1
    k
    w_n = x(:,k);
    u_n = u(:,k,:);

    iteration = 0;
    x_ = f(w_n,u_n);
    x1_guess = w_n + (1/2 - sqrt(3)/6)*h*x_;
    x2_guess = w_n + (1/2 + sqrt(3)/6)*h*x_;
    k1 = f(x1_guess,u_n);
    k2 = f(x2_guess,u_n);
    
    while (err(k1,k2,w_n,u_n) > tol && iteration < max_iterations)
        iteration = iteration + 1;
        df_n1 = df(w_n+a11*h*k1+a12*h*k2,u_n);
        df_n2 = df(w_n+a21*h*k1+a22*h*k2,u_n);
        d_GL = [eye(L)-df_n1*a11*h, -df_n1*a12*h; -df_n2*a21*h, eye(L)-df_n2*a22*h];
        GL = GL_iterate(k1,k2,w_n,u_n);
        k_next = [k1;k2] - pinv(d_GL)*GL;
        k1 = k_next(1:L);
        k2 = k_next(L+(1:L));
    end   
    if err(k1,k2,w_n,u_n) > tol
        error('Newton did not converge at timestep %d.', k);
        
    end
    new_x = w_n + h / 2 * (k1 + k2);
    if any(isnan(new_x))
        error('NaN at timestep %d',k)
    else
        x(:,k+1) = new_x;
    end
    
end