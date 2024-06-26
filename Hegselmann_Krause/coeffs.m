function u = coeffs(x,w,L,rho0)

a = zeros(L,1);
b = zeros(L,1);
for i = 1:L
    cosi = cos(i*2*pi*x);
    sini = sin(i*2*pi*x);
    a(i) = 2*(w.*sini)'*rho0;
    b(i) = 2*(w.*cosi)'*rho0;
end
a0 = 2*sum(w.*rho0);
u = [a0/2;a;b];