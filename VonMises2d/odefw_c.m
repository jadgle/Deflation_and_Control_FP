function dadt = odefw_c(a,u,A,C)
global L D1 D2 M1 B
    Baa = zeros(L,1);
    uD1a = zeros(L,1);
    uD2a = zeros(L,1);
    u1 = u(:,:,1);
    u2 = u(:,:,2);
    for m = 1:L
        Baa(m) = a'*B(:,:,m)*a;
        uD1a(m) = u1'*D1(:,:,m)*a;
        uD2a(m) = u2'*D2(:,:,m)*a;
    end
    dadt = M1*(-(A+C)*a-Baa-uD1a-uD2a);
end