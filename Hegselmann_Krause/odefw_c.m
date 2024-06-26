function dadt = odefw_c(a,u,A,C)
global L D M1 B
    Baa = zeros(L,1);
    uDa = zeros(L,1);
    for m = 1:L
        Baa(m) = a'*B(:,:,m)*a;
        uDa(m) = u'*D(:,:,m)*a;
    end
    dadt = M1*(-(A+C)*a-Baa-uDa);
end