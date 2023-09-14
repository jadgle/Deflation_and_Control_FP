function dadt = odefw_c(a,u,A,C)
global L a_infty D M1 B
    Baa = zeros(L,1);
    uDa = zeros(L,1);
    uDa_infty = zeros(L,1);
    for m = 1:L
        Baa(m) = a'*B(:,:,m)*a;
        uDa(m) = u'*D(:,:,m)*a;
        uDa_infty(m) = u'*D(:,:,m)*a_infty;
    end
    dadt = M1*(-(A+C)*a-Baa-uDa-uDa_infty);
end