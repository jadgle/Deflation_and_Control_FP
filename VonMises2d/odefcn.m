function dudt = odefcn(t,u,A,B,C,M1)
    L = size(A,2);
    Buu = zeros(L,1);
    for m = 1:L
        Buu(m) = u'*B(:,:,m)*u;
    end
    dudt = M1*(-(A+C)*u-Buu);
end