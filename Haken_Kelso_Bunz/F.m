function fun = F(u,A,B,C)
    L = length(u);
    Buu = zeros(L,1);
    for m = 1:L
        Buu(m) = u'*B(:,:,m)*u;
    end
    Buu = [Buu;-1];
    fun = -A*u-Buu-C*u;
    
end