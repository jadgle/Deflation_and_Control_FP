function dfun = dF(u,A,B,C)
    L = length(u);
    dB = zeros(L,L);
    for n = 1:L
        for m = 1:L
            dB(n,m) = u'*B(:,m,n)+B(m,:,n)*u;
        end
    end 
    dB = [dB;zeros(1,L)];
    dfun = -A-C-dB;
end