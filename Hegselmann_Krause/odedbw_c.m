function dadt = odedbw_c(a,p,u,A,C,M)
    global L D M1 B a_infty
    dB = zeros(L,L);
    dD = zeros(L,L);
    for n = 1:L
        dD(:,n) = u'*D(:,:,n);
        dB(:,n) = a'*(B(:,:,n)+B(:,:,n)');
    end 
    dadt = (M1*(-A-C-dB-dD))';
end