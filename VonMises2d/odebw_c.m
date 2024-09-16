function dadt = odebw_c(a,p,u,A,C,M,gamma)
    global L D1 D2 M1 B a_infty
    dB = zeros(L,L);
    dD1 = zeros(L,1);
    dD2 = zeros(L,1);
    u1 = u(:,:,1);
    u2 = u(:,:,2);
    for n = 1:L
        dD1(:,n) = u1'*D1(:,:,n);
        dD2(:,n) = u2'*D2(:,:,n);
        dB(:,n) = a'*(B(:,:,n)+B(:,:,n)');
    end 
    dadt = ((M1*(-A-C-dB-dD1-dD2))'*p+2*gamma*M*(a-a_infty));
end

