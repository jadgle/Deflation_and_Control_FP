function dadt = odebw_c(a,p,u,A,C)
    global L D M1 B
    dB = zeros(L,L);
    dD = zeros(L,L);
    for n = 1:L
        dD(:,n) = D(:,:,n)'*u;
        dB(:,n) = (B(:,:,n)+B(:,:,n)')*a;
    end 
    dadt = M1*(-A-C-dB-dD)*p-a;
end