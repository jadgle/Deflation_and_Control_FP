function dF = odedfw_c(a,u,A,C)
    global L D M1 B
    dB = zeros(L,L);
    dD = zeros(L,L);
    for n = 1:L
        dD(:,n) = D(:,:,n)'*u;
        dB(:,n) = (B(:,:,n)+B(:,:,n)')*a;
    end 
    dF = M1*(-A-C-dB-dD);
end
