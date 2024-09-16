function dF = odedfw_c(a,u,A,C)
    global L D1 D2 M1 B
    dB = zeros(L,L);
    dD1 = zeros(L,L);
    dD2 = zeros(L,L);
    u1 = u(:,:,1);
    u2 = u(:,:,2);
    for n = 1:L
        dD1(:,n) = u1'*D1(:,:,n);
        dD2(:,n) = u2'*D2(:,:,n);
        dB(:,n) = a'*(B(:,:,n)+B(:,:,n)');
    end 
    dF = M1*(-A-C-dB-dD1-dD2);
end
