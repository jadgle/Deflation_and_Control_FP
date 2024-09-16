function [ga1,ga2] = g(a)
    global L D1 D2 M1
    Da1 = zeros(L,L);
    Da2 = zeros(L,L);
    for m = 1:L
        Da1(m,:) = a'*(D1(:,:,m));
        Da2(m,:) = a'*(D2(:,:,m));
    end
    ga1 = M1*(-Da1);
    ga2 = M1*(-Da2);
end