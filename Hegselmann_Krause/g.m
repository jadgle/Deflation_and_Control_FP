function ga = g(a)
    global L D M1
    Da = zeros(L,L);
    for m = 1:L
        Da(m,:) = a'*D(:,:,m)';
    end
    ga = M1*(-Da);
end