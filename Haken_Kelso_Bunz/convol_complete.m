function c = convol_complete(x,w,m1,m2,rho,Z)
c = zeros(size(x));
for i = 1:length(x)
    x_bar = x(i);
    c(i) = sum(cos(x-x_bar).*rho(x,m1,m2)*rho(x_bar,m1,m2).*w); 
end
c = sum(c.*w);