function F = transfer_fun(n)
y(1:4096,1) = [0];
x = 4096;
d = 0.527 * 10^-3;
c = 3 * 10^8;
for i=1:x
    
  F(i,1) = [(E_silicon(i,1))-4*n*(e^(-1i*w(i,1)*d*(n-1)/c)*(1/(1+n)^2))]
end
