[t_ZnS,E_ZnS] = open_picotd('ZnS.picotd',4096);
t_ZnS(4092,:) = [ ];

%trunc and padding
t_ZnS_trunc = 0:0.078125:(4095*0.078125);
E_ZnS_trunc(1:4096,1)=[0];
E_ZnS_trunc(1857:2172,1) = E_ZnS(1857:2172,1);

%ffts and phase
F_ZnS = fft(E_ZnS_trunc); 
F_ZnS = F_ZnS ./F_reference;
%F_ZnS = fft(E);
F_ZnS_abs = abs(F_ZnS);
F_ZnS_phase = unwrap(angle(F_ZnS));
w_ZnS = 0:(1/(4095*0.078125)):(1/(0.078125));
w_ZnS = w_ZnS'.* 10^12;

%Solving for epsilon
d_ZnS = 0.527 * 10^-3;
c = 3 * 10^8;
x0_ZnS= 12* ones(4096,1);
p_ZnS = zeros(4096,1);
%x = zeros(4096,1);
for i=1:4096
    options = optimset('Display','off');  % Turn off display
    f=@(x) (F_ZnS(i,1))-4*x*(exp(-1i*w_ZnS(i,1)*d_ZnS*(x-1)/c))*(1/(1+x)^2);
    p_ZnS(i) = fsolve(f, x0_ZnS(i),options);
end

plot(w_ZnS,real(p_ZnS))
