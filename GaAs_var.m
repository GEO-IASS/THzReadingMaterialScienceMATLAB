[t_GaAs,E_GaAs] = open_picotd('GaAs.picotd',4096);
t_GaAs(4092,:) = [ ];

%trunc and padding
t_GaAs_trunc = 0:0.078125:(4095*0.078125);
E_GaAs_trunc(1:4096,1)=[0];
E_GaAs_trunc(1857:2172,1) = E_GaAs(1857:2172,1);

%ffts and phase
F_GaAs = fft(E_GaAs_trunc); 
F_GaAs = F_GaAs ./F_reference;
%F_GaAs = fft(E);
F_GaAs_abs = abs(F_GaAs);
F_GaAs_phase = unwrap(angle(F_GaAs));
w_GaAs = 0:(1/(4095*0.078125)):(1/(0.078125));
w_GaAs = w_GaAs'.* 10^12;

%Solving for epsilon
d_GaAs = 0.527 * 10^-3;
c = 3 * 10^8;
x0_GaAs= 12* ones(4096,1);
p_GaAs = zeros(4096,1);
%x = zeros(4096,1);
for i=1:4096
    options = optimset('Display','off');  % Turn off display
    f=@(x) (F_GaAs(i,1))-4*x*(exp(-1i*w_GaAs(i,1)*d_GaAs*(x-1)/c))*(1/(1+x)^2);
    p_GaAs(i) = fsolve(f, x0_GaAs(i),options);
end

plot(w_GaAs,real(p_GaAs))
