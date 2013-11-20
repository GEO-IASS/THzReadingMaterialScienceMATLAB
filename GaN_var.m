[t_GaN,E_GaN] = open_picotd('GaN.picotd',4096);
t_GaN(4092,:) = [ ];

%trunc and padding
t_GaN_trunc = 0:0.078125:(4095*0.078125);
E_GaN_trunc(1:4096,1)=[0];
E_GaN_trunc(1857:2172,1) = E_GaN(1857:2172,1);

%ffts and phase
F_GaN = fft(E_GaN_trunc); 
F_GaN = F_GaN ./F_reference;
%F_GaN = fft(E);
F_GaN_abs = abs(F_GaN);
F_GaN_phase = unwrap(angle(F_GaN));
w_GaN = 0:(1/(4095*0.078125)):(1/(0.078125));
w_GaN = w_GaN'.* 10^12;

%Solving for epsilon
d_GaN = 0.527 * 10^-3;
c = 3 * 10^8;
x0_GaN= 12* ones(4096,1);
p_GaN = zeros(4096,1);
%x = zeros(4096,1);
for i=1:4096
    options = optimset('Display','off');  % Turn off display
    f=@(x) (F_GaN(i,1))-4*x*(exp(-1i*w_GaN(i,1)*d_GaN*(x-1)/c))*(1/(1+x)^2);
    p_GaN(i) = fsolve(f, x0_GaN(i),options);
end

plot(w_GaN,real(p_GaN))
