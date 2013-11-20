[t,E_ref] = open_picotd('Reference.picotd',4096);
t(4092,:) = [ ];
t_reference = t;
E_reference = E_ref;
%trunc and padding
t_silicon_trunc = 0:0.078125:(4095*0.078125);
E_reference_trunc(1:4096,1)=[0];
E_reference_trunc(1857:2172,1) = E_reference(1857:2172,1);

%ffts and phase
F_reference = fft(E_reference_trunc);
%F_reference = fft(E_reference);
F_reference_abs = abs(F_reference);
F_reference_phase = unwrap(angle(F_reference));

w_reference = 0:(1/(4095*0.078125)):(1/(0.078125));
w_reference = w_reference'.* 10^12;


[t_silicon,E_silicon] = open_picotd('Silicon.picotd',4096);
t_silicon(4092,:) = [ ];

%trunc and padding
t_silicon_trunc = 0:0.078125:(4095*0.078125);
E_silicon_trunc(1:4096,1)=[0];
E_silicon_trunc(1857:2172,1) = E_silicon(1857:2172,1);

%ffts and phase
F_silicon = fft(E_silicon_trunc); 
F_silicon = F_silicon ./F_reference;
%F_silicon = fft(E);
F_silicon_abs = abs(F_silicon);
F_silicon_phase = unwrap(angle(F_silicon));
w_silicon = 0:(1/(4095*0.078125)):(1/(0.078125));
w_silicon = w_silicon'.* 10^12;

%Solving for epsilon
d_silicon = 0.527 * 10^-3;
c = 3 * 10^8;
x0_silicon= 3.42* ones(4096,1);
p_silicon = zeros(4096,1);
%x = zeros(4096,1);
for i=1:4096
    options = optimset('Display','off');  % Turn off display
    f=@(x) (F_silicon(i,1))-4*x*(exp(-1i*w_silicon(i,1)*d_silicon*(x-1)/c))*(1/(1+x)^2);
    p_silicon(i) = fsolve(f, x0_silicon(i),options);
end
subplot(2,2,1);
plot(w_silicon,real(p_silicon));
title('Silicon_n_real vs w');

subplot(2,2,2);
plot(w_silicon(1:76,1),F_silicon_phase(1:76,1));
title('Silicon phaseVsW');

subplot(2,2,3);
plot(w_silicon(1:76,1),F_silicon_abs(1:76,1));
title('Silicon absVsW');

