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