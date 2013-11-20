
E_transfer = E_silicon./E_reference;
%trunc and padding

%ffts and phase
F_transfer = fft(E_transfer,3);
F_transfer_abs = abs(F_transfer);
F_transfer_phase = unwrap(angle(F_transfer));

plot(t_reference,E_transfer);
plot(w_reference,F_transfer);
plot(w_reference,F_transfer_phase);