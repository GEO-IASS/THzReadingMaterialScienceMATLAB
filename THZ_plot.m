 
%plot(t_reference, E_reference,t_silicon,E_silicon);
%plot(w_reference(1:2046,1), F_reference_phase(1:2046,1),'r',w_silicon(1:2046,1),F_silicon_phase(1:2046,1))
%plot(w_silicon(1:76,1),F_silicon_phase(1:76,1));
%plot(w_silicon(1:76,1),F_silicon_abs(1:76,1));

%plot(w_reference(1:2046,1), F_reference_abs(1:2046,1),w_silicon(1:2046,1),F_silicon_abs(1:2046,1))
%semilogy(w_reference(1:2045,1), F_reference_abs(1:2045,1),w_silicon(1:2045,1),F_silicon_abs(1:2045,1))
%semilogy(w_reference(1:2045,1), F_reference_phase(1:2045,1),w_silicon(1:2045,1),F_silicon_phase(1:2045,1))
%grid on

 plot(w_reference, F_reference_abs,w_silicon, F_silicon_abs)
%plot(w_silicon(1:2046,1),F_silicon_phase(1:2046,1))