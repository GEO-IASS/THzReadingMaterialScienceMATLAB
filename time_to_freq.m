function [f,E_f]=time_to_freq(t, E_t, N,m)
Tm=(t(3)-t(2))*1e-12;
fm=1/Tm;
E_f=fft(E_t,N*m);
f_pod=round(length(E_f)/2);
E_f=E_f(1:f_pod);
f=(0:f_pod-1);
f=(fm/N).*f;
f=f/(m*1e12);
end