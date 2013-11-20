 clear all;
 clc;
%width of sample
d=527e-6;
cc=3e8;%speed of light
%reading the results
[t,E]=open_picotd('Silicon.picotd',4101);
[t_ref,E_ref]=open_picotd('Reference.picotd',4101);
%
%plotting the original results
figure(1)
plot(t_ref(1:4096),E_ref,'b','Linewidth',1);
hold on;
plot(t(1:4096),E,'r','Linewidth',1);
xlabel('time(psecond)','FontWeight','bold');
ylabel('strenght','FontWeight','bold');
title('original results','FontWeight','bold');
legend('Reference','Silicon');
grid on;
L=2;%fft parameter
width=200;%the cutting parameter
%cutting the results
[x,peak]=min(E);
[x,peak_ref]=min(E_ref);
dt=t(peak)-t(peak_ref);
dt=dt*1e-12;
start=peak-width/2;
stop=peak+width/2;
start_ref=peak_ref-width/2;
stop_ref=peak_ref+width/2;
E_trunc=E(start:stop);
E_trunc_ref=E_ref(start_ref:stop_ref);
t_trunc=t(start:stop);
t_trunc_ref=t_ref(start_ref:stop_ref);
%plott cutted information
figure(2);
plot(t_trunc_ref,E_trunc_ref,'b','Linewidth',2);
hold on;
plot(t_trunc,E_trunc,'r','Linewidth',2);
xlabel('time(psecond)','FontWeight','bold');
ylabel('strenght','FontWeight','bold');
title('original results','FontWeight','bold');
legend('Reference','Silicon');
grid on;
%
%adding zeros
pad(start:stop)=E_trunc;
pad_ref(start_ref:stop_ref)=E_trunc_ref;

% fft
[f_data,E_f]=time_to_freq(t_trunc,pad,4101,L);
[f_data,E_f_ref]=time_to_freq(t_trunc_ref,pad_ref,4101,L);
%division
H=E_f./E_f_ref;
% confining the frequnecey 
freq_max=1;
while db(E_f_ref(freq_max))>-40
      freq_max=freq_max+1;
end
freq_min=1;
while(f_data(freq_min)<0.2)
      freq_min=freq_min+1;
end

% cutting the fft of results
E_f_trunc=E_f(1:freq_max);
E_f_ref_trunc=E_f_ref(1:freq_max);
f_data_trunc=f_data(1:freq_max);
H_trunc=H(1:freq_max);

% unwraping the phase
p_raw=unwrap(phase(E_f_trunc));
p_ref=unwrap(phase(E_f_ref_trunc));
p_norm=unwrap(phase(H_trunc));
range=length(f_data_trunc);
%plotting the 
figure(3)
plot(f_data_trunc(50:range),abs(H_trunc(50:range)),'b','Linewidth',2);
 xlabel('Frequency(THz)','FontWeight','bold')
 ylabel('Abs(H)','FontWeight','bold')
 title('magnitude for H','FontWeight','bold')
grid on
% plotting phase
 figure(4)
plot(f_data_trunc(120:freq_max),p_norm(120:freq_max),'b')
 xlabel('frequncey (THz)','FontWeight','bold')
 ylabel('Phase"degree"','FontWeight','bold')
 title('Phase','FontWeight','bold')
 grid on
 %plotting the sample ad refrence db
 figure(5)
plot(f_data_trunc(50:freq_max), db(E_f_ref_trunc(50:freq_max)),'b','Linewidth',2)
 xlabel('Frequency (THz)','FontWeight','bold')
 ylabel('Filed Strength (db)','FontWeight','bold')
 hold on
plot(f_data_trunc(50:freq_max),db(E_f_trunc(50:freq_max)),'r','Linewidth',2)
  title('FFT of reference and Sample','FontWeight','bold')
  legend('reference', 'sample','FontWeight','bold');
  grid on
%finding n,k
 size=length(f_data);
n=1-cc/(2*pi*d*1e12)*phase(H)./f_data;
alpha=-2/d*log(abs(H).*(n+1).*(n+1)./n/4);
k=alpha*cc/4/1e12/pi./f_data;
eps_real=n.*n-k.*k;
eps_imag=2*n.*k;
%plotting n,k
figure(6)
 plot(f_data(20:1:1000),n(20:1:1000),'b','Linewidth',2)
 xlabel('frequencey(THz)','FontWeight','bold')
 ylabel('index of refraction(real part)','FontWeight','bold')
figure(7)
 plot(f_data(200:1000), k(200:1000), 'r','Linewidth',2)
 xlabel('Frequencey(THz)','FontWeight','bold')
 ylabel('index of reflaction(Imaginary part)','FontWeight','bold')
%plotting eps1 and eps 2
figure(8)
 plot(f_data(200:1:1000),eps_real(200:1:1000),'b','Linewidth',2)
 xlabel('frequencey(THz)','FontWeight','bold')
 ylabel('epsilon 1','FontWeight','bold')
 figure(9)
 plot(f_data(200:1000), eps_imag(200:1000), 'r','Linewidth',2)
 xlabel('Freq (THz)','FontWeight','bold')
 ylabel('epsilon 2','FontWeight','bold')
 %
 ninital=cc*dt/(d-cc*dt);
 %
 con=n.*k.*f_data;
 figure(10)
 plot(f_data(200:1:1000),con(200:1:1000),'b','Linewidth',2)
 xlabel('Freq (THz)','FontWeight','bold')
 ylabel('Conductivity','FontWeight','bold')
 grid on
 
%H=@(nsample) 4.*nsample./(1+nsample).^2.*exp(-1i*2*pi.*f_data./cc.*(d*nsample-d))-H_trunc(1,a123);
%a=fsolve(H,q(1,a123));
%end

  q=n+1i*k;
for a123=200:1000;
H1=@(nsample) (4*nsample/(1+nsample)^2*exp(-1i*2*pi*f_data(1,a123)*1e12/cc*(d*nsample-d)))-H(1,a123);
a(1,a123)=fsolve(H1,q(1,a123));
end
%
for b123=200:1000;
    H1=@(nsample) (4*nsample/(1+nsample)^2*exp(-1i*2*pi*f_data(1,b123)*1e12/cc*(d*nsample-d)))-H(1,b123);
    p(1,b123)=fsolve(H1,a(1,b123));
end
n1=real(p);
k1=-imag(p);
    
figure(11)
plot(f_data_trunc(200:1:1000),n1(200:1:1000),'b','linewidth',2)
hold on
plot(f_data_trunc(200:1:1000),k1(200:1:1000),'r','linewidth',2)
grid on