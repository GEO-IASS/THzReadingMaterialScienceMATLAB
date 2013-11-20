 clear all;
 clc;

[t_sample,E_sample] = open_picotd('Silicon.picotd',4101);
[t_reference,E_reference] = open_picotd('Reference.picotd',4101);


%time to reference 
%%t_silicon_trunc = 0:0.078125:(4095*0.078125);
%trunc and padding
width = 1000;  
c = 3e8;
d=690e-6;
% the cutting parameter
[x,peak] = min(E_sample);                      %finds the minimum peak in the chart
[x,peak_reference] = min(E_reference);         %finds the minimum peak in the chart
delta_t = t_sample(peak) - t_sample(peak_reference);        %finds the time difference in the peaks
delta_t = delta_t*10^-12;    
start = peak - width/2;
stop = peak + width/2;
start_reference = peak_reference - width/2;
stop_reference = peak_reference + width/2;

E_sample_trunc = E_sample(start:stop);
E_reference_trunc = E_reference(start_reference:stop_reference);

t_sample_trunc = t_sample(start:stop);
t_reference_trunc = t_reference(start_reference:stop_reference);
%E vs t plot
figure(1)
plot(t_reference(1:4096),E_reference,'b','Linewidth',1);
hold on;
plot(t_sample(1:4096),E_sample,'r','Linewidth',1);
grid on;
%E vs t truncated plot
figure(2);
plot(t_reference_trunc,E_reference_trunc,'b','Linewidth',2);
hold on;
plot(t_sample_trunc,E_sample_trunc,'r','Linewidth',2);
grid on;

%adding zeros

E_reference_pad(1:4096,1)=[0];
E_reference_pad(start_reference:stop_reference,1) = E_reference(start_reference:stop_reference,1);
E_sample_pad(1:4096,1)=[0];
E_sample_pad(start:stop,1) = E_sample(start:stop);

Ts = 0.078125e-12;       %time period
Fs = 1/Ts;          % Sampling frequency interval
L = 4096 ;         % number of sample
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
F_sample = fft(E_sample_pad,NFFT*2);         % fft of sample
F_sample = F_sample(1:NFFT); 
F_reference = fft(E_reference_pad,NFFT*2);    %fft of reference
F_reference = F_reference(1:NFFT);
f = ((Fs/2)*linspace(0,1,NFFT))'*1e-12;         %trapose to keep in same direction     

%{
%ffts and phase
F_sample = fft(E_sample_pad); 
F_reference = fft(E_reference_pad);
%}
H_sample = F_sample./F_reference;

% confining the frequency
freq_max=1;
while db(F_reference(freq_max))>-40
      freq_max=freq_max+1;
end

F_sample_trunc = F_sample(1:freq_max);
F_reference_trunc = F_reference(1:freq_max);
H_trunc = H_sample(1:freq_max);

% unwraping the phase
p_raw=unwrap(phase(E_sample_pad));
p_ref=unwrap(phase(E_reference_pad));
p_norm=unwrap(phase(H_sample));

%plot for transfer function,phase and fft
figure(4)
plot(f(50:freq_max),abs(H_trunc(50:freq_max)),'b','Linewidth',2);
 xlabel('frequency(THz)','FontWeight','bold');
 ylabel('Abs(H)','FontWeight','bold');
 title('magnitude for H','FontWeight','bold');
 figure(5)
plot(f(120:freq_max),p_norm(120:freq_max),'b');
 xlabel('frequncey (THz)','FontWeight','bold');
 ylabel('Phase','FontWeight','bold');
 title('Phase','FontWeight','bold');
 figure(6)
plot(f(50:freq_max), db(F_reference_trunc(50:freq_max)),'b','Linewidth',2);
 hold on;
plot(f(50:freq_max),db(F_sample_trunc(50:freq_max)),'r','Linewidth',2);
  title('FFT of reference and Sample','FontWeight','bold');
  legend('reference', 'sample');
  
   H_sample = H_sample';
  f = f';
  

 %finding n,k
 n = 1+c/(2*pi*d*1e12)*phase(H_sample)./f;
 alpha = -2/d*log(abs(H_sample).*(n+1).*(n+1)./n/4);
 k = alpha*c/4/1e12/pi./f;
 epsilon_real = n.*n-k.*k;
 epsilon_imag = 2*n.*k;
 
 figure(7)
 plot(f(20:1:1000),n(20:1:1000),'b','Linewidth',2);
 xlabel('frequency(THz)','FontWeight','bold');
 ylabel('n','FontWeight','bold');
 
 figure(8)
 plot(f(200:1000), k(200:1000), 'r','Linewidth',2)
 xlabel('frequency(THz)','FontWeight','bold')
 ylabel('k','FontWeight','bold')
%plotting eps1 and eps 2
figure(8)
 plot(f(200:1:1000),epsilon_real(200:1:1000),'b','Linewidth',2)
 xlabel('frequency(THz)','FontWeight','bold')
 ylabel('epsilon 1','FontWeight','bold')
 figure(9)
 plot(f(200:1000), epsilon_imag(200:1000), 'r','Linewidth',2)
 xlabel('frequency(THz)','FontWeight','bold')
 ylabel('epsilon 2','FontWeight','bold')
 
  n_inital=c*delta_t/(d-c*delta_t);
 %
 conductivity =n.*k.*f;
 figure(10)
 plot(f(200:1:1000),conductivity(200:1:1000),'b','Linewidth',2)
 xlabel('Freq (THz)','FontWeight','bold')
 ylabel('Conductivity','FontWeight','bold')
 grid on
 

%Solving for epsilon

x0_sample = n + 1i*k ;                           

for i=130:1000
    options = optimset('Display','off');  % Turn off display
    func=@(x) 4*x*(exp(-1i*2*pi*f(1,i)*1e12*d*(x-1)/c))*(1/(1+x)^2)-H_sample(1,i);
    %(4*x/(1+x)^2*exp(-1i*2*pi*f(1,i)*1e12/c*(d*x-d)))-H_sample(1,i);
    p_sample(1,i) = fsolve(func, x0_sample(1,i),options); 
end


n1=real(p_sample);
k1=-imag(p_sample);
    
figure(11)
plot(f(200:1:1000),n1(200:1:1000),'b','linewidth',2)
hold on
plot(f(200:1:1000),k1(200:1:1000),'r','linewidth',2)
grid on

