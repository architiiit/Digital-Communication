%NAME-ARCHIT VASHIST
%ROLL -2021BEC0022
%DELTA MODULATION IN MATLAB
clc;
clear all;
close all;
%message frequency
fm=1;
%message amplitude
A_m=2;
%sampling frequency higher
fs1=200*fm;
%sampling time
Ts1=1/fs1;
%time interval
t1=0:Ts1:1-Ts1;
%message signal
m_t=A_m*sin(2*pi*fm*t1);
%step size initialisation
step_size=0.4;
%Higher sampling time for taking samples from the given message signal
t2=0:2*Ts1:1-Ts1;
%samples
sample=A_m*sin(2*pi*fm*t2);
figure(1)
subplot(4,1,1)
plot(m_t);
title('message Signal');
grid on;
subplot(4,1,2)
stem(sample);
title('sampled message Signal');
grid on;
% Initialize delta modulation output signal
delta_modulated_signal = zeros(size(sample));
% Initialize delta modulation output signal
delta_modulated_signal(1) = sample(1); % Set the first sample
quantized_signal=[];
% Perform delta modulation
for i = 2:length(sample)
 if sample(i) >=delta_modulated_signal(i-1)%if sample is greater then add step size
 delta_modulated_signal(i) = delta_modulated_signal(i-1) + step_size;
 quantized_signal(i)=1;
 else %if sample is less than quantized level then decrease
 delta_modulated_signal(i) = delta_modulated_signal(i-1) - step_size;
 quantized_signal(i)=0;
 end
end
% Plot the delta-modulated signal
subplot(4,1,3)
stairs(delta_modulated_signal);
title('Delta Modulated Signal');
hold on;
plot(sample);
hold off;
grid on;
subplot(4,1,4)
stem(quantized_signal);
title('Quantized signal');
grid on;
%%DEMODULATION
% Initialize the demodulated signal
demodulated_signal = zeros(size(delta_modulated_signal));
demodulated_signal(1)=0;
% Perform delta demodulation
for i = 2:length(quantized_signal)
 
 if quantized_signal(i) == 1 %means increase by step size
 demodulated_signal(i) = demodulated_signal(i-1) + step_size;
 else %means decrease by step size
 demodulated_signal(i) = demodulated_signal(i-1) - step_size;
 end
end
% Plot the demodulated signal
figure(2);
subplot(3,1,1);
stairs(demodulated_signal);
title('Demodulated Signal');
xlabel('Time');
ylabel('Amplitude');
grid on;
%demodulated signal
subplot(3,1,2);
plot(demodulated_signal);
title('Demodulated Signal');
xlabel('Time');
ylabel('Amplitude');
grid on;
%low pass filter 
filtered_signal=lowpass(demodulated_signal,fm,(1/(2*Ts1)));
% Plot the filtered signal
subplot(3,1,3);
plot(filtered_signal);
title('Filtered Signal');
xlabel('Time');
ylabel('Amplitude');
grid on;
hold on;
plot(sample);
title('Original Message Signal');
xlabel('Time');
ylabel('Amplitude');
legend('recovered signal','original message signal');
hold off
grid on;