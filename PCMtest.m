%%PULSE CODE MODULATION
%NAME:ARCHIT VASHIST
%ROLL:2021BEC0022

clc;
clear all;
close all;

%message frequency
fm=2;

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

%quantization levels
L=8;

%step size
dell=(2*A_m)/(L-1);

%deta by 2
interval=dell/2;

%partition
quantization_levels = A_m:-dell:-A_m;

%set threshold level
threshold_level_1=A_m-(dell/2);

%set other threshold levels
threshold_levels=threshold_level_1 :-dell:-A_m;

%Higher sampling time for taking samples from the given message signal
t2=0:2*Ts1:1-Ts1;

%samples
sample=A_m*sin(2*pi*fm*t2);

figure(1)
subplot(4,1,1)
plot(t1,m_t);
title('message Signal');
grid on;


subplot(4,1,2)
stem(t2,sample);
title('sampled message Signal');
grid on;

%setting quantization levels from zero to seven
encoded_quantization_levels=0:L-1;


% quantized samples using for loop
quantized_sample = zeros(1, length(sample));

for i = 1:length(sample)
    for j = 1:length(threshold_levels)
        
        if(sample(i) >= threshold_levels(j) && j==1)
            quantized_sample(i) = quantization_levels(j);
            encoded_signal(i)=encoded_quantization_levels(j);

        elseif(sample(i) <= threshold_levels(j) && j==length(threshold_levels))

            quantized_sample(i) = quantization_levels(j+1);
            encoded_signal(i)=encoded_quantization_levels(j+1);

        elseif(sample(i) >= threshold_levels(j) && sample(i) <= threshold_levels(j-1))

            quantized_sample(i) = quantization_levels(j);
            encoded_signal(i)=encoded_quantization_levels(j);
        end
    end
end

subplot(4,1,3)
stem(t2,quantized_sample);
yline(quantization_levels,'-');
title('Quantized Sampled Signal');
grid on


subplot(4,1,4)
stairs(t2, quantized_sample);
title('Encoded Signal');
grid on


figure(2)
% plot the signals overlapping
subplot(2,1,1)
hold on;
plot(t2,sample,'LineWidth',1.5);
yline(quantization_levels,'-.');
stairs(t2, quantized_sample, ...
    'LineWidth',1.5);
hold off;
legend('Sampled Message Signal','Encoded Signal');
title('Sampled Message Signal and Encoded Signal');
grid on;


%error in both signals
error=quantized_sample-sample;

subplot(2,1,2);
plot(t2,error);
title('error');
grid on;


disp('sum of error is :');
disp(sum(error));

figure(3)
stem(encoded_signal);
xlabel('samples');
ylabel('quantization levels');

coded_signal=[];
for i=1:length(encoded_signal)
code=dec2bin(encoded_signal(i),3);
coded_signal=[coded_signal code];
end
final_codes=[];
for j=1:length(coded_signal)
coded=str2num(coded_signal(j));
final_codes=[final_codes coded];
end

figure(4)
stem(final_codes);
xlabel('PULSE CODED BIT PATTERN');