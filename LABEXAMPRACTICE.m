% clc;
% clear all;
% close all;
% 
% A_m=25;
% 
% f_m=3;
% 
% f_s=200*f_m;
% 
% T_s=1/f_s;
% 
% t=0:T_s:1-T_s;
% 
% m_t=A_m*sin(2*pi*f_m*t);
% 
% t1=0:(2*T_s):1-T_s;
% 
% sample=A_m*sin(2*pi*f_m*t1);
% 
% delta=5;
% %delta_modulated_signal=zeros(size(sample));
% delta_modulated_signal(1)=0;
% quantized_signal(1)=0;
% step_size(1)=delta;
% d_error(1)=1;
% 
% for i=2:length(sample)
%     if sample(i)>delta_modulated_signal(i-1)
%         d_error(i)=+1;
%         step_size(i)=d_error(i)*abs(step_size(i-1))+d_error(i-1)*delta;
%         delta_modulated_signal(i)=delta_modulated_signal(i-1)+step_size(i);
%         quantized_signal(i)=1;
%     else
% 
%         d_error(i)=-1;
%         step_size(i)=d_error(i)*abs(step_size(i-1))+d_error(i-1)*delta;
%         delta_modulated_signal(i)=delta_modulated_signal(i-1)+step_size(i);
%         quantized_signal(i)=0;
%     end
% end
% figure(1)
% subplot(4,1,1)
% plot(m_t)
% 
% subplot(4,1,2)
% stem(sample)
% 
% subplot(4,1,3)
% stairs(delta_modulated_signal)
% hold on
% plot(sample)
% hold off
% 
% 
% subplot(4,1,4)
% stem(quantized_signal)
% 
% 
% 
% demodulated_signal(1)=0;
% d_error(1)=1;
% d_step_size(1)=delta;
% 
% for i=2:length(quantized_signal)
%     if quantized_signal(i)==1
%         d_error(i)=+1;
%         d_step_size(i)=d_error(i)*abs(d_step_size(i-1))+d_error(i-1)*delta;
%         demodulated_signal(i)=demodulated_signal(i-1)+d_step_size(i);
%     else
%         d_error(i)=-1;
%         d_step_size(i)=d_error(i)*abs(d_step_size(i-1))+d_error(i-1)*delta;
%         demodulated_signal(i)=demodulated_signal(i-1)+d_step_size(i);
%     end
% end
% 
% figure(2)
% subplot(4,1,1)
% stairs(demodulated_signal)
% 
% subplot(4,1,2)
% plot(demodulated_signal)
% 
% filtered_signal=lowpass(demodulated_signal,f_m,1/(2*T_s));
% subplot(4,1,3)
% plot(filtered_signal)
% hold on
% plot(sample)
% hold off
% legend('filtered','sampled ')


clc;
clear all;
close all;

A_m=2;
f_m=1;
f_s=200*f_m;
T_s=1/f_s;
t=0:T_s:1-T_s;
m_t=A_m*sin(2*pi*f_m*t);
L=8;
dell=2*A_m/(L-1);
interval=dell/2;
quantization_levels=A_m:-dell:-A_m;
threshold_level_1=A_m-dell/2;
threshold_levels=threshold_level_1:-dell:-A_m;
t2=0:2*T_s:1-T_s;
sample=A_m*sin(2*pi*f_m*t2);

figure(1)
subplot(4,1,1)
plot(t,m_t);
title('message Signal');
grid on;

subplot(4,1,2)
stem(t2,sample);
title('sampled message Signal');
grid on;

encoded_quantization_levels=0:L-1;
quantized_sample=[];
for i=1:length(sample)
    for j=1:length(threshold_levels)
        if(sample(i)>=threshold_levels(j)&&j==1)
            quantized_sample(i)=quantization_levels(j);
        elseif(sample(i)<=threshold_levels(j)&&j==length(threshold_levels))
            quantized_sample(i)=quantization_levels(j);
