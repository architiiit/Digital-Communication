clc;
clear all;
close all;

fc=100;
fm=fc/10;
fs=100*fc;
Ts=1/fs;
Tm=1/fm;
t=0:Ts:1-Ts;
m_t=cos(2*pi*fm*t);
c_t=0.5*square(2*pi*fc*t)+0.5;
s_t=m_t.*c_t;

t_t=[];
for i=1:length(s_t)
    if(s_t(i)==0)
        t_t=[t_t s_t(i)];
    else
        t_t=[t_t s_t(i)+1];
    end
end
figure(1)
subplot(4,1,1);
plot(t,m_t);
title('message signal');
xlabel('timeperiod');
ylabel('amplitude');
grid on;

subplot(4,1,2);
plot(t,c_t);
title('carrier signal');
xlabel('timeperiod');
ylabel('amplitude');
grid on;

subplot(4,1,3);
plot(t,s_t);
xlabel('timeperiod');
ylabel('amplitude');
title('modulated signal');
grid on;

subplot(4,1,4);
plot(t,t_t);
title('single side band modulated signal');
xlabel('timeperiod');
ylabel('amplitude');
grid on;


d_m=s_t.*c_t;
figure(2);
subplot(3,1,1);
plot(t,d_m);
title('demodulated signal');
xlabel('timeperiod');
ylabel('amplitude');
grid on;


