%CODE FOR DOUBLE SIDE BAND FULL CAREER(DSB-FC)

clc;
clear all;
close all;


%INPUTS
Ac=input('Enter the amplitude of carrier signal :');
fc=input('Enter the frequency of carrier signal :');
Am=input('Enter the amplitude of message signal :');
fm=input('Enter the frequency of message signal :');
k_a=input('Enter the value of modulation sensitivity:');

interval=(1/(20*fc));%interval steps
end_t=0.050;%ending time
t=0:interval:end_t; %time(fs>>2fm) by sampling theorem 

if(abs(k_a*Am)>1)
 disp("Phase reversal will take place Over Modulation");
elseif(abs(k_a*Am)<1)
 disp("No Phase reversal Under Modulation");
else
 disp("Critical Modulation");
end


%signals
m_t=Am*cos(2*pi*fm*t);%message signal
c_t=Ac*cos(2*pi*fc*t);%career signal
%s_t=Ac*(1+(Am/Ac)*cos(2*pi*fm*t)).*cos(2*pi*fc*t);%MODULATED SIGNAL
s_t=c_t.*(1+(k_a*m_t));


%PLOTS

figure(1)
plot(t,s_t)
grid on
xlabel('time(sec)')
ylabel('Amplitude(V)')
title('MODULATED SIGNAL')


figure(2)
plot(t,c_t)
grid on
xlabel('time(sec)')
ylabel('Amplitude(V)')
title('CARRIER SIGNAL')


figure(3)
plot(t,m_t)
grid on
xlabel('time(sec)')
ylabel('Amplitude(V)')
title('MESSAGE SIGNAL')