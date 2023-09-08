%CODE FOR SINGLE SIDE BAND SUPRESSED CAREER(SSB-SC)


clc;
clear all;
close all;

%INPUTS
Ac=input('Enter the amplitude of carrier signal :');
fc=input('Enter the frequency of carrier signal :');
Am=input('Enter the amplitude of message signal :');
fm=input('Enter the frequency of message signal :');

interval=(1/(20*fc));%interval steps
end_t=0.05;%ending time
t=0:interval:end_t; %time(fs>>2fm) by sampling theorem 

%signals
m_t=Am*cos(2*pi*fm*t);%message signal
c_t=Ac*cos(2*pi*fc*t);%career signal
s_t=(Am*Ac/2)*cos(2*pi*(fc+fm)*t);%MODULATED SIGNAL


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