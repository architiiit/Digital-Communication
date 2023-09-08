%% This code is about BASK implementation in matlab
clc;
clear all;
close all;


%generate binary bits
%bit_pattern=round(rand(1,5));
%Ts=1/fs; sampling time
%fs>=2fm or fc; fails at fm=1 means fs=2;Ts=0.5;
%Ts=1/fs;fs=32fc;fc=100;
%t=0:Ts:100;
x=round(rand(1,5));%bit pattern as input
fc=100;%career frequency
t=0:1/(32*fc):0.1;%sampling time
phi_t=sqrt(2/1)*cos(2*pi*fc*t);%career signal
stem(phi_t)
plot(phi_t)
for j=1:1:length(x)
    if x(j)==1
        x_t=repmat(x(j),1,length(t))/length(x);
    else
        x_t=repmat(x(j),1,length(t))/length(x);
    end
    x1_t=cat(2,x1_t,x_t);
end

%plot(x1_t);
s_t=phi_t.*x1_t;
plot(s_t);

%demodulation process
y_t=phi_t.*s_t;
