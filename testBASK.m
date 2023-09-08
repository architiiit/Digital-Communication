%% This code is about BASK implementation in matlab

%Name:Archit Vashist
%Roll no.2021BEC0022

clc;
clear all;
close all;

%bit generation
b_k=round(rand(1,10));
disp('Input signal bits')
disp(b_k)
%career frequency
f_c=10;

%sampling frequency
f_s=20*f_c;

 %sampling time period
T_s=1/f_s;

%time interval for sampling
t=0:T_s:1-T_s;

%bit length
T_b=length(t)/length(b_k);

%career wave
c_t=sqrt(2/T_b)*cos(2*pi*f_c*t);

UNRZ=[];%UNRZ storing variable

%uniploal NRZ creation
for j=1:1:length(b_k)
    if b_k(j)==1
        a=repmat(b_k(j),1,T_b);
    else
        a=repmat(0,1,T_b);
    end
    UNRZ=[UNRZ a];
end



%modulated signal
s_t=c_t.*UNRZ;

figure(1);
%UNRZ signal
subplot(3,1,1)
plot(UNRZ);
title('unipolar NRZ signal')
xlabel('time')
ylabel('Amplitude')
grid on;

%CARRRER SIGNAL PLOT
subplot(3,1,2)
plot(c_t);
title('career signal')
xlabel('time')
ylabel('Amplitude')
grid on;

%MODULATED SIGNAL PLOT
subplot(3,1,3)
plot(s_t);
title('modulated signal')
xlabel('time')
ylabel('Amplitude')
grid on;


%%DEMODULATION STEPS
u_t=s_t.*c_t;
figure(2);
subplot(3,1,1);
plot(u_t);
title('Product modulator output')
xlabel('time')
ylabel('Amplitude')
grid on;


bit_s = []; % recovered bits

for i = 1:length(b_k)
    % calculate the start and end index of the current bit length in y_t
    start = (i-1)*T_b+1;
    ending = i*T_b;

    % sum the values of y_t over the current bit length
    b = sum(u_t(start:ending));
   bit_s=[bit_s b];
end

final_bits=[];

for i=1:length(b_k)
if bit_s(i)>0.5
    final_bits=[final_bits 1];
else
    final_bits=[final_bits 0];
end
end



%display recovered bits
disp('Recovered bits:')
disp(final_bits)


UNRZ1=[];%UNRZ REPRESENTATION OF RECOVERED BITS
for j=1:1:length(b_k)
    if final_bits(j)==1
        a=repmat(final_bits(j),1,T_b);
    else
        a=repmat(0,1,T_b);
    end
   UNRZ1 =[UNRZ1 a];
end

%PLOTTING OF RECOVERED BITS
subplot(3,1,2);
plot(UNRZ1);
title('RECOVERED BITS UNRZ')
xlabel('time')
ylabel('Amplitude')
grid on;

if biterr(b_k,final_bits)
    disp('BIT ERROR PLEASE CHECK THE CODE ')
else
    disp('NO BIT ERROR')
end