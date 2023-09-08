%NAME:ARCHIT VASHIST
%ROLL NO:2021BEC0022


%% BPSK IMPLEMENTATION IN MATLAB
clc;
clear all;
close all;


%bit generation
bit_pattern=round(rand(1,10));
disp('Input signal bits');
disp(bit_pattern);

% career frequency
fc=10;

%sampling frequency
fs=20*fc;

%sampling time
Ts=1/fs;

%time interval for sampling
t=0:Ts:1-Ts;

%bit length
T_b=length(t)/length(bit_pattern);

%career wave
phi_t=sqrt(2/T_b)*cos(2*pi*fc*t);

PNRZ=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(bit_pattern)
    if bit_pattern(j)==1
        x_t=repmat(bit_pattern(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ=cat(2,PNRZ,x_t);
end

figure(1)
%UNRZ signal
subplot(3,1,1);
plot(PNRZ);
title('POLAR NRZ Signal')
grid on;


%plot career signal
subplot(3,1,2);
plot(phi_t);
title('Career Signal');
grid on;


%modulated signal
s_t=phi_t.*PNRZ;


%plot modulated signal
subplot(3,1,3);
plot(s_t);
title('Modulated Signal');
grid on;


%DEMODULATION SECTION

y_t=s_t.*phi_t; %STEP1->product modulator

figure(2)
subplot(3,1,1);
plot(y_t);
title('Demodulation Step 1');
grid on;


%INTEGRATOR

bit_s = []; %integrator values

for j = 1:length(bit_pattern)

    % calculate the start and end index for integration interval
    start_index = (j-1)*T_b+1;
    end_index = j*T_b;

    % integrator
    bit_sum = sum(y_t(start_index:end_index));
    
    % storing integrator values in an array for message signal creation
    bit_s = [bit_s bit_sum];

end


%DECISION DEVICE
bits=[];
for j=1:1:length(bit_pattern)
    if bit_s(j)>0 %THRESHOLD VALUE
        a=1;
    else
        a=0;
    end
   bits =[bits a];
end

disp('output signal bits');
disp(bits);


PNRZ_R=[];%POLAR NON RETURN TO ZERO REPRESENTATION OF RECOVERED SIGNAL

for j=1:1:length(bits)
    if bit_pattern(j)==1
        x_t=repmat(bits(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_R=cat(2,PNRZ_R,x_t);
end


%BIPOLAR NRZ OF RECOVERED SIGNAL
subplot(3,1,2);
plot(PNRZ_R);
title('recovered bits after demodulation(recovered signal)');
grid on;

%MESSAGE SIGNAL
subplot(3,1,3);
plot(PNRZ);
title('actual bits (message signal)');
grid on;


%CHECK ERROR IN SIGNAL
error=biterr(bit_pattern,bits);
if error==0
    disp('NO BIT ERROR');
else
    disp('BIT ERROR EXIST');
end

