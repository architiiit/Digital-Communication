%NAME:ARCHIT VASHIST
%ROLL NO:2021BEC0022

%% DPSK IMPLEMENTATION IN MATLAB
clc;
clear all;
close all;

%bit generation
b_k=round(rand(1,9));
disp('Input signal bits');
disp(b_k);


% career frequency
fc=10;

%sampling frequency
fs=20*fc;

%sampling time
Ts=1/fs;

%time interval for sampling
t=0:Ts:1-Ts;

%bit length
T_b=length(t)/(length(b_k)+1);

%career signal
phi_t=sqrt(2/T_b)*cos(2*pi*fc*t);


PNRZ_bk=[]; %PNRZ storing variable of message signal
%POLAR NRZ creation
for j=1:1:length(b_k)
    if b_k(j)==1
        x_t=repmat(b_k(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_bk=cat(2,PNRZ_bk,x_t);
end
figure(1)
subplot(4,1,1);
plot(PNRZ_bk);
title('POLAR NRZ Signal of actual bits')
grid on;
d_km=[1];%delayed bits
%d_k=[];%actual bits to be transmitted


%perform logical network
for j=1:1:length(b_k)
    a=not(xor(b_k(j),d_km(j)));
    d_km=[d_km a];
end

%for j=2:1:length(d_km)
%   d_k=[d_k d_km(j)];
%end

PNRZ=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(d_km)
    if d_km(j)==1
        x_t=repmat(d_km(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ=cat(2,PNRZ,x_t);
end




%modulated signal
s_t=phi_t.*PNRZ;


subplot(4,1,2);
plot(PNRZ);
title('POLAR NRZ Signal of xnor signal d_k-1')
grid on;


%plot career signal
subplot(4,1,3);
plot(phi_t);
title('Career Signal');
grid on;


%plot modulated signal
subplot(4,1,4);
plot(s_t);
title('Modulated Signal');
grid on;



%%DEMODULATION

y_t=s_t.*phi_t; %STEP1->product modulator

figure(2)
subplot(4,1,1);
plot(y_t);
title('PRODUCT MODULATOR');
grid on;

%INTEGRATOR

bit_s = []; %integrator values

for j = 1:length(d_km)

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
for j=1:1:length(d_km)
    if bit_s(j)>0 %THRESHOLD VALUE
        a=1;
    else
        a=0;
    end
   bits =[bits a];
end



final_bits=[];
for k=1:1:length(bits)-1
    b=not(xor(bits(k),bits(k+1)));
    final_bits=[final_bits b];
end


PNRZ_bits=[]; %PNRZ storing variable of message signal
%POLAR NRZ creation
for j=1:1:length(bits)
    if bits(j)==1
        x_t=repmat(bits(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_bits=cat(2,PNRZ_bits,x_t);
end
subplot(4,1,2);
plot(PNRZ_bits);
title('recovered bits after integration');
grid on;



PNRZ_final_bits=[]; %PNRZ storing variable of message signal
%POLAR NRZ creation
for j=1:1:length(final_bits)
    if final_bits(j)==1
        x_t=repmat(final_bits(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_final_bits=cat(2,PNRZ_final_bits,x_t);
end
subplot(4,1,3);
plot(PNRZ_final_bits);
title('recovered bits after xnor logic');
grid on;

subplot(4,1,4);
plot(PNRZ_bk);
title('actual sent bits');
grid on;

disp('output signal bits');
disp(final_bits);

%CHECK ERROR IN SIGNAL
error=biterr(b_k,final_bits);
if error==0
    disp('NO BIT ERROR');
else
    disp('BIT ERROR EXIST');
end