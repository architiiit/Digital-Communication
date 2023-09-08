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

%career wave inphase component
phi_t1=sqrt(2/T_b)*cos(2*pi*fc*t);


%career wave quadrature component
phi_t2=sqrt(2/T_b)*sin(2*pi*fc*t);

even_bits=[];
odd_bits=[];
for j=1:1:length(bit_pattern)
    if(mod(j,2)==0)
       even_bits=[even_bits bit_pattern(j)];
    else
        odd_bits=[odd_bits bit_pattern(j)];
    end
end

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

EPNRZ=[]; %PNRZ storing variable of message signal of even bits

% POLAR NRZ creation
for j=1:1:length(even_bits)
    if even_bits(j)==1
        x_t=repmat(even_bits(j),1,2*T_b);   
    else
        x_t=repmat(-1,1,2*T_b);
    end
    EPNRZ=cat(2,EPNRZ,x_t);
end


OPNRZ=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(odd_bits)
    if odd_bits(j)==1
        x_t=repmat(odd_bits(j),1,2*T_b);   
    else
        x_t=repmat(-1,1,2*T_b);
    end
    OPNRZ=cat(2,OPNRZ,x_t);
end

s_t1=OPNRZ.*phi_t1;
s_t2=EPNRZ.*phi_t2;

figure(1)

subplot(3,2,1)
plot(PNRZ)
title('MESSAGE SIGNAL BITS')

grid on;
subplot(3,2,2);
plot(EPNRZ);
title('Even POLAR NRZ Signal')
grid on;

subplot(3,2,3);
plot(OPNRZ);
title('Odd POLAR NRZ Signal')
grid on;

subplot(3,2,4);
plot(s_t1);
title('s_t1')
grid on;

subplot(3,2,5);
plot(s_t2);
title('s_t2')
grid on;

s_t=s_t1+s_t2;

subplot(3,2,6);
plot(s_t);
title('s_t')
grid on;


%%DEMODULATION
y1_t=s_t.*phi_t1;
y2_t=s_t.*phi_t2;

figure(2)
subplot(2,1,1)
plot(y1_t)
title('PRODUCT MODULATOR 1')
grid on;

subplot(2,1,2)
plot(y2_t)
title('PRODUCT MODULATOR 2')
grid on;
final_bit=[];

for j = 1:2:length(bit_pattern)

    % calculate the start and end index for integration interval
    start_index = (j-1)*T_b+1;
    end_index = j*T_b;

    % integrator
    bit_sum1 = sum(y1_t(start_index:end_index));
    bit_sum2 = sum(y2_t(start_index:end_index));
    
    if(bit_sum1>0 && bit_sum2>0)
        final_bit=[final_bit 1 1];
    elseif(bit_sum1>0 && bit_sum2<0)
        final_bit=[final_bit 1 0];
    elseif(bit_sum1<0 && bit_sum2>0)
        final_bit=[final_bit 0 1];
    elseif(bit_sum1<0 && bit_sum2<0)
        final_bit=[final_bit 0 0];
    end
end

FPNRZ=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(final_bit)
    if final_bit(j)==1
        x_t=repmat(final_bit(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    FPNRZ=cat(2,FPNRZ,x_t);
end

figure(3)
subplot(2,1,1)
plot(FPNRZ)
title('FINAL RECOVERED BITS IN PNRZ')
grid on;

subplot(2,1,2)
plot(PNRZ)
title('MESSAGE SIGNAL BITS')
grid on;

disp('output signal bits');
disp(final_bit);

%CHECK ERROR IN SIGNAL
error=biterr(bit_pattern,final_bit);
if error==0
    disp('NO BIT ERROR');
else
    disp('BIT ERROR EXIST');
end