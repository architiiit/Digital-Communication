%% This code is about BFSK implementation in matlab

%Name:Archit Vashist
%Roll no.2021BEC0022

clc;
clear all;
close all;


%bit generation
bit_pattern=round(rand(1,5));
disp('Input signal bits');
disp(bit_pattern);

%first career frequency
fc1=10;

%second career frequency
fc2=5;

fs=20*fc1;
Ts=1/fs;

%time interval for sampling
t=0:Ts:1-Ts;

%bit length
T_b=length(t)/length(bit_pattern);

%first career wave
phi_t1=sqrt(2/T_b)*cos(2*pi*fc1*t);

%second career wave
phi_t2=sqrt(2/T_b)*cos(2*pi*fc2*t);

UNRZ=[];%UNRZ storing variable

%uniploal NRZ creation
for j=1:1:length(bit_pattern)
    if bit_pattern(j)==1
        x_t=repmat(bit_pattern(j),1,T_b);
    else
        x_t=repmat(bit_pattern(j),1,T_b);
    end
    UNRZ=cat(2,UNRZ,x_t);
end

IUNRZ=[];%UNRZ storing variable

%INVERTED uniploal NRZ creation
for j=1:1:length(bit_pattern)
    if bit_pattern(j)==1
        x1_t=repmat(0,1,T_b);
    else
        x1_t=repmat(1,1,T_b);
    end
    IUNRZ=cat(2,IUNRZ,x1_t);
end

s_t1=phi_t1.*UNRZ;
s_t2=phi_t2.*IUNRZ;
s_t=s_t1+s_t2;

figure(1)
%UNRZ signal
subplot(4,1,1)
plot(UNRZ);
title('unipolar NRZ signal')
grid on;

%Inverted NRZ
subplot(4,1,2)
plot(IUNRZ);
title('inverted NRZ signal')
grid on;

%MODULATED SIGNAL PLOT
subplot(4,1,3)
plot(s_t);
title('modulated signal')
grid on;


%FIRST CAREER SIGNAL
figure(2);
subplot(3,1,1);
plot(phi_t1);
title('career 1')
grid on;

%SECOND CAREER SIGNAL
subplot(3,1,2);
plot(phi_t2);
title('career 2')
grid on;



%%DEMODULATION STEPS
x_t1=s_t.*phi_t1;
x_t2=s_t.*phi_t2;


figure(3);
subplot(3,1,1);
plot(x_t1);
title('x1')
grid on;

%SECOND CAREER SIGNAL
subplot(3,1,2);
plot(x_t2);
title('x2')
grid on;

bits = []; % recovered bits

for j = 1:length(bit_pattern)
    % calculate the start and end index of the current bit length in y_t
    start_index = (j-1)*T_b+1;
    end_index = j*T_b;

    % sum the values of y_t over the current bit length
    bit_sum1 = sum(x_t1(start_index:end_index));
    bit_sum2 = sum(x_t2(start_index:end_index));

        % determine the bit value based on the threshold
    if bit_sum1 > bit_sum2
        bit_value = 1;
    else
        bit_value = 0;
    end
     bits = [bits bit_value];
end

disp('Recovered bits:')
disp(bits)


figure(4)
subplot(2,1,1);
plot(bits);
title('RECOVERED BITS (AFTER INTEGRATION)');
grid on;

BITS=[];%BINARY REPRESENTATION OF RECOVERED BITS
for j=1:1:length(bits)
    if bits(j)==1
        a=repmat(bits(j),1,T_b);
    else
        a=repmat(bits(j),1,T_b);
    end
   BITS =cat(2,BITS,a);
end


subplot(2,1,2);
plot(BITS);
title('RECOVERED BITS');
grid on;
biterr(bit_pattern,bits)