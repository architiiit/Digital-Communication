%% This code is about BFSK implementation in matlab

%Name:Archit Vashist
%Roll no: 2021BEC0022

clc;
clear all;
close all;


%bit generation
bit_pattern=round(rand(1,5));
disp('Input signal bits');
disp(bit_pattern);

%higher career frequency
fc1=10;

%lower career frequency
fc2=5;

%sampling frequency
fs=20*fc1;

%sampling time
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

IUNRZ=[];%Inverted UNRZ storing variable

%INVERTED uniploal NRZ creation
for j=1:1:length(bit_pattern)
    if bit_pattern(j)==1
        x1_t=zeros(1,T_b);
    else
        x1_t=ones(1,T_b);
    end
    IUNRZ=cat(2,IUNRZ,x1_t);
end

s_t1=phi_t1.*UNRZ;%product modulator with phi_1
s_t2=phi_t2.*IUNRZ;%product modulator with phi_2
s_t=s_t1+s_t2;%modulated signal

figure(1)
%UNRZ signal
subplot(4,1,1)
plot(UNRZ);
title('unipolar NRZ signal')
grid on;

%Inverted UNRZ
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
x_t1=s_t.*phi_t1;%upper side product modulator
x_t2=s_t.*phi_t2;%lower side product modulator


figure(3);
%x1 plot
subplot(3,1,1);
plot(x_t1);
title('x1')
grid on;

%x2 plot
subplot(3,1,2);
plot(x_t2);
title('x2')
grid on;


bit_s1 = []; % upper side integrator
bit_s2 = []; % lower side integrator 

for j = 1:length(bit_pattern)
    % calculate the start and end index 
    start_index = (j-1)*T_b+1;
    end_index = j*T_b;

    % integrator for x1 and x2
    bit_sum1 = sum(x_t1(start_index:end_index));
    bit_sum2 = sum(x_t2(start_index:end_index));
    
    % integrator values in an array
    bit_s1 = [bit_s1 bit_sum1];
    bit_s2 = [bit_s2 bit_sum2]; 

end

figure(4)
%subplot(3,1,1);
%plot(bit_s1);
%title('RECOVERED BITS (AFTER INTEGRATION) UPPER HALF');
%grid on;


%subplot(3,1,2);
%plot(bit_s2);
%title('RECOVERED BITS (AFTER INTEGRATION) LOWER HALF');
%grid on;

bits=[];  
%RECOVERED BITS STORAGE (DECISION DEVICE)
for j=1:1:length(bit_pattern)
    if bit_s1(j)>bit_s2(j) %COMPARISON 
        a=1;
    else
        a=0;
    end
   bits =[bits a];
end

disp('output signal bits');
disp(bits);


BITS=[];
%REPRESENTATION OF RECOVERED BITS IN UNRZ(LINE CODED FORM)
for j=1:1:length(bits)
    if bits(j)==1
        x=repmat(bits(j),1,T_b);
    else
        x=repmat(bits(j),1,T_b);
    end
   BITS =cat(2,BITS,x);
end

subplot(2,1,1);
plot(UNRZ);
title('INPUT BITS IN UNRZ');
grid on;

subplot(2,1,2);
plot(BITS);
title('RECOVERED BITS IN UNRZ');
grid on;


check=biterr(bit_pattern,bits);
if check==0
    disp('no bit error')
else
    disp('bit error')
end