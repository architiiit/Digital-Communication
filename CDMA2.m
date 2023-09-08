%NAME:ARCHIT VASHIST
%ROLL NO:2021BEC0022


%% CDMA IMPLEMENTATION IN MATLAB
clc;
clear all;
close all;

tic
%bit generation
b_k1=round(rand(1,5));
b_k2=round(rand(1,5));
disp('Input signal bits 1');
disp(b_k1);
disp('Input signal bits 2');
disp(b_k2);


% career frequency
fc=100;

%sampling frequency
fs=20*fc;

%sampling time
Ts=1/fs;

%time interval for sampling
t=0:Ts:1-Ts;

%bit length
T_b=length(t)/length(b_k1);
T_c=T_b/4;

%career wave
phi_t=sqrt(2/T_b)*cos(2*pi*fc*t);

%%orthogonal codes allocated to each user
user1=[1 1 0 0 ];
user2=[1 0 0 1 ];



PNRZ_b_k1=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(b_k1)
    if b_k1(j)==1
        x_t=repmat(b_k1(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_b_k1=cat(2,PNRZ_b_k1,x_t);
end


PNRZ_b_k2=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(b_k2)
    if b_k2(j)==1
        x_t=repmat(b_k2(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_b_k2=cat(2,PNRZ_b_k2,x_t);
end


PNRZ_c_k1=[]; %PNRZ storing variable of carrier signal

%POLAR NRZ creation
for i=1:length(b_k1)
for j=1:1:length(user1)
    if user1(j)==1
        x_t=repmat(user1(j),1,T_c);   
    else
        x_t=repmat(-1,1,T_c);
    end
    PNRZ_c_k1=cat(2,PNRZ_c_k1,x_t);
end
end

PNRZ_c_k2=[]; %PNRZ storing variable of carrier signal
%POLAR NRZ creation
for i=1:length(b_k2)
for j=1:1:length(user2)
    if user2(j)==1
        x_t=repmat(user2(j),1,T_c);   
    else
        x_t=repmat(-1,1,T_c);
    end
    PNRZ_c_k2=cat(2,PNRZ_c_k2,x_t);
end
end



figure(1)
%PLOT MESSAGE SIGNAL
subplot(9,1,1);
plot(PNRZ_b_k1,'MarkerSize',3, ...
    'LineWidth',3);
title('MESSAGE SIGNAL b_k1')
grid on;


%plot career signal
subplot(9,1,2);
plot(PNRZ_b_k2,'MarkerSize',3, ...
    'LineWidth',3);
title('MESSAGE SIGNAL b_k2');
grid on;
%plot career signal
subplot(9,1,3);
plot(PNRZ_c_k1,'MarkerSize',3, ...
    'LineWidth',3);
title('ORTHOGONAL SIGNAL c_k1');
grid on;

%MODULATED SIGNAL
m_t1=PNRZ_c_k1.*PNRZ_b_k1;
m_t2=PNRZ_c_k2.*PNRZ_b_k2;

%plot modulated signal
subplot(9,1,4);
plot(PNRZ_c_k2,'MarkerSize',3, ...
    'LineWidth',3);
title('ORTHOGONAL SIGNAL c_k2');
grid on;

%plot modulated signal
subplot(9,1,5);
plot(m_t1,'MarkerSize',3, ...
    'LineWidth',3);
title('PRODUCT SIGNAL m1(t)');
grid on;

%plot modulated signal
subplot(9,1,6);
plot(m_t2,'MarkerSize',3, ...
    'LineWidth',3);
title('PRODUCT SIGNAL m2(t)');
grid on;





s_t1=m_t1.*phi_t;
s_t2=m_t2.*phi_t;

subplot(9,1,7);
plot(s_t1,'MarkerSize',1, ...
    'LineWidth',1);
title('MODULATED SIGNAL s_t1');
grid on;

subplot(9,1,8);
plot(s_t2,'MarkerSize',1, ...
    'LineWidth',1);
title('MODULATED SIGNAL s_t2');
grid on;

s_t=s_t1+s_t2;

%plot modulated signal
subplot(9,1,9);
plot(s_t,'MarkerSize',1, ...
    'LineWidth',1);
title('MODULATED SIGNAL S(t)');
grid on;

%%DEMODULATION SECTION

r_t=s_t.*phi_t;

figure(2)
subplot(8,1,1)
plot(r_t,'MarkerSize',1, ...
    'LineWidth',1);
title('recovered signal (product modulator)')
grid on;



bit_s = []; %integrator values
%TAKE t_c INTERVAL AND DO INTEGRATION FOR DECISION RULE 
for j = 1:length(user1)*length(b_k1)
    % calculate the start and end index for integration interval
    start_index = (j-1)*T_c+1;
    end_index = j*T_c;

    % integrator
    bit_sum = sum(r_t(start_index:end_index));
    
    % storing integrator values in an array for message signal creation
    bit_s = [bit_s bit_sum];

end



%DECISION RULE
bits=[];
for j=1:length(user1)*length(b_k1)
    if bit_s(j)>0 %THRESHOLD VALUE
        a=1;
    elseif bit_s(j)==0
        a=0;
    else
        a=-1;
    end
   bits =[bits a];
end

PNRZ_bits=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(bits)
    if bits(j)==1
        x_t=repmat(bits(j),1,T_c);   
    elseif bits(j)==-1
        x_t=repmat(-1,1,T_c);
    else
        x_t=zeros(1,T_c);
    end
    PNRZ_bits=cat(2,PNRZ_bits,x_t);
end



subplot(8,1,2)
plot(PNRZ_bits,'MarkerSize',3, ...
    'LineWidth',3);
title('recovered signal BITS')
grid on;


u_1=PNRZ_bits.*PNRZ_c_k1;
u_2=PNRZ_bits.*PNRZ_c_k2;


subplot(8,1,3)
plot(u_1,'MarkerSize',1, ...
    'LineWidth',3);
title('PNRZ USER1(PRODUCT MODULATOR WITH PN SEQUENCE OF USER1)')
grid on;

subplot(8,1,4)
plot(u_2,'MarkerSize',1, ...
    'LineWidth',3);
title('PNRZ USER2(PRODUCT MODULATOR WITH PN SEQUENCE OF USER2)')
grid on;

bit_s1=[];
%TAKE t_b INTERVAL AND DO INTEGRATION FOR DECISION RULE 
for j = 1:length(b_k1)
    % calculate the start and end index for integration interval
    start_index = (j-1)*T_b+1;
    end_index = j*T_b;

    % integrator
    bit_sum = sum(u_1(start_index:end_index));
    
    % storing integrator values in an array for message signal creation
    bit_s1 = [bit_s1 bit_sum];

end



%DECISION RULE
bits1=[];
for j=1:length(bit_s1)
    if bit_s1(j)>0 %THRESHOLD VALUE
        a=1;
    else
        a=0;
    end
   bits1 =[bits1 a];
end




bit_s2=[];
%TAKE t_b INTERVAL AND DO INTEGRATION FOR DECISION RULE 
for j = 1:length(b_k2)
    % calculate the start and end index for integration interval
    start_index = (j-1)*T_b+1;
    end_index = j*T_b;

    % integrator
    bit_sum = sum(u_2(start_index:end_index));
    
    % storing integrator values in an array for message signal creation
    bit_s2 = [bit_s2 bit_sum];

end



%DECISION RULE
bits2=[];
for j=1:length(bit_s2)
    if bit_s2(j)>0 %THRESHOLD VALUE
        a=1;
    else
        a=0;
    end
   bits2 =[bits2 a];
end


PNRZ_bits1=[];
%POLAR NRZ creation
for j=1:1:length(bits1)
    if bits1(j)==1
        x_t=repmat(bits1(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_bits1=cat(2,PNRZ_bits1,x_t);
end

PNRZ_bits2=[];
%POLAR NRZ creation
for j=1:1:length(bits2)
    if bits2(j)==1
        x_t=repmat(bits2(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_bits2=cat(2,PNRZ_bits2,x_t);
end



subplot(8,1,5)
plot(PNRZ_bits1,'MarkerSize',1, ...
    'LineWidth',3);
title('PNRZ USER1(RECOVERED SIGNAL)')
grid on;

subplot(8,1,6)
plot(PNRZ_b_k1,'MarkerSize',1, ...
    'LineWidth',3);
title('PNRZ USER1(ACTUAL SENT SIGNAL)')
grid on;


subplot(8,1,7)
plot(PNRZ_bits2,'MarkerSize',1, ...
    'LineWidth',3);
title('PNRZ USER2(RECOVERED SIGNAL)')
grid on;

subplot(8,1,8)
plot(PNRZ_b_k2,'MarkerSize',1, ...
    'LineWidth',3);
title('PNRZ USER2(ACTUAL SENT SIGNAL)')
grid on;

if biterr(b_k1,bits1)==0
    disp('NO BIT ERROR FOR USER1')
else 
    disp('BIT ERROR FOR USER1')
end

if biterr(b_k2,bits2)==0
    disp('NO BIT ERROR FOR USER2')
else 
    disp('BIT ERROR FOR USER2')
end

toc


























