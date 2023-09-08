%NAME:ARCHIT VASHIST
%ROLL NO:2021BEC0022


%% CDMA IMPLEMENTATION IN MATLAB
clc;
clear all;
close all;


%bit generation
%b_k=round(rand(1,5));
b_k=[0 1 1 0 1];
disp('Input signal bits');
disp(b_k);

% career frequency
fc=30;

%sampling frequency
fs=20*fc;

%sampling time
Ts=1/fs;

%time interval for sampling
t=0:Ts:1-Ts;

%bit length
T_b=length(t)/length(b_k);
T_c=T_b/4;

%career wave
phi_t=sqrt(2/T_b)*cos(2*pi*fc*t);

%select user
user_selection=menu('SELECT USER' ,'USER1','USER2','USER3','USER4');

%%orthogonal codes allocated to each user
user1=[1 1 0 0 ];
user2=[1 0 0 1 ];
user3=[1 1 1 1 ];
user4=[1 0 1 0 ];

%orthogonal codes allocated to the user
if user_selection==1
        user=user1
else if user_selection==2;
        user=user2
else if user_selection==3;
        user=user3
else if user_selection==4;
        user=user4
end
end
end
end

PNRZ_b_k=[]; %PNRZ storing variable of message signal

%POLAR NRZ creation
for j=1:1:length(b_k)
    if b_k(j)==1
        x_t=repmat(b_k(j),1,T_b);   
    else
        x_t=repmat(-1,1,T_b);
    end
    PNRZ_b_k=cat(2,PNRZ_b_k,x_t);
end

PNRZ_c_k=[]; %PNRZ storing variable of carrier signal

%POLAR NRZ creation
for i=1:length(b_k)
for j=1:1:length(user)
    if user(j)==1
        x_t=repmat(user(j),1,T_c);   
    else
        x_t=repmat(-1,1,T_c);
    end
    PNRZ_c_k=cat(2,PNRZ_c_k,x_t);
end
end

figure(1)
%PLOT MESSAGE SIGNAL
subplot(5,1,1);
plot(PNRZ_b_k);
title('MESSAGE SIGNAL b_k(t)')
grid on;


%plot career signal
subplot(5,1,2);
plot(PNRZ_c_k);
title('ORTHOGONAL SIGNAL c_k(t)');
grid on;


%MODULATED SIGNAL
m_t=PNRZ_c_k.*PNRZ_b_k;


%plot modulated signal
subplot(5,1,3);
plot(m_t);
title('PRODUCT SIGNAL m(t)');
grid on;

x_t=phi_t.*m_t;


subplot(5,1,4);
plot(phi_t);
title('CARRIER SIGNAL c_t');
grid on;

subplot(5,1,5);
plot(x_t);
title('MODULATED SIGNAL x(t)');
grid on;



%%DEMODULATION SECTION

y_t=x_t.*phi_t; %STEP1->product modulator

figure(2)
subplot(5,1,1);
plot(y_t);
title('PRODUCT MODULATOR (x_t*phi_t)');
grid on;


%%INTEGRATOR
bit_s = []; %integrator values


%TAKE t_c INTERVAL AND DO INTEGRATION FOR DECISION RULE 
for j = 1:length(user)*length(b_k)
    % calculate the start and end index for integration interval
    start_index = (j-1)*T_c+1;
    end_index = j*T_c;

    % integrator
    bit_sum = sum(y_t(start_index:end_index));
    
    % storing integrator values in an array for message signal creation
    bit_s = [bit_s bit_sum];

end



%DECISION RULE
bits=[];
for j=1:length(user1)*length(b_k)
    if bit_s(j)>0 %THRESHOLD VALUE
        a=1;
    else
        a=-1;
    end
   bits =[bits a];
end

PNRZ_bits=[]; 
%PNRZ storing variable of bits after decision rule

%POLAR NRZ creation
for j=1:length(bits)
    if bits(j)==1
        x_t=repmat(bits(j),1,T_c);   
    else
        x_t=repmat(-1,1,T_c);
    end
    PNRZ_bits=cat(2,PNRZ_bits,x_t);
end

subplot(5,1,2);
plot(PNRZ_c_k);
title('PNRZ of Orthogonal signal');
grid on;

subplot(5,1,3);
plot(PNRZ_bits);
title('PNRZ of recovered signal after first product modulator');
grid on;

%%MULTIPLYING WITH ORTHOGONAL SIGNAL AFTER DECISION MAKING
u_t=PNRZ_c_k.*PNRZ_bits;

subplot(5,1,4);
plot(u_t);
title('recovered signal after decision rule u_t');
grid on;

subplot(5,1,5);
plot(PNRZ_b_k);
title('ACTUAL b_k signal sent');
grid on;

bit_s1=[];
%TAKE t_b INTERVAL AND DO INTEGRATION FOR DECISION RULE 
for j = 1:length(b_k)
    % calculate the start and end index for integration interval
    start_index = (j-1)*T_b+1;
    end_index = j*T_b;

    % integrator
    bit_sum = sum(u_t(start_index:end_index));
    
    % storing integrator values in an array for message signal creation
    bit_s1 = [bit_s1 bit_sum];

end



%DECISION RULE
bits1=[];
for j=1:length(b_k)
    if bit_s1(j)>0 %THRESHOLD VALUE
        a=1;
    else
        a=0;
    end
   bits1 =[bits1 a];
end
if(biterr(b_k,bits1))
    disp('BIT ERROR ')
else
    disp('NO BIT ERROR')
end



