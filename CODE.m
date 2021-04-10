%% ELEC 4700 - Assignment 4
% Saifuddin Mohammed, #101092093
set(0,'DefaultFigureWindowStyle','docked')
set(0,'defaultaxesfontsize',20)
set(0,'defaultaxesfontname','Times New Roman')
set(0,'DefaultLineLineWidth',2);
close all;
clear all; 
clc;

%% A linear fit was performed in order to obtain the value of R3. The value of R3 obtained was 25 Ohms.

% Setting up Parameters
% The MNA PA was utilized for to perform the remainder of the assignment. 

% Parameters
R_1 = 1;
C_1 = 0.25;
R_2 = 2;
L_1 = 0.2;
aL = 100;
R_3 = 25;
R_4 = 0.1;
R_0 = 1000;
V_INP=1;

% The first step will be to set up the G and C matrices as per the KCL
% equations.

%Initializing a 7X7 matrix for 7 equations
G=zeros(7,7);
% Setting up the G matrix.
G(1,:)=[1 0 0 0 0 0 0];
G(2,:)=[(-1/R_1) (1/R_2+1/R_1) 0 0 0 1 0];
G(3,:)=[0 0 1/R_3 0 0 -1 0]; 
G(4,:)=[0 0 -1*aL/R_3 1 0 0 -1]; 
G(5,:)=[0 0 0 -1/R_4 (1/R_4+1/R_0) 0 0];
G(6,:)=[0 -1 1 0 0 0 0]; 
G(7,:)=[0 0 0 1 0 0 -aL];

% G = [ 1 0 0 0 0 0 0; 
%     (-1/R_1) (1/R_2+1/R_1) 0 0 0 1 0;
%     0 0 1/R_3 0 0 -1 0; 
%     0 0 -1*aL/R_3 1 0 0 -1;
%     0 0 0 -1/R_4 (1/R_4+1/R_0) 0 0;
%     0 -1 1 0 0 0 0;
%     0 0 0 1 0 0 -aL]; 




% Initializing a 7x7 matrix for the Conductances
C=zeros(7,7);
% Filling in the values as per the equations  
C(1,:)=[0 0 0 0 0 0 0]; 
C(2,:)=[-C_1 C_1 0 0 0 0 0];
C(3,:)=[0 0 0 0 0 0 0]; 
C(4,:)=[0 0 0 0 0 0 0 ]; 
C(5,:)=[0 0 0 0 0 0 0];
C(6,:)=[0 0 0 0 0 L_1 0 ]; 
C(7,:)=[0 0 0 0 0 0 0]; 

% To get G and C matrices. 
C
G



% C = [0 0 0 0 0 0 0;
%     -C_1 C_1 0 0 0 0 0;
%     0 0 C_n 0 0 0 0;
%     0 0 0 0 0 0 0;
%     0 0 0 0 0 0 0;
%     0 0 0 0 0 L_1 0;
%     0 0 0 0 0 0 0];
    
% Q3 -  Part b i)

% Objective is to perform a  DC sweep. 

Ts=20;

DC= zeros(3,Ts);

% The DC sweep will be performed from -10V to 10V
DC(1,:)=linspace(-10,10,Ts);

for sweep=1:Ts
    
   % Solving using GV = F
   
    F=[DC(1,sweep); 0; 0; 0; 0; 0; 0];   % Setting up the F array.
    
    V=G\F;
    
    DC(2,sweep)=V(5);
    
    DC(3,sweep)=V(3);
    
end


% Plotting the DC sweep of V1 and V_0 and V_3 along with it. 
figure(1)
hold on;
plot(DC(1,:),DC(2,:), 'r');
plot(DC(1,:),DC(3,:), 'b');
hold off;
title('DC Sweep Plot,  SM101092039');
ylabel('Voltage [V]');
xlabel('Input Voltage [V]');
legend('V_0','V_3');



% Q3-  Part b ii)
% The objective of this part was to perform AC analysis of the circuit.
% Similar proecdure was followed as DC sweep, however omega was taken into
% consideration as it is an AC analysis.

Ts=2000;
Data2=zeros(2,Ts);
Data2(1,:)=linspace(0,500,Ts);


for sweep=1:Ts
    
    F=[1; 0; 0; 0; 0; 0; 0];  % Setting up the F Matrix
    W=Data2(1,sweep); % Getting the value for Omega,w
    
   % Solving using GV=F, but the frequency domain is considered.
   
    V=(G+1j*W*C)\F;
    
    Data2(2,sweep)=V(5);
   
end


% Plot the output voltage V_0 as the function of Angular frequency
figure(2)
plot(Data2(1,:),real(Data2(2,:)), 'r');
title(' Plot of V_0 as a Function of \omega, SM101092039');
xlabel('Angular frequency, \omega [rad/s]');
ylabel('V_O (V)');


% Plot the output voltage V_0 as the function of Angular frequency
% To get the gain, Just divide V_0/V_I and conver it to dB.

figure(3)
A_v = 20*log10(real(Data2(2,:))./V_INP);
semilogx(Data2(1,:),A_v, 'r');
title('Plot Gain as a Function of \omega, SM101092039');
xlabel('\omega [rad/s]');
ylabel('Gain [dB]')



% Q3- Part b iii)
% The objective of this part to generate a histogram of gain perturbations
% utilizng the gain obtained in the previous part

W=pi;

samples=10000;

Perturb=zeros(1,samples);

for sweep=1:samples
    
    C_new=C_1+randn()*0.05; % Congifuration of random capacitances with standard deviation of 0.05
    
    C(2,:)=[-C_new +C_new 0 0 0 0 0];
    
    % Setting up the F matrix
    F=[1; 0; 0; 0; 0; 0; 0];
    
    % GV = F is being uztilized to solve
    V=(G+1j*W*C)\F;
    
    Perturb(1,sweep)=V(5);
    
end

% Generating a historgram depicting the gain perturbations.
figure(4);
histogram(real(Perturb),30, 'FaceColor', 'g');
title('Gain for perturbations, SM101092039');
xlabel('Gain, A_V');
ylabel('Total Samples')



%% Q4 Part d i)
% In this part 3 different input configurations will be tested.


C(2,:)=[-C_1 C_1 0 0 0 0 0];

Sim_Set=1;
StepS=1000;

%Setting up the step size
St_Size=Sim_Set/StepS;


Freq_Sim1=1/St_Size;

Sim_param=StepS+1;

F_sim =(-(Sim_param-1)/2:(Sim_param-1)/2)*(Freq_Sim1/Sim_param);


Range_F=zeros(3,3,StepS+1);

V_X=[0; 0; 0; 0; 0; 0; 0];

Sig_Gen=zeros(3,StepS+1);

Ratio=1/0.03;

% Initating a For Loop for impelementing the different input signals.

for In_Sig=1:3
    
    for sweep=1:StepS
        
        % For step input that transitions from 0 to 1 at 0.03s
        
        if(In_Sig==1)
            
            if(sweep*St_Size<0.03)
                
                Sig_Gen(In_Sig,sweep)=0;
                
            else
                
                Sig_Gen(In_Sig,sweep)=1;
                
            end
            
            
        % For a sine function input    
        elseif(In_Sig==2)
            
            Sig_Gen(In_Sig,sweep)=sin(2*pi*Ratio*sweep*St_Size);
       
            
        % For the Gaussian Input 
        
        else
            
            
            Sig_Gen(In_Sig,sweep)=exp(-(sweep*St_Size-0.1)^2/(2*0.03^2));
        end
        
        
        
    end
    
    
end


for In_Sig=1:3
    
    V_X=[0; 0; 0; 0; 0; 0; 0];
    
    
    for sweep=1:StepS
        
        
        V_INP=Sig_Gen(In_Sig,sweep);
        F=[V_INP; 0; 0; 0; 0; 0; 0];
        
        
        Range_F(In_Sig,1,sweep+1)=sweep*St_Size;
        
        Range_F(In_Sig,2,sweep+1)=V_INP;
        
        
        
        % For transient Analysis, the eqaution V.dV/dt + GV = F must be
        % solved.
        % Solving the equation gives
        A=C/St_Size+G;
        
        V=(A)\(C*V_X/St_Size+F);
        
        Range_F(In_Sig,3,sweep+1)=V(5);
        

        V_X=V;
        
    end
    
    
end


% Finally plot all the responses by checking the type of the input in the
% time and the frequency domains.

% Generating a plot for the step input

Response=1;
In_1=zeros(3,StepS+1);
In_1(:,:)=Range_F(Response,:,:);

figure(5)
subplot(2,1,1)
hold on;
plot(In_1(1,:),In_1(2,:), 'r');
plot(In_1(1,:),In_1(3,:), 'b');
hold off;
legend('V_I','V_0');
title('Time Domain Response for Step Input, SM101092039');
xlabel('t [s]');
ylabel('Voltage [V]');

% Getting the Frequency Domain Response
% The frequency domain response can be obtained through FFT. 

subplot(2,1,2)
% Perform FFT for conversion from time to frequency
F1=fft(In_1(2,:));
F2=fft(In_1(3,:));
hold on;
plot(F_sim,fftshift(abs(F1)), 'g');
plot(F_sim,fftshift(abs(F2)), 'b');
hold off;
legend('V_I','V_0');
title('Frequency Domain Response from the Step Input, SM101092039');
xlabel('F[Hz]');
ylabel('Mag');




% Generating a plot for the sine function input
Response=2;
In_2=zeros(3,StepS+1);
In_2(:,:)=Range_F(Response,:,:);

figure(6)
subplot(2,1,1)
hold on;
plot(In_2(1,:),In_2(2,:), 'r');
plot(In_2(1,:),In_2(3,:),'b');
hold off;
legend('V_I','V_0');
title('Response from Sinusoidal Input, SM101092039');
xlabel('t [s]');
ylabel('Voltage [V]');

% Getting the Frequency Domain Response
% The frequency domain response can be obtained through FFT. 

subplot(2,1,2)
% Perform FFT for conversion from time to frequency
F1=fft(In_2(2,:));
F2=fft(In_2(3,:));
hold on;
plot(F_sim,fftshift(abs(F1)), 'g');
plot(F_sim,fftshift(abs(F2)), 'b');
hold off;
legend('V_i','V_O');
title('Frequency Domain Response for Sine Function Input, SM101092039');
xlabel('F[Hz]');
ylabel('Mag');


% Generating a plot for the gaussian pulse input
Response=3;
In_3=zeros(3,StepS+1);
In_3(:,:)=Range_F(Response,:,:);

figure(7)
subplot(2,1,1)
hold on;
plot(In_3(1,:),In_3(2,:), 'r');
plot(In_3(1,:),In_3(3,:),'b');
hold off;
legend('V_I','V_0');
title('Response from Gaussian Pulse Input, SM101092039');
xlabel('t [s]');
ylabel('Voltage [V]');

% Getting the Frequency Domain Response
% The frequency domain response can be obtained through FFT. 

subplot(2,1,2)
% Perform FFT for conversion from time to frequency
F1=fft(In_3(2,:));
F2=fft(In_3(3,:));
hold on;
plot(F_sim,fftshift(abs(F1)), 'g');
plot(F_sim,fftshift(abs(F2)), 'b');
hold off;
legend('V_I','V_0');
title('Frequency Domain Response for the Gaussian Pulse Input, SM101092039');
xlabel('F[Hz]');
ylabel('Mag');




%% Q5 Cicruit with Noise
% a) The circuit was updated to have a current source in parallel with R3
% b) A new capacitor C_n was added to the circuit and the C amtrix was
% expected to be modified.

C_n = 0.00001;
I_n = 0.001; 


% Q5 Part c i) 
% The implemenration procedure is similar to the previous part where noise
% was not added. The only addition is that noise was modelled and will be considered in simulations
% now.

C(2,:)=[-C_1 C_1 0 0 0 0 0];

StepS=1000;

Sim_Set=1;

St_Size=Sim_Set/StepS;

Freq_Sim1=1/St_Size;

Sim_param=StepS+1;

F_sim=(-(Sim_param-1)/2:(Sim_param-1)/2)*(Freq_Sim1/Sim_param);



Range_F=zeros(3,3,StepS+1);

V_X=[1; 0; 0; 0; 0; 0; 0];


Sig_Gen=zeros(3,StepS+1);

Ratio=1/0.03;


% Initating a For Loop for impelementing different values of the capacitor

for In_Sig=1:3
    
    for iter=1:StepS
        
        % A gaussian pulse with noise
        Sig_Gen(In_Sig,iter)=exp(-(iter*St_Size-0.1)^2/(2*0.03^2));
        
    end
    
end


for In_Sig=1:3

    if(In_Sig==1)
        
        Cn_mod=C_n;
        
    elseif(In_Sig==2)
        
        
        Cn_mod=C_n*1000;
        
        
    else
        
        
        Cn_mod=C_n*10000;
    end
    
    C(3,:)=[0 0 Cn_mod 0 0 0 0]; 
    
    V_X=[0; 0; 0; 0; 0; 0; 0];

    
    
    for iter=1:StepS
        
        V_INP=Sig_Gen(In_Sig,iter);

        % The new current source I_n
        I_n=0.001*randn();
        
        F=[V_INP; 0; I_n; 0; 0; 0; 0];
        Range_F(In_Sig,1,iter+1)=iter*St_Size;
        Range_F(In_Sig,2,iter+1)=V_INP;
        
        A=C/St_Size+G;
        V=(A)\(C.*V_X/St_Size+F);
        Range_F(In_Sig,3,iter+1)=V(5);

        V_X=V;
    end
    
end


% Updating the C matrix to include the new capacitor
C(3,:)=[0 0 C_n 0 0 0 0]; 



% Plotting the response of the Gaussian pulse with noise
Response=1;
In_1=zeros(3,StepS+1);
In_1(:,:)=Range_F(Response,:,:);
figure(8)
subplot(2,1,1)
hold on;
plot(In_1(1,:),In_1(2,:), 'g');
plot(In_1(1,:),In_1(3,:), 'b');
hold off;
legend('V_I','V_0');
title(' Time Domain Response of Gaussian Input pulse with Noise, SM101092039');
xlabel('t [s]');
ylabel('Voltage [V]');

%Plotting the Freqeuncy Response
subplot(2,1,2)
F1=fft(In_1(2,:));
F2=fft(In_1(3,:));
hold on;
plot(F_sim,fftshift(abs(F1)), 'r');
plot(F_sim,fftshift(abs(F2)), 'b');
hold off;
legend('V_I','V_0');
title(' Frequency Domain Response for Gaussian Input pulse with Noise, SM101092039');
xlabel('F [Hz]');
ylabel('Mag');



% Testing the response of the Cicruit with different values of Cout
Response=1;
In_1=zeros(3,StepS+1);
In_1(:,:)=Range_F(Response,:,:);
figure(9)
hold on;
plot(In_1(1,:),In_1(2,:), 'r');
plot(In_1(1,:),In_1(3,:), 'b');
hold off;
legend('V_I','V_O');
title('Response for Cout = 0.00001, SM101092039')
xlabel('t [s]');
ylabel('Voltage [V]');



Response=2;
In_2=zeros(3,StepS+1);
In_2(:,:)=Range_F(Response,:,:);
figure(10)
hold on;
plot(In_2(1,:),In_2(2,:), 'r');
plot(In_2(1,:),In_2(3,:), 'b');
hold off;
legend('V_i','V_O');
title(' Time Domain Response for Cout= 10 mF, SM101092039');
xlabel('t [s]');
ylabel('Voltage [V]');



Response=3;
In_3=zeros(3,StepS+1);
In_3(:,:)=Range_F(Response,:,:);
figure(11)
hold on;
plot(In_3(1,:),In_3(2,:), 'r');
plot(In_3(1,:),In_3(3,:), 'b');
hold off;
legend('V_i','V_O');
title('Time Domain Response for Cout = 100 mF, SM101092039 ');
xlabel('t [s]');
ylabel('Voltage [V]');




% Testing the Response of the circuit with noise but different time step

C(2,:)=[-C_1 +C_1 0 0 0 0 0];

New_Ts=10000;

Sim_Set=1;

New_StS=Sim_Set/New_Ts;

F_Sim1=1/New_StS;

Sim_Final=New_Ts+1;

Fsim_Final=(-(Sim_Final-1)/2:(Sim_Final-1)/2)*(F_Sim1/Sim_Final);

New_Sim=zeros(3,3,New_Ts+1);

V_X=[1; 0; 0; 0; 0; 0];

Sig_Gen=zeros(3,New_Ts+1);

Ratio=1/0.03;

for In_Sig=1:1
    
    for iter=1:New_Ts
        
        Sig_Gen(In_Sig,iter)=exp(-(iter*New_StS-0.1)^2/(2*0.03^2));
        
    end
    
end


for In_Sig=1:1
    
    
    C(3,:)=[0 0 C_n 0 0 0 0];            % Similar implemenation process as the previous parts

    V_X=[0; 0; 0; 0; 0; 0; 0];

    for iter=1:New_Ts
        
        V_INP=Sig_Gen(In_Sig,iter);

        I_n=0.001*randn();
        
        F=[V_INP; 0; -I_n; 0; 0; 0; 0];
        
        New_Sim(In_Sig,1,iter+1)=iter*New_StS;
        
        New_Sim(In_Sig,2,iter+1)=V_INP;
        
        A=C/New_StS+G;
        
        V=(A)\(C*V_X/New_StS+F);
        
        New_Sim(In_Sig,3,iter+1)=V(5);

        V_X=V;
    end
    
end

C(3,:)=[0 0 C_n 0 0 0 0]; 

C


% Plotting the responses with different Time steps

% Response for 1000 timesteps
figure(12)
Response=1;
In_1=zeros(3,StepS+1);
In_1(:,:)=Range_F(Response,:,:);
hold on;
plot(In_1(1,:),In_1(2,:), 'r');
plot(In_1(1,:),In_1(3,:), 'b');
hold off;
legend('V_I','V_0');
title('Time Domain Response for 1000 timesteps, SM101092039');
xlabel('t [s]');
ylabel('Voltage [V]');


% Response for 10,000 timesteps
figure(13)
Response=1;
Fsim1=zeros(3,New_Ts+1);
Fsim1(:,:)=New_Sim(Response,:,:);
hold on;
plot(Fsim1(1,:),Fsim1(2,:), 'r');
plot(Fsim1(1,:),Fsim1(3,:), 'b');
hold off;
legend('V_I','V_0');
title('Time Domain Response for 10000 timesteps, SM101092039');
xlabel('t [s]');
ylabel('Voltage [V]');

