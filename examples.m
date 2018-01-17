% This MATLAB program checks the feasibility of LMIs from Theorems 1, 2, 3 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data implementation of
% derivative-dependent control using artificial delays," IEEE Transactions on Automatic Control, 2018. 
%% Section II.C: chain of three integrators 
A=diag(ones(1,2),1);        % |
B=[0; 0; 1];                % | Parameters from (15)
C=[1 0 0];                  % |
Kbar={-2e-4, -.06, -.342};  % Gains of (4)
q=30*(1:2);                 % Discrete-time delays
alpha=1e-3;                 % Decay rate 

% Delayed sampled-data controller (6)
h=.044;                     % Sampling 
[K0,K]=LMI_TAC18_th1(A,B,C,Kbar,h,q,alpha); 
if isempty(K0)
    disp('Theorem 1: not feasible'); 
else
    disp('Theorem 1: feasible'); 
    disp(['K0=' mat2str(K0)]); 
    for i=1:length(K)
        disp(['K' num2str(i) '=' mat2str(K{i})]); 
    end
end
disp('-----------------------'); 

% Event-triggered controller (6), (13), (14)
sigma=2e-3;     % Event-triggering parameter 
h=.042;         % Sampling

[K0,K,Omega]=LMI_TAC18_th2(A,B,C,Kbar,h,q,alpha,sigma); 
if isempty(K0)
    disp('Theorem 2: not feasible'); 
else
    disp('Theorem 2: feasible'); 
    disp(['K0=' mat2str(K0)]); 
    for i=1:length(K)
        disp(['K' num2str(i) '=' mat2str(K{i})]); 
    end
    disp(['Omega=' mat2str(Omega)]); 
end
disp('-----------------------'); 

%% Section III.B: PID control 
a1=8.4; a2=0; b=35.71;              % Parameters of the system (16) 
kpbar=-10; kibar=-40; kdbar=-.65;   % Gains of continuous PID 
q=7;                                % Discrete-time delays
alpha=5;                            % Decay rate 
tf=10;                              % Time of simulations 

% Delayed sampled-data controller (18)
h=4.7e-3;     % Sampling 
[kp,ki,kd]=LMI_TAC18_th3(a1,a2,b,kpbar,kibar,kdbar,h,q,alpha,0); 
if isempty(kp)
    disp('Theorem 3 (sigma=0): not feasible'); 
else
    disp('Theorem 3 (sigma=0): feasible'); 
    disp(['kp=' mat2str(kp)]); 
    disp(['ki=' mat2str(ki)]); 
    disp(['kd=' mat2str(kd)]); 
end
disp('-----------------------'); 

% Event-triggered controller (18), (21), (22)
sigma=9e-3;     % Event-triggering parameter 
h=4e-3;         % Sampling
[kp,ki,kd]=LMI_TAC18_th3(a1,a2,b,kpbar,kibar,kdbar,h,q,alpha,sigma); 
if isempty(kp)
    disp(['Theorem 3 (sigma=' num2str(sigma) '): not feasible']); 
else
    disp(['Theorem 3 (sigma=' num2str(sigma) '): feasible']); 
    disp(['kp=' mat2str(kp)]); 
    disp(['ki=' mat2str(ki)]); 
    disp(['kd=' mat2str(kd)]); 
end
disp('-----------------------'); 