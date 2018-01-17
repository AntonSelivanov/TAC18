function [kp,ki,kd]=LMI_TAC18_th3(a1,a2,b,kpbar,kibar,kdbar,h,q,alpha,sigma)
% This MATLAB program checks the feasibility of LMIs from Theorem 3 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data implementation of derivative-dependent control using artificial delays," IEEE Transactions on Automatic Control, 2018.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% a1, a2, b             - the parameters of the system (16) 
% kpbar, kibar, kdbar   - controller gains from (17)
% h                     - sampling
% q                     - discrete delay
% sigma                 - event-triggering threshold from (22) 
% alpha                 - decay rate 

% Output: 
% kp, ki, kd            - controller gains, are empty if not feasible

%% Calculation of the controller gains using (19)
kp=kpbar+kdbar/(q*h); 
ki=kibar; 
kd=-kdbar/(q*h); 
%% Decision variables 
P=sdpvar(3); 
S=sdpvar(3); % Psi<=0 implies S, W, R, omega >= 0
sdpvar W R omega 
%% Notations (25) and G
A=[0 1 0; -a2+b*(kp+kd), -a1-q*h*b*kd, b*ki; 1 0 0]; 
Av=[0 0 0; b*kp 0 b*ki; 1 0 0]; 
B=[0; b; 0]; 
G=h^2*exp(2*alpha*h)*S+[0 0 0; 0 1 0; 0 0 0]*R*kd^2*(q*h)^2; 
%% The LMIs
Psi=blkvar; 
Psi(1,1)=P*A+A'*P+2*alpha*P+[0 0 0; 0 1 0; 0 0 0]*W*kd^2*h^2*exp(2*alpha*h); 
Psi(1,2)=P*Av*sqrt(h);
Psi(1,3)=P*B; 
Psi(1,4)=P*B; 
Psi(1,5)=P*B; 
Psi(1,6)=[kp+kd; -q*h*kd; ki]*omega*sigma; 
Psi(1,7)=A'*G; 
Psi(2,2)=-pi^2/4*S*h; 
Psi(2,6)=[kp; 0; ki]*omega*sigma*sqrt(h); 
Psi(2,7)=Av'*G*sqrt(h); 
Psi(3,3)=-W*pi^2/4*exp(-2*alpha*q*h); 
Psi(3,6)=omega*sigma; 
Psi(3,7)=B'*G; 
Psi(4,4)=-R*4/(q*h)^2*exp(-2*alpha*q*h); 
Psi(4,6)=omega*sigma; 
Psi(4,7)=B'*G; 
Psi(5,5)=-omega; 
Psi(5,7)=B'*G; 
Psi(6,6)=-omega*sigma; 
Psi(7,7)=-G; 
Psi=sdpvar(Psi); 

LMIs=[P>=0,Psi<=0];

%% Solution of LMIs
options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

if sol.problem == 0
    primal=check(LMIs); 
    if min(primal)<=0
        kp=[]; ki=[]; kd=[];  
    end
else
    kp=[]; ki=[]; kd=[];  
    yalmiperror(sol.problem) 
end