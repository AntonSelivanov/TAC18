function [K0,K]=LMI_TAC18_th1(A,B,C,Kbar,h,q,alpha)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov and E. Fridman, "Sampled-data implementation of derivative-dependent control using artificial delays," IEEE Transactions on Automatic Control, 2018.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)

% Input: 
% A, B, C  - the parameters of the LTI system (1) 
% Kbar     - cell array of \bar{K}_0,...,\bar{K}_{r-1} from (4)
% h        - sampling
% q        - vector of discrete delays
% alpha    - decay rate 

% Output: 
% K0, K - controller gains, K0 is empty if not feasible

%% Dimensions
[n,m]=size(B); 
l=size(C,1); 
r=length(Kbar); % relative degree 
%% Calculation of K_i using (11)
M=kron(fliplr(vander([0 -q(:)'*h]))*(diag(1./factorial(0:r-1))),eye(l)); % M from (8)
Km=[Kbar{:}]*M^(-1); % (11) 
K0=Km(:,1:l); 
K=cell(1,r-1); 
for i=1:r-1
    K{i}=Km(:,l*i+1:l*i+l); % K{i}=K_i
end
%% Decision variables 
P=sdpvar(n); 
W0=sdpvar(m); % Phi<=0 implies W0, W{i}, R{i} >= 0
W=cell(1,r-1); 
R=cell(1,r-1); 
for i=1:r-1
    W{i}=sdpvar(m); 
    R{i}=sdpvar(m); 
end
%% Notations 
Cbar=obsv(A,C);         % observability matrix 
Cbar=Cbar(1:r*l,:);     % \bar{C} from (8)
D=A+B*[K0,K{:}]*M*Cbar; % D from (9)

H=zeros(m); 
for i=1:r-1
    H=H+(q(i)*h)^r*K{i}'*R{i}*K{i}; 
end 

Wsum=K0'*W0*K0; 
for i=1:r-1
    Wsum=Wsum+K{i}'*W{i}*K{i}; 
end
%% LMI
Phi=blkvar; 
Phi(1,1)=P*D+D'*P+2*alpha*P+h^2*exp(2*alpha*h)*(C*A)'*Wsum*(C*A); 
Phi(1,2)=kron(ones(1,r),P*B); 
Phi(1,3)=kron(ones(1,r-1),P*B); 
Phi(1,4)=(C*A^(r-1)*D)'*H; 
Phi(2,2)=-pi^2/4*blkdiag(W0,W{:})*kron(diag(exp(-2*alpha*[0 q(:)']*h)),eye(m)); 
Phi(2,4)=kron(ones(r,1),(C*A^(r-1)*B)'*H); 
Phi(3,3)=-(factorial(r))^2*blkdiag(R{:})*kron(diag(exp(-2*alpha*q*h)./(q*h).^r),eye(m)); 
Phi(3,4)=kron(ones(r-1,1),(C*A^(r-1)*B)'*H); 
Phi(4,4)=-H; 
Phi=sdpvar(Phi); 
%% Solution of LMIs
LMIs=[P>=0, Phi<=0]; 
options=sdpsettings('solver','lmilab','verbose',0); 
sol=optimize(LMIs,[],options); 

if sol.problem == 0
    primal=check(LMIs); 
    if min(primal)<=0 
        K0=[]; 
    end
else
    K0=[]; 
    yalmiperror(sol.problem)
end