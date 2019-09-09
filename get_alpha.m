function [alpha_hat]=get_alpha(z,x,p,b,m_y,options)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%reg param
r_rr_grid=[0,0.05,0.1];
%%
Nl=length(z);
N_train = floor(0.8*Nl); N_valid = Nl - N_train;

%at data
B=zeros(Nl,p);
M=zeros(p,Nl);

for i=1:Nl
   B(i,:)=b(z(i),x(i,:))';
   M(:,i)=m_y(x(i,:),b);

end
%estimate alpha by RRR 
% M_hat=mean(M,2); M_hat_train=mean(M(:,1:N_train),2); M_hat_valid=mean(M(:,(N_train+1):Nl),2);
% G_hat=B'*B./Nl; G_hat_train = B(1:N_train,:)'*B(1:N_train,:)./N_train; G_hat_valid = B((N_train+1):Nl,:)'*B((N_train+1):Nl,:)./N_train;
% 
% %init
% rho0=zeros(p,1);
% 
% rho_hat = lasso_grid(rho0, G_hat_train,G_hat_valid,M_hat_train,M_hat_valid,r_rr_grid,options);
% alpha_hat=@(z,x) b(z,x)'*rho_hat;

%estimate alpha by plugging in pi(x) via Matlab lasso + logit 
B = B(:,1:(p/2)); % take out the interaction terms between Z and b(X)
b_x = @(x) [ones(length(x),1); x; x.^2; x.^3; x.^4]; % need to fix later but can't double index in matlab...
[Beta,FitInfo] = lassoglm(B,z,'binomial','NumLambda',5,'CV',2);
indx = FitInfo.Index1SE;
Beta0 = Beta(:,indx);
cnst = FitInfo.Intercept(indx);
Beta1 = [cnst;Beta0];
pi_hat = @(x) glmval(Beta1,b_x(x)','logit');
alpha_hat = @(z,x) (z - pi_hat(x))./(pi_hat(x).*(1-pi_hat(x)));



end

