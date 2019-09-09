function [gamma_hat]=get_gamma(y,z,x,p,b,options)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% %reg param
% r_cef_y_grid=[0,0.05,0.1];
% 
% 
Nl=length(y);
% N_train = floor(0.8*Nl); N_valid = Nl - N_train;
% 
%at data
B=zeros(Nl,p);
% N_y=zeros(p,Nl);

for i=1:Nl
   B(i,:)=b(z(i),x(i,:))';
%    N_y(:,i)=y(i).*b(z(i),x(i));
end
% 
% N_y_hat=mean(N_y,2);  N_y_hat_train=mean(N_y(:,1:N_train),2); N_y_hat_valid=mean(N_y(:,(N_train+1):Nl),2);
% G_hat=B'*B./Nl; G_hat_train = B(1:N_train,:)'*B(1:N_train,:)./N_train; G_hat_valid = B((N_train+1):Nl,:)'*B((N_train+1):Nl,:)./N_train;
% 
% %init
% beta0=zeros(p,1);
% 
% beta_hat = lasso_grid(beta0, G_hat_train,G_hat_valid,N_y_hat_train,N_y_hat_valid,r_cef_y_grid,options);
% gamma_hat=@(z,x) b(z,x)'*beta_hat;

%kappa reweighting via Matlab lasso + logit
warning('off','stats:lasso:MaxIterReached');
[Beta,FitInfo] = lassoglm(B,y,'binomial','NumLambda',5,'CV',2);
indx = FitInfo.Index1SE;
Beta0 = Beta(:,indx);
cnst = FitInfo.Intercept(indx);
Beta1 = [cnst;Beta0];
% lambda = 2.2*sqrt(Nl)*norminv(1-(.1/log(Nl))/(2*(2*p)*Nl)); % use lambda level from Victor's paper
% [Beta,FitInfo] = lassoglm(B,y,'binomial','Lambda',lambda);
% cnst = FitInfo.Intercept;
% Beta1 = [cnst;Beta];

gamma_hat =@(z,x) glmval(Beta1,b(z,x)','logit');

end

