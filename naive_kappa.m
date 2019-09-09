function theta_tilde=naive_kappa(y,d,z,x,p,b,options)
%%
N = length(y);
N_train = floor(0.8*N); N_valid = N - N_train;
%reg param
r_prop_grid = [0, 0.05, 0.1,0.2];


%at data
B=zeros(N,p);
B_train=zeros(N_train,p); B_valid=zeros(N_valid,p);
N_z_train=zeros(p,N_train); N_z_valid=zeros(p,N_valid);

for i=1:N
    if i <= N_train
        B_train(i,:)=b(x(i))';
        N_z_train(:,i)=z(i).*b(x(i));
    else
        j = i-N_train;
        B_valid(j,:)=b(x(i))';
        N_z_valid(:,j)=z(i).*b(x(i));
    end
    B(i,:)=b(x(i))';
end

G_hat_train=B_train'*B_train./N_train; G_hat_valid=B_valid'*B_valid./N_valid;
N_z_hat_train=mean(N_z_train,2); N_z_hat_valid=mean(N_z_valid,2);
  
%kappa reweighting via Matlab logit and no selection
Beta1 = glmfit(B,z,'binomial','link','logit');          
pi_hat = glmval(Beta1,B,'logit');
alpha_hat = (z - pi_hat)./(pi_hat.*(1-pi_hat));
num_tilde = y.* alpha_hat;
denom_tilde = d.*alpha_hat;


%kappa reweighting via Matlab lasso + logit 
% [Beta,FitInfo] = lassoglm(B,z,'binomial','NumLambda',5,'CV',2);
% indx = FitInfo.Index1SE;
% Beta0 = Beta(:,indx);
% cnst = FitInfo.Intercept(indx);
% Beta1 = [cnst;Beta0];
% pi_hat = glmval(Beta1,B,'logit');
% alpha_hat = (z - pi_hat)./(pi_hat.*(1-pi_hat));
% num_tilde = y.* alpha_hat;
% denom_tilde = d.*alpha_hat;
%kappa reweighting via manual lasso + linear
% %init
% pi0=zeros(p,1);
% pi_hat = lasso_grid(pi0, G_hat_train,G_hat_valid,N_z_hat_train,N_z_hat_valid,r_prop_grid,options);
% alpha_hat=@(z,x) (z - b(x)'*pi_hat)/(b(x)'*pi_hat)/(1-b(x)'*pi_hat);
% kappa_hat=@(w,z,x,alpha_hat) w * alpha_hat(z,x); % w is the appropriate r.v
% 
% num_tilde=zeros(1,N);
% denom_tilde=zeros(1,N);
% 
% for i=1:N
%    num_tilde(i)=kappa_hat(y(i),z(i),x(i),alpha_hat);
%    denom_tilde(i)=kappa_hat(d(i),z(i),x(i),alpha_hat);
% 
% end
%%
% theta_tilde=mean(num_tilde)/mean(denom_tilde);

%create an index for trimming
keep = (pi_hat <= 0.9 & pi_hat >= 0.1 );
theta_tilde = mean(num_tilde(keep))/mean(denom_tilde(keep));

end


