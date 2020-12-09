%% set up

clear;


% choose f, simulate y
fy=@(z,x) 0 + 2*x.^2.*z;  
fd=@(z,x) 0 + (1*x).*z;  
% px=@(x) 0.2*(x<0.5)+0.8*(x>0.5);
% px=@(x) 0.2 + 0.5*x;
% px=@(x) 0.5;
% px=@(x) 0.1*(x<0.5)+0.9*(x>0.5);
% px=@(x) 0.5./exp(x);
% px=@(x) 0.1 + 0.5*x + 0.3*x.^2;
px=@(x) 0.05*(x<0.5)+0.95*(x>0.5);
filename = 'dr_plugin_2_1_005095x_test';
 
%% universal
rng('default')
options = optimset('Display','off','Algorithm','interior-point');
warning('off','all')
%dict for simulation
b_zx=@(z,x) [ones(length(x),1); x; x.^2; x.^3; x.^4;...
                z; x.*z; x.^2.*z;  x.^3.*z; x.^4.*z];

p_zx=2 * (4 + 1);
b_x = @(x) [ones(length(x),1); x; x.^2; x.^3; x.^4];
p_x = 5; %dim of dict
 
%functional
m_y=@(x,gamma_y) gamma_y(1,x)-gamma_y(0,x);
m_d=@(x,gamma_d) gamma_d(1,x)-gamma_d(0,x);

%influence
psi=@(y,d,z,x,m_y,m_d,alpha,gamma_y,gamma_d) (m_y(x,gamma_y)+alpha(z,x)*(y-gamma_y(z,x)))/(m_d(x,gamma_d)+alpha(z,x)*(d-gamma_d(z,x)));

%general lasso objective
lasso_obj_param=@(rho,G,M,r) rho'*G*rho-2.*M'*rho+2.*r.*norm(rho,1);
obj=lasso_obj_param;

%%set up sample splitting ids
N=1000;
% N = n; % for 401(k)
L=5;
Nl=N/L;

id=zeros(N,L);
for l=1:5
    id((l-1)*Nl+1:(l*Nl),l)=1;
end
%%
idn=1-id;
idn=idn.*(1:N)';
idn=nonzeros(idn);
%idn=vec2mat(idn,800)';
idn=reshape(idn,N-Nl,5);

id=id.*(1:N)';
id=nonzeros(id);
%id=vec2mat(id,200)';
id=reshape(id,Nl,5);

%% simulations

N_sim=500;

%% for counterfaction distribution
grid1 = [-1,0,1,2,3,4];
grid0 = [-2,-1,0,1,2,3];
theta1_rrr=zeros(N_sim,length(grid1));
theta0_rrr=zeros(N_sim,length(grid0));
theta1_plugin_trim1=zeros(N_sim,length(grid1));
theta0_plugin_trim1=zeros(N_sim,length(grid0));
theta1_plugin_trim2=zeros(N_sim,length(grid1));
theta0_plugin_trim2=zeros(N_sim,length(grid0));

theta1_naive_kappa=zeros(N_sim,length(grid1));
theta0_naive_kappa=zeros(N_sim,length(grid0));
%%
for n=1:N_sim
    [n/N_sim]
    [y,d,z,x]=sim_rr(fy,fd,px,N,n);
    %%
% kappa-weight
    parfor k=1:length(grid1)
        ind = y <= grid1(k);
        theta1_naive_kappa(n,k) = naive_kappa((d).*ind,d,z,x,p_x,b_x,options);%manual kappa reweighting
    end
    parfor k=1:length(grid0)
        ind = y <= grid0(k);
        theta0_naive_kappa(n,k) = naive_kappa((d-1).*ind,d-1,z,x,p_x,b_x,options);%manual kappa reweighting
    end
% DML-RRR (need to change get_alpha function to use the RRR approach)
%     [theta1_rrr(n,:),theta0_rrr(n,:),~,~,~] = cntr_dist_dml_rr(grid1,grid0,y,d,z,x,p_zx,b_zx,m_y,m_d,L,id,idn,options);
     
% DML-plugin with different trimming set (need to change get_alpha
% function to use the plugin approach)
%     [theta1_plugin_trim1(n,:),theta0_plugin_trim1(n,:),~,~,~,...
%       theta1_plugin_trim2(n,:),theta0_plugin_trim2(n,:),~,~,~] = cntr_dist_dml_rr_trim(grid1,grid0,y,d,z,x,p_zx,b_zx,m_y,m_d,L,id,idn,options);
end
%%
% filename = 'dr_plugin_2_1_005095x_cntr1_dist_logit_trim';
theta1_hat = theta1_naive_kappa;
% theta1_hat = [theta1_rrr, theta1_naive_kappa]; %concatenate
% theta1_hat = theta1_rrr;
% theta1_hat = [theta1_plugin_trim1, theta1_plugin_trim2]; %concatenate

dlmwrite(strcat(filename,'.csv'), theta1_hat,'-append','delimiter',',');

% filename = 'dr_plugin_2_1_005095x_cntr0_dist_logit_trim';
theta0_hat = theta0_naive_kappa;
% theta0_hat = [theta0_rrr, theta0_naive_kappa]; %concatenate
% theta0_hat = theta0_rrr;
% theta0_hat = [theta0_plugin_trim1, theta0_plugin_trim2]; %concatenate

dlmwrite(strcat(filename,'.csv'), theta0_hat,'-append','delimiter',',');

