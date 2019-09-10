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

%dict for 401(k) example of Belloni et al (2017)
% b_zx=@(z,x) [1; x'; z; x'.*z];
% p_zx=2 * (size(x,2) + 1);
% b_x = @(x) [1; x'];
% p_x = size(x,2) + 1;
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
%% 401(k) example
% grid1 = pt; grid0 = pt; % pt is grid defined in the 401(k) example of Belloni et al (2017)
% [theta1_hat,theta0_hat,j_hat,omega1_hat, omega0_hat] = cntr_dist_dml_rr(grid1,grid0,y,d,z,x,p_zx,b_zx,m_y,m_d,L,id,idn,options);
% %% joint confidence band
% mu = zeros(length(grid1),1);
% % sigma = eye(length(grid));
% V = j_hat^-1 * omega1_hat*j_hat^-1; theta1_hat_var = diag(V)/N;
% sigma = diag(diag(V))^-0.5 * V * diag(diag(V))^-0.5;
% rng('default')  % For reproducibility
% R = mvnrnd(mu,sigma,10000);
% R_max = max(abs(R),[],2); % take row-wise max
% c1 = quantile(R_max, 0.95);
% %%
% mu = zeros(length(grid0),1);
% V = j_hat^-1 * omega0_hat*j_hat^-1; theta0_hat_var = diag(V)/N;
% sigma = diag(diag(V))^-0.5 * V * diag(diag(V))^-0.5;
% rng('default')  % For reproducibility
% R = mvnrnd(mu,sigma,10000);
% R_max = max(R,[],2); % take row-wise max
% c0 = quantile(R_max, 0.95);
% 
% %% plot simultaneous confidence band
% theta1_l = theta1_hat - c1 * theta1_hat_var.^(0.5); 
% theta1_u = theta1_hat + c1 * theta1_hat_var.^(0.5); 
% 
% theta0_l = theta0_hat - c0 * theta0_hat_var.^(0.5); 
% theta0_u = theta0_hat + c0 * theta0_hat_var.^(0.5); 
% 
% filename = 'ExampleNetTFA_rrr_kappa_cntr_dist_spec4';
% fig = figure('Name',strcat(filename,' CB'))
% hold on;
% 
% for idx=1:length(grid1)
%     plot([grid1(idx) grid1(idx)], [theta1_l(idx) theta1_u(idx)],'-k.');
% end
% hold on
% h1=scatter(grid1,theta1_hat,40,'x');
% 
% hold on
% 
% for idx=1:length(grid0)
%     plot([grid0(idx) grid0(idx)], [theta0_l(idx) theta0_u(idx)],'-k.');
% end
% hold on
% h2=scatter(grid0,theta0_hat,40,'*');
% hold on
% ax = gca;
% ax.FontSize = 15;
% xlabel('net financial assets','FontSize',15);
% % ylabel('CDF of counterfactual outcomes');
% legend1 = 'Estimated CDF of Y^{(1)}';
% legend2 = 'Estimated CDF of Y^{(0)}';
% legend([h1,h2],{legend1, legend2},'Location','southeast','FontSize',15);
% %%
% dir = './simulation results/';
% figurename = strcat(dir, filename, '_CB','.eps');
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,figurename,'-depsc','-r0')
% filename = strcat(dir, filename, '_CB','.csv');
% %%
% CB = [grid1' theta1_hat theta1_l theta1_u grid0' theta0_hat theta0_l theta0_u]; %concatenate
% dlmwrite(filename, CB,'delimiter',',');
% 
% %%
% CB = csvread('./simulation results/ExampleNetTFA_rrr_kappa_cntr_dist_spec4_CB.csv');
% grid1 = CB(:,1); theta1_hat = CB(:,2);  theta1_l = CB(:,3); theta1_u =  CB(:,4); 
% grid0 = CB(:,5); theta0_hat = CB(:,6); theta0_l =  CB(:,7); theta0_u =  CB(:,8);

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

