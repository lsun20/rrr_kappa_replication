function [theta1_hat_trim1,theta0_hat_trim1,j_hat_trim1,omega1_hat_trim1, omega0_hat_trim1,...
    theta1_hat_trim2,theta0_hat_trim2,j_hat_trim2,omega1_hat_trim2, omega0_hat_trim2]=...
    cntr_dist_dml_rr_trim(grid1,grid0,y,d,z,x,p,b,m_y,m_d,L,id,idn,options)
%This function (ideally) loops over the grid to calculate the components of
%DML moments for cntr distribution estimates; cntr controls whether treated
%or untreated

numerator1_full_trim1 = []; numerator0_full_trim1 = []; denominator_full_trim1 = [];
numerator1_full_trim2 = []; numerator0_full_trim2 = []; denominator_full_trim2 = [];
alpha_hat_trim1 = []; alpha_hat_trim2 = [];
%%
for l=1:L
 %%   
    id_nl=idn(:,l);
    id_l=id(:,l);

    %%%stage 1 for data NOT in I_l i.e. rows indexed by idn_l
    y_nl=y(id_nl);
    d_nl=d(id_nl);
    z_nl=z(id_nl);
    x_nl=x(id_nl,:);
    y_l=y(id_l);
    d_l=d(id_l);
    z_l=z(id_l);
    x_l=x(id_l,:);
   
    Nl=length(y_l);
%%
    alpha_hat=get_alpha(z_nl,x_nl,p,b,m_y,options);
    gamma_d_hat=get_gamma(d_nl,z_nl,x_nl,p,b,options);

%%
    %%%stage 2 for data in I_l i.e. rows indexed by id_l

    alpha_hat_vec=zeros(1,Nl);
%%
    for i=1:Nl
        alpha_hat_vec(:,i) = alpha_hat(z_l(i),x_l(i,:));
    end
%% trimming [0.01,0.99]
    alpha_hat_vec_trim1 = alpha_hat_vec;
    alpha_hat_vec_trim1(alpha_hat_vec<1/(1-10^-6)) & alpha_hat_vec>0)  = 1/(1-10^-6));
    alpha_hat_vec_trim1(alpha_hat_vec>(1/10^-6))  = (1/10^-6);
    alpha_hat_vec_trim1(alpha_hat_vec>-1/(1-10^-6)) & alpha_hat_vec<0)  = -1/(1-10^-6));
    alpha_hat_vec_trim1(alpha_hat_vec<-(1/10^-6))  = -(1/10^-6);
    
%% trimming [0.05,0.95]   
    alpha_hat_vec_trim2 = alpha_hat_vec;
    alpha_hat_vec_trim2(alpha_hat_vec<1/(1-10^-12) & alpha_hat_vec>0)  = 1/(1-10^-12);
    alpha_hat_vec_trim2(alpha_hat_vec>(1/10^-12))  = (1/10^-12);
    alpha_hat_vec_trim2(alpha_hat_vec>-1/(1-10^-12) & alpha_hat_vec<0)  = -1/(1-10^-12);
    alpha_hat_vec_trim2(alpha_hat_vec<-(1/10^-12))  = -(1/10^-12);
%%

    numerator1_trim1=zeros(length(grid1),Nl);     numerator0_trim1=zeros(length(grid0),Nl);
	denominator_trim1=zeros(1,Nl);
    numerator1_trim2=zeros(length(grid1),Nl);     numerator0_trim2=zeros(length(grid0),Nl);
	denominator_trim2=zeros(1,Nl);
   for i=1:Nl
        denominator_trim1(:,i) =(m_d(x_l(i,:),gamma_d_hat)+alpha_hat_vec_trim1(:,i)*(d_l(i)-gamma_d_hat(z_l(i),x_l(i,:))));
        denominator_trim2(:,i) =(m_d(x_l(i,:),gamma_d_hat)+alpha_hat_vec_trim2(:,i)*(d_l(i)-gamma_d_hat(z_l(i),x_l(i,:))));
    end
  %%         
   parfor k = 1:length(grid1)
       %%
        ind_nl = y_nl <= grid1(k);
        ind_l = y_l <= grid1(k);
        gamma_y1_hat=get_gamma(d_nl.*ind_nl,z_nl,x_nl,p,b,options);

        for i=1:Nl
            if floor(i/100) == i/100
                disp(i);
            end
           numerator1_trim1(k,i) = (m_y(x_l(i,:),gamma_y1_hat)+alpha_hat_vec_trim1(:,i)*(d_l(i)*ind_l(i)-gamma_y1_hat(z_l(i),x_l(i,:))));
           numerator1_trim2(k,i) = (m_y(x_l(i,:),gamma_y1_hat)+alpha_hat_vec_trim2(:,i)*(d_l(i)*ind_l(i)-gamma_y1_hat(z_l(i),x_l(i,:))));
        end
   end
  %%  
   parfor k = 1:length(grid0)
        ind_nl = y_nl <= grid0(k);
        ind_l = y_l <= grid0(k);
        gamma_y0_hat=get_gamma(-(d_nl-1).*ind_nl,z_nl,x_nl,p,b,options); % if doing logit, need to flip the sign to be binary and then flip back as below
        gamma_y0_hat = @(z,x) -gamma_y0_hat(z,x);

        for i=1:Nl
           numerator0_trim1(k,i) = (m_y(x_l(i,:),gamma_y0_hat)+alpha_hat_vec_trim1(:,i)*((d_l(i)-1)*ind_l(i)-gamma_y0_hat(z_l(i),x_l(i,:))));
           numerator0_trim2(k,i) = (m_y(x_l(i,:),gamma_y0_hat)+alpha_hat_vec_trim2(:,i)*((d_l(i)-1)*ind_l(i)-gamma_y0_hat(z_l(i),x_l(i,:))));
        end
   end

    numerator1_full_trim2 = [numerator1_full_trim2, numerator1_trim2]; %append till k x N
    numerator0_full_trim2 = [numerator0_full_trim2, numerator0_trim2]; %append till k x N
    denominator_full_trim2 = [denominator_full_trim2, denominator_trim2]; %append till k x N

    numerator1_full_trim1 = [numerator1_full_trim1, numerator1_trim1]; %append till k x N
    numerator0_full_trim1 = [numerator0_full_trim1, numerator0_trim1]; %append till k x N
    denominator_full_trim1 = [denominator_full_trim1, denominator_trim1]; %append till k x N

    alpha_hat_trim1 = [alpha_hat_trim1, alpha_hat_vec_trim1]; %append till k x N
    alpha_hat_trim2 = [alpha_hat_trim2, alpha_hat_vec_trim2]; %append till k x N
   
end
%%
N = Nl*L;

numerator1_full_trim1_final = ...
    numerator1_full_trim1(:,abs(alpha_hat_trim1)<(1/10^-6) & abs(alpha_hat_trim1)>1/(1-10^-6));
numerator0_full_trim1_final = ...
    numerator0_full_trim1(:,abs(alpha_hat_trim1)<(1/10^-6) & abs(alpha_hat_trim1)>1/(1-10^-6));
denominator_full_trim1_final = ...
    denominator_full_trim1(:,abs(alpha_hat_trim1)<(1/10^-6) & abs(alpha_hat_trim1)>1/(1-10^-6));

numerator1_full_trim2_final = ...
    numerator1_full_trim2(:,abs(alpha_hat_trim2)<(1/10^-12) & abs(alpha_hat_trim2)>1/(1-10^-12));
numerator0_full_trim2_final = ...
    numerator0_full_trim2(:,abs(alpha_hat_trim2)<(1/10^-12) & abs(alpha_hat_trim2)>1/(1-10^-12));
denominator_full_trim2_final = ...
    denominator_full_trim2(:,abs(alpha_hat_trim2)<(1/10^-12) & abs(alpha_hat_trim2)>1/(1-10^-12));
%%
N = size(denominator_full_trim1_final,2); %size after trimming
theta1_hat_trim1=mean(numerator1_full_trim1_final,2)./mean(denominator_full_trim1_final,2); % k x 1
theta0_hat_trim1=mean(numerator0_full_trim1_final,2)./mean(denominator_full_trim1_final,2); % k x 1
j_hat_trim1 = mean(denominator_full_trim1_final,2); % 1 x 1
psi1_hat_trim1 = numerator1_full_trim1_final - repmat(theta1_hat_trim1,1,N) .* denominator_full_trim1_final; % k x N
omega1_hat_trim1 = psi1_hat_trim1 * psi1_hat_trim1' / N;
psi0_hat_trim1 = numerator0_full_trim1_final - repmat(theta0_hat_trim1,1,N) .* denominator_full_trim1_final; % k x N
omega0_hat_trim1 = psi0_hat_trim1 * psi0_hat_trim1' / N;

N = size(denominator_full_trim2_final,2); %size after trimming
theta1_hat_trim2=mean(numerator1_full_trim2_final,2)./mean(denominator_full_trim2_final,2); % k x 1
theta0_hat_trim2=mean(numerator0_full_trim2_final,2)./mean(denominator_full_trim2_final,2); % k x 1
j_hat_trim2 = mean(denominator_full_trim2_final,2); % 1 x 1
psi1_hat_trim2 = numerator1_full_trim2_final - repmat(theta1_hat_trim2,1,N) .* denominator_full_trim2_final; % k x N
omega1_hat_trim2 = psi1_hat_trim2 * psi1_hat_trim2' / N;
psi0_hat_trim2 = numerator0_full_trim2_final - repmat(theta0_hat_trim2,1,N) .* denominator_full_trim2_final; % k x N
omega0_hat_trim2 = psi0_hat_trim2 * psi0_hat_trim2' / N;
end

