function [theta1_hat,theta0_hat,j_hat,omega1_hat, omega0_hat]=cntr_dist_dml_rr(grid1,grid0,y,d,z,x,p,b,m_y,m_d,L,id,idn,options)
%This function (ideally) loops over the grid to calculate the components of
%DML moments for cntr distribution estimates; cntr controls whether treated
%or untreated

numerator1_full = []; numerator0_full = []; denominator_full = [];
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

    numerator1=zeros(length(grid1),Nl);     numerator0=zeros(length(grid0),Nl);
	denominator=zeros(1,Nl);
    alpha_hat_vec=zeros(1,Nl);
%%
    for i=1:Nl
        alpha_hat_vec(:,i) = alpha_hat(z_l(i),x_l(i,:));
        denominator(:,i) =(m_d(x_l(i,:),gamma_d_hat)+alpha_hat_vec(:,i)*(d_l(i)-gamma_d_hat(z_l(i),x_l(i,:))));
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
           numerator1(k,i) = (m_y(x_l(i,:),gamma_y1_hat)+alpha_hat_vec(:,i)*(d_l(i)*ind_l(i)-gamma_y1_hat(z_l(i),x_l(i,:))));
        end
   end
  %%  
   parfor k = 1:length(grid0)
        ind_nl = y_nl <= grid0(k);
        ind_l = y_l <= grid0(k);
        gamma_y0_hat=get_gamma(-(d_nl-1).*ind_nl,z_nl,x_nl,p,b,options); % if doing logit, need to flip the sign to be binary and then flip back as below
        gamma_y0_hat = @(z,x) -gamma_y0_hat(z,x);

        for i=1:Nl
           numerator0(k,i) = (m_y(x_l(i,:),gamma_y0_hat)+alpha_hat_vec(:,i)*((d_l(i)-1)*ind_l(i)-gamma_y0_hat(z_l(i),x_l(i,:))));
        end
   end

    numerator1_full = [numerator1_full, numerator1]; %append till k x N
    numerator0_full = [numerator0_full, numerator0]; %append till k x N
    denominator_full = [denominator_full, denominator]; %append till k x N

end
N = Nl*L;
theta1_hat=mean(numerator1_full,2)./mean(denominator_full,2); % k x 1
theta0_hat=mean(numerator0_full,2)./mean(denominator_full,2); % k x 1
j_hat = mean(denominator_full,2); % 1 x 1
psi1_hat = numerator1_full - repmat(theta1_hat,1,N) .* denominator_full; % k x N
omega1_hat = psi1_hat * psi1_hat' / N;
psi0_hat = numerator0_full - repmat(theta0_hat,1,N) .* denominator_full; % k x N
omega0_hat = psi0_hat * psi0_hat' / N;

end

