
function rho_hat=lasso_grid(rho0, G_train,G_valid,M_train,M_valid,r_grid,options)
%searches over r_grid for the optimal l1 regulariation parameter

r_grid_len = length(r_grid);
p = length(rho0); %dim of dict

mse_valid = @(rho) rho'*G_valid*rho-2.*M_valid'*rho;
    
%init
mse = zeros(r_grid_len,1);
rho_hat_grid = zeros(p,r_grid_len);
for i=1:r_grid_len
    
    r = r_grid(i);

    %empirical lasso objective
    lasso_obj =  @(rho) rho'*G_train*rho-2.*M_train'*rho+2.*r.*norm(rho,1);

    %lasso estimators
    [rho_hat,~]= fminunc(lasso_obj,rho0,options);
    rho_hat_grid(:,i) = rho_hat;
    mse(i) = mse_valid(rho_hat);
    
end

[~,I]=min(mse); rho_hat = rho_hat_grid(:,I); %optimize over reg parameter

end

