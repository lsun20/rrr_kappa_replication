install.packages('latex2exp')
library(latex2exp)
## numerical integration to get true cntr distribution
set.seed(666)
x = runif(1000);
#grid1 = c(-1,0,1,2,3,4);
#grid0 = c(-2,-1,0,1,2,3);
grid1 = c(-2,-1.5,-1,-0.5,0,0.5); grid0 = c(-2,-1.5,-1,-0.5,0,0.5);
cntr1_theo = rep(0,length(grid1));
cntr0_theo = rep(0,length(grid0));

 
for (k in 1:length(grid1)) {
  cntr1_theo[k] = mean(pnorm(grid1[k]-2*x^2)*x)*2;
}

for (k in 1:length(grid0)) {
  cntr0_theo[k] = mean(pnorm(grid0[k]-2*x^2)*(x-1) + pnorm(grid0[k]))*2;
}
 
cntr1_theo[cntr1_theo>1] = 1; cntr0_theo[cntr0_theo>1] = 1;


algo_legend = c( 'Auto-DML', 'DML', TeX('$\\kappa$-weight'), 'truth');
folds_legend = c('5 folds','2 folds','10 folds', 'truth');
lambda_legend = c (TeX('$\\lambda$'), TeX('1/2 $\\lambda$'),TeX('2 $\\lambda$'), 'truth')

# vary the truth
cntr_theo = cntr0_theo; grid = grid0; ylabel_text = TeX('Estimated CDF of $Y^{(0)}$');
 
# file configuration
temp = read.csv(paste0(dir,filename,'.csv'));
dir = './results/dml_rrr_sim_cntr_211004/';
filename = 'untreated_lambda_L5_keep_full';  
# table configuration
filename = 'untreated_keep'
tablename = paste0(dir, filename, '.csv');

# figure configuration
legend_text = algo_legend
filename = 'untreated_lambda_L5_keep'; 
figurename = paste0(dir, filename, '.pdf');

# read in the simulation results
theta_rrr_keep = temp[,1:(length(grid0))];
theta_rrr_keep_lb = temp[,(1+length(grid0)):(2*length(grid0))];
theta_rrr_keep_ub = temp[,(1+2*length(grid0)):(3*length(grid0))];
coverage = t(t(theta_rrr_keep_lb) <= cntr_theo & cntr_theo <= t(theta_rrr_keep_ub)); 

theta_dr_plugin_keep  = temp[,(1+4*length(grid0)):(5*length(grid0))];
theta_kappa_keep = temp[,(1+5*length(grid0)):(6*length(grid0))];

# export summary statistics
results = matrix(NA,10,length(grid)+1);
results[1,] = c(colMeans(theta_rrr_keep) - cntr_theo, 
                    mean(colMeans(theta_rrr_keep) - cntr_theo));
results[2,] = c(colMeans(theta_dr_plugin_keep) - cntr_theo,
                    mean(colMeans(theta_dr_plugin_keep) - cntr_theo));
results[3,] = c(colMeans(theta_kappa_keep) - cntr_theo,
                    mean(colMeans(theta_kappa_keep) - cntr_theo));
results[4,] = c(diag(var(theta_rrr_keep)), 
                    mean(diag(var(theta_rrr_keep))));
results[5,] = c(diag(var(theta_dr_plugin_keep)), 
                          mean(diag(var(theta_dr_plugin_keep))));
results[6,] = c(diag(var(theta_kappa_keep)), 
                         mean(diag(var(theta_kappa_keep))));
rmse <- function(var) {
  colMeans(sweep(var, 2, cntr_theo)^2)^0.5
}
results[7,] = c(rmse(theta_rrr_keep), mean(rmse(theta_rrr_keep)));
results[8,] = c(rmse(theta_dr_plugin_keep), mean(rmse(theta_dr_plugin_keep)));
results[9,] = c(rmse(theta_kappa_keep), mean(rmse(theta_kappa_keep)));
results[10,] =c(colMeans(theta_rrr_keep_ub - theta_rrr_keep_lb),
                mean(apply(coverage,1,min)));

write.csv(round(results,3),tablename)

# export quantiles of simulations
iqr_lb <- function(var) {
  quantile(var,0.1)
}
iqr_ub <- function(var) {
  quantile(var,0.9)
}
ymin = min(apply(theta_rrr_keep, 2,iqr_lb),
           apply(theta_dr_plugin_keep, 2,iqr_lb),
           apply(theta_kappa_keep, 2,iqr_lb))  

ymax = max(apply(theta_rrr_keep, 2,iqr_ub),
           apply(theta_dr_plugin_keep, 2,iqr_ub),
           apply(theta_kappa_keep, 2,iqr_ub))  

dev.off()

pdf(figurename)
cex <- 1.5; sp <- 0.1;
par(mfrow=c(1,1),mar=c(5,5,1,1) ,  cex.lab=cex, cex.axis=cex, cex.main=cex)
plot(1, type="n",xlim=c(min(grid)-0.5,max(grid)+0.5),
     xaxt="n",ylim = c(ymin, ymax),xlab="y",ylab=ylabel_text)

segments(grid+2*sp,apply(theta_dr_plugin_keep, 2,iqr_lb),grid+2*sp,
         apply(theta_dr_plugin_keep, 2,iqr_ub),lty=2, lwd=3,col="blue")

segments(grid+sp ,apply(theta_rrr_keep, 2,iqr_lb),grid+sp ,
         apply(theta_rrr_keep, 2,iqr_ub),lty=1, lwd=3,col="red")

segments(grid+3*sp,apply(theta_kappa_keep, 2,iqr_lb),grid+3*sp,
         apply(theta_kappa_keep, 2,iqr_ub),lty=3, lwd=3,col="darkgreen")

points(grid,cntr_theo, cex = cex,pch =15)
points(grid+sp,apply(theta_rrr_keep, 2,median), cex = cex,pch =16,col="red")
points(grid+2*sp,apply(theta_dr_plugin_keep, 2,median), cex = cex,pch =17,col="blue")
points(grid+3*sp,apply(theta_kappa_keep, 2,median), cex = cex,pch =19,col="green")

axis(1, at=grid, labels=grid)
legend("topleft",legend = legend_text ,lty=c(1,2,3,NA),lwd = 4,
       pch=c(NA,NA,NA,19), col=c('red','blue','darkgreen',1),cex=cex)
dev.off()
