# load data
dir <- './results/sipp_cntr_191116/'
spec <- 'treated'
filename <- paste(dir,spec,'_CB','.csv',sep = "")
cntr1 = read.csv(filename)
spec <- 'untreated'
filename <- paste(dir,spec,'_CB','.csv',sep = "")
cntr0 = read.csv(filename)
# plotting
dir <- './results/sipp_cntr_191116/'
spec <- 'cntr'
figurename <- paste(dir,spec,'_CB','.pdf',sep = "")
dev.off()
pdf(figurename)
# plot(cntr1$perc,cntr1$theta_hat, xlab="net financial assets",ylab="CDF of counterfactual outcome",)

plot(range(cntr1$perc), range(0,1),  type="n",col=1, xlab="net financial assets", lwd=2,
     ylab="CDF of counterfactual outcomes", 
     main="",
     sub=" ");

polygon(c(cntr1$perc,rev(cntr1$perc)),c(cntr1$theta_u,rev(cntr1$theta_l)), density=100, border=F, 
        col='light blue', lty = 1, lwd = 1);
lines(cntr1$perc, cntr1$theta_hat, lty = 1, lwd=1, col = 4  );

polygon(c(cntr0$perc,rev(cntr0$perc)),c(cntr0$theta_u,rev(cntr0$theta_l)), density=100, border=F, 
        col='light green', lty = 1, lwd = 1);
lines(cntr0$perc, cntr0$theta_hat,  lty = 2, lwd=1, col = 4  );


legend(min(cntr0$perc), 1, c(' ', ' '), col = c('light blue','light green'), lwd = c(4,4), horiz = F, bty = 'n');
legend(min(cntr0$perc), 1, c('W/ 401(k)', 'W/o 401(k)'), col = c('light blue','dark green'), lty = c(1,2), lwd = c(1,1), horiz = F, bty = 'n');

dev.off()