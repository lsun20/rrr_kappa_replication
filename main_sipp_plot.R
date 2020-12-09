library(latex2exp)
# load data
dir <- './results/sipp_cntr_200216/'
# cntr_CB_nn_full 200216 uses full sample and learns gamma by nn 
# cntr_CB_nn 200106 pre-processes and learns gamma by nn 
# cntr_CB 200105 pre-processes and learns gamma  lasso

spec <- 'treated'
filename <- paste(dir,spec,'_CB_nn','.csv',sep = "")
cntr1 = read.csv(filename)
spec <- 'untreated'
filename <- paste(dir,spec,'_CB_nn','.csv',sep = "")
cntr0 = read.csv(filename)
# plotting
spec <- 'cntr'
figurename <- paste(dir,spec,'_CB_nn_full','.pdf',sep = "")
dev.off()
pdf(figurename)
cex <- 2
par(mfrow=c(1,1),mar=c(5,5,1,1) ,  cex.lab=cex, cex.axis=cex, cex.main=cex)

# plot(cntr1$perc,cntr1$theta_hat, xlab="net financial assets",ylab="CDF of counterfactual outcome",)

plot(range(cntr1$perc), range(0,1),  type="n",col=1, xlab="net financial assets", lwd=2, cex.lab=2,
     ylab="", 
     main="",
     sub=" ", xaxt = "n");

polygon(c(cntr1$perc,rev(cntr1$perc)),c(cntr1$theta_u,rev(cntr1$theta_l)), density=100, border=F, 
        col='light blue', lty = 1, lwd = 1);
lines(cntr1$perc, cntr1$theta_hat, lty = 1, lwd=1, col = 4 , xaxt = "n" );

polygon(c(cntr0$perc,rev(cntr0$perc)),c(cntr0$theta_u,rev(cntr0$theta_l)), density=100, border=F, 
        col='light green', lty = 1, lwd = 1);
lines(cntr0$perc, cntr0$theta_hat,  lty = 2, lwd=1, col = 4  , xaxt = "n");

axis(1, at=seq(0,100000,20000), labels=c('0','20k','40k','60k','80k','100k'))

# legend('bottomright', c(' ', ' '), col = c('light blue','light green'), lwd = c(4,4), horiz = F, bty = 'n',
#        cex=cex);
legend('bottomright',
       c(TeX('Estimated CDF of $\\Y^{(1)}$'), TeX('Estimated CDF of $\\Y^{(0)}$')), 
       col = c('light blue','dark green'), lty = c(1,2), lwd = c(3,3), horiz = F, bty = 'n',
       cex=cex);

dev.off()

# check whether our estimated CDFs are monotone:
all(cntr1$theta_hat == cummax(cntr1$theta_hat))
all(cntr0$theta_hat == cummax(cntr0$theta_hat))
all(cntr0$theta_l == cummax(cntr0$theta_l))
all(cntr0$theta_u == cummax(cntr0$theta_u))
all(cntr1$theta_u == cummax(cntr1$theta_u))
all(cntr1$theta_l == cummax(cntr1$theta_l))