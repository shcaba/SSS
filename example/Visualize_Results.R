#########################################################################################
#
#           VISUALIZE THE RESULTS FROM SIMPLE STOCK SYNTHESIS (SSS)
#
#########################################################################################

#load results from file
load("C:/SSS/sssexample_BH/SSS_out.DMP") 
load("C:/SSS/sssexample_RickPow/SSS_out.DMP")

# plot the results
names(SSS.out)

# histogram of parameter draws
par(mfrow = c(2,2))
hist(SSS.out$Priors$h, main = "h", xlim = c(0.2, 1))
hist(SSS.out$Priors$Dep, main = "Depletion", xlim = c(0, 1))
hist(SSS.out$Priors$M_f, main = "M (females)", xlim = c(0.02,0.08))
hist(SSS.out$Priors$M_m, main = "M (males)", xlim = c(0.02,0.08))

# Only plots if the Ricker Power function is used:
par(mfrow = c(1,2))
hist(SSS.out$Priors$`FMSY/M`, main = "Fmsy/M", xlim = c(0.5, 1.5))
hist(SSS.out$Priors$`BMSY/B0`/100, main = "Bmsy/B0", xlim = c(0.2, 0.8))


# Check FMSY/M and BMSY/B0 draws and outputs
par(mfrow = c(1,2))
plot(SSS.out$Priors$"FMSY/M", SSS.out$Posteriors$"FMSY"/ SSS.out$Posteriors$"M_f",pch=21,bg="gray",xlab=expression(paste(F[MSY]/M," draws",sep="")),ylab=expression(paste(F[MSY]/M," SS out",sep="")),main="generalized Ricker")
abline(a=0,b=1)       
plot(SSS.out$Priors$"BMSY/B0"/100, SSS.out$Posteriors$"BMSY/B0",pch=21,bg="gray",
     xlab=expression(paste(B[MSY]/B[0]," draws",sep="")),ylab=expression(paste(B[MSY]/B[0]," SS out",sep="")),xlim=c(0,1),ylim=c(0,1), main="generalized Ricker")
abline(a=0,b=1)    

# boxplots of the projected OFL values
par(mfrow = c(1,2))
ymax = max(SSS.out$Posteriors$OFL_2018, SSS.out$Posteriors$OFL_2019, SSS.out$Posteriors$AdjCatch_2018, SSS.out$Posteriors$AdjCatch_2019)
ymin = min(SSS.out$Posteriors$OFL_2018, SSS.out$Posteriors$OFL_2019, SSS.out$Posteriors$AdjCatch_2018, SSS.out$Posteriors$AdjCatch_2019)

boxplot(cbind(SSS.out$Posteriors$OFL_2018, SSS.out$Posteriors$OFL_2019), ylim = c(ymin, ymax), main = "OFL", xlab = "")
boxplot(cbind(SSS.out$Posteriors$AdjCatch_2018, SSS.out$Posteriors$AdjCatch_2019), ylim = c(ymin, ymax),main = "ACL", xlab = "")


# plot the summary biomass trajectories
par(mfrow = c(1,2))
yr = as.numeric(colnames(SSS.out$Summary_Biomass))
n = dim(SSS.out$Summary_Biomass)[1]
plot(yr, SSS.out$Summary_Biomass[1,], lwd = 2, type = 'l', 
     ylim = c(0, max(SSS.out$Summary_Biomass)), ylab = "Summary Biomass", xlab = "Year")
for (nn in 1:n){ lines(yr, SSS.out$Summary_Biomass[nn,], col = 'grey')}
lines(yr, apply(SSS.out$Summary_Biomass,2,median), lwd = 2)

# plot the relative stock status
yr = as.numeric(colnames(SSS.out$Rel_Stock_status_series))
n = dim(SSS.out$Rel_Stock_status_series)[1]
plot(yr, SSS.out$Rel_Stock_status_series[1,], lwd = 2, type = 'l', ylim = c(0, 1), 
     xlab = "Year", ylab = "Relative Stock Status")
for (nn in 1:n){ lines(yr, SSS.out$Rel_Stock_status_series[nn,], col = 'grey')}
lines(yr, apply(SSS.out$Rel_Stock_status_series,2,median), lwd = 2)