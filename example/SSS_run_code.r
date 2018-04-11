#RUN SSS
Dir.in<-"C:/Users/"
source(paste(Dir.in,"/SSS_code_newSRs.r",sep=""))
#load(paste(Dir.in,"/Dep_sdmod_spp_in.DMP",sep=""))

POP.SSS.BH<-SSS(filepath=paste(Dir.in,"/sssexample_BH",sep=""),
                      file.name=c("simple_pop.dat","simple_pop.ctl"),
                      reps=10,
                      seed.in=19,
                      M.in=c(3,0.035,0.4,3,0.037,0.4),
                      Dep.in=c(2,0.32,0.2),
                      SR_type=3,
                      h.in=c(2,0.779,0.152),
                      FMSY_M.in=c(-30,0.8,0.2),
                      BMSY_B0.in=c(-2,0.4,0.05),
                      L1.in=c(0,0,0,0),
                      Linf.in=c(0,0,0,0),
                      k.in=c(0,0,0,0),
                      Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),
                      R_start=c(0,9),
                      sum_age=0,
                      sb_ofl_yrs=c(2017,2017,2018),
                      f_yr=2016,
                      year0=1918,
                      genders=T)


POP.SSS.gR<-SSS(filepath=paste(Dir.in,"/sssexample_RickPow",sep=""),
                    file.name=c("simple_pop.dat","simple_pop.ctl"),
                    reps=10,
                    seed.in=19,
                    M.in=c(3,0.035,0.4,3,0.037,0.4),
                    Dep.in=c(2,0.32,0.2),
                    SR_type=9,
                    h.in=c(2,0.779,0.152),
                    FMSY_M.in=c(30,0.8,0.2),
                    BMSY_B0.in=c(2,0.4,0.05),
                    L1.in=c(0,0,0,0),
                    Linf.in=c(0,0,0,0),
                    k.in=c(0,0,0,0),
                    Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),
                    R_start=c(0,9),
                    sum_age=0,
                    sb_ofl_yrs=c(2017,2017,2018),
                    f_yr=2016,
                    year0=1918,
                    genders=T)

#load results from file
load(paste(Dir.in,"/sssexample_BH/SSS.DMP",sep="")) 

#Check FMSY/M and BMSY/B0 draws and outputs
plot(POP.SSS.gR$Priors$"FMSY/M",POP.SSS.gR$Posteriors$"FMSY"/POP.SSS.gR$Posteriors$"M_f",pch=21,bg="gray",xlab=expression(paste(F[MSY]/M," draws",sep="")),ylab=expression(paste(F[MSY]/M," SS out",sep="")),main="generalized Ricker")
abline(a=0,b=1)       
plot(POP.SSS.gR$Priors$"BMSY/B0"/100,POP.SSS.gR$Posteriors$"BMSY/B0",pch=21,bg="gray",xlab=expression(paste(B[MSY]/B[0]," draws",sep="")),ylab=expression(paste(B[MSY]/B[0]," SS out",sep="")),xlim=c(0,1),ylim=c(0,1),main="generalized Ricker")
abline(a=0,b=1)       
