library(msm)
#################
### Functions ###
#################
#1 = par file
#2 = forecastfile
#3 = data file
#4 = ctl file
#5 = SS.exe
Change.SSfiles<-function(main.dir,sppname.in,file.make=c(1,1,1,1,1),starter.par="N")
  
{
  for(i in 1:length(sppname.in))
  {
    if(file.make[1]==1)
    {
      if(starter.par=="N"){file.copy(from=paste(main.dir,"/baseSSfiles/starterfiles/",sppname.in[i],"_starter.ss",sep=""),to=paste(main.dir,"/",sppname.in[i],"/starter.SS",sep=""),overwrite=TRUE)}
      if(starter.par=="Y"){file.copy(from=paste(main.dir,"/baseSSfiles/starterfilespar/",sppname.in[i],"_starter.ss",sep=""),to=paste(main.dir,"/",sppname.in[i],"/starter.SS",sep=""),overwrite=TRUE)}
    }
    if(file.make[2]==1){file.copy(from=paste(main.dir,"/baseSSfiles/forecastfiles/",sppname.in[i],"_forecast.ss",sep=""),to=paste(main.dir,"/",sppname.in[i],"/forecast.SS",sep=""),overwrite=TRUE)}
    if(file.make[3]==1){file.copy(from=paste(main.dir,"/baseSSfiles/datfiles/",sppname.in[i],".dat",sep=""),to=paste(main.dir,"/",sppname.in[i],"/",sppname.in[i],".dat",sep=""),overwrite=TRUE)}
    if(file.make[4]==1){file.copy(from=paste(main.dir,"/baseSSfiles/ctlfiles/",sppname.in[i],".ctl",sep=""),to=paste(main.dir,"/",sppname.in[i],"/",sppname.in[i],".ctl",sep=""),overwrite=TRUE)}
    if(file.make[5]==1){file.copy(from=paste(main.dir,"/baseSSfiles/SS3.exe",sep=""),to=paste(main.dir,"/",sppname.in[i],"/SS3.exe",sep=""),overwrite=TRUE)}
  }
}

rtlnorm <- function (n, mean = 0, sd = 1, lower = -Inf, upper = Inf)
{
  ret <- numeric()
  if (length(n) > 1)
    n <- length(n)
  while (length(ret) < n) {
    y <- rlnorm(n - length(ret), mean, sd)
    y <- y[y >= lower & y <= upper]
    ret <- c(ret, y)
  }
  stopifnot(length(ret) == n)
  ret
}

rbeta.ab <- function(n, m, s, a, b)
{
  # calculate mean of corresponding standard beta dist
  mu.std <- (m-a)/(b-a)
  
  # calculate parameters of std. beta with mean=mu.std and sd=s
  alpha <- (mu.std^2 - mu.std^3 - mu.std*s^2) / s^2
  beta  <- (mu.std - 2*mu.std^2 + mu.std^3 - s^2 + mu.std*s^2) / s^2
  
  # generate n draws from standard beta
  b.std <- rbeta(n, alpha, beta)
  
  # linear transformation from beta(0,1) to beta(a,b)
  b.out <- (b-a)*b.std + a
  
  return(b.out)
}

RUN.SS<-function(path, ss.exe="SS3",ss.cmd=" -nohess -nox")
{
  navigate <- paste("cd ",path,sep="")
  command <- paste(navigate," & ",ss.exe,ss.cmd,sep="")
  shell(command,invisible=TRUE,translate=TRUE)
}

SSS.jr.exSSS<-function(filepath,file.name,reps=1000,Dep.in=c(1,0.4,0.1),M.in=c(3,0.1,0.4,3,0.1,0.4),h.in=c(1,0.6,0.2),L1.in=c(0,0,0,0),Linf.in=c(0,0,0,0),k.in=c(0,0,0,0),Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),R_start=6,sum_age=0,sb_ofl_yrs=c(2010,2011,2012),f_yr=2012,year0=1916,genders=F)
{
  start.time<-Sys.time()
  Spp.quant.out<-list()
  if(h.in[1]>=0){Input.draws<-as.data.frame(matrix(NA,nrow=reps,ncol=6))}
  if(h.in[1]<0){Input.draws<-as.data.frame(matrix(NA,nrow=reps,ncol=7))}
  Input.draws.M<-as.data.frame(matrix(NA,nrow=reps,ncol=4))
  sb.years<-c(year0:sb_ofl_yrs[1])
  Quant.out<-as.data.frame(matrix(NA,nrow=reps,ncol=24))
  SB.out<-as.data.frame(matrix(NA,nrow=reps,ncol=length(sb.years)))
  colnames(SB.out)<-sb.years
  
  starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
  sum_age_line<-strsplit(starter.new[grep("summary biomass",starter.new)], " ")[[1]]
  sum_age_line[1]<-sum_age
  starter.new[grep("summary biomass",starter.new)]<-paste(sum_age_line, collapse=" ")
  write(starter.new,paste(filepath,"/starter.ss",sep=""))
  i<-1
  ii<-1
  while(i<reps+1)
  {
    print(i)
    print(ii)
    if(ii%%10==0){print(file.name[1])}
    ### Draw random values ###
    #-99 = no estimation
    #0 = normal
    #10 = truncated normal
    #1 = symmetric beta
    #2 = full beta (rbeta)
    #3 = lognormal
    #30 = truncated lognormal
    #4 = uniform
    if(Dep.in[1]==2){Dep.draw<-round(1-rbeta.ab(1,1-Dep.in[2],Dep.in[3],0.01,0.99),2)}
    if(Dep.in[1]==10){Dep.draw<-round(rtnorm(1,Dep.in[2],Dep.in[3],0.01,1),2)}
    if(Dep.in[1]==10){Input.draws[i,2]<-Dep.draw}
    #Steenpess
    if(h.in[1]>=0)
    {
      if(h.in[1]==2){h.draw<-round(rbeta.ab(1,h.in[2],h.in[3],0.25,0.99),2)}
      if(h.in[1]==10){h.draw<-round(rtnorm(1,h.in[2],h.in[3],0.25,0.99),2)}
      if(h.in[1]==30){h.draw<-round(rtlnorm(1,h.in[2],h.in[3],0.25,0.99),2)}
      if(h.in[1]==4){h.draw<-round(runif(1,h.in[2],h.in[3]),2)}      
      Input.draws[i,1]<-h.draw
    }
    #Three parameter S-R model
    if(Zfrac.Beta.in[1]>=0)
    {
      if(Zfrac.Beta.in[1]==10){Zfrac.draw<-round(rtnorm(1,Zfrac.Beta.in[2],Zfrac.Beta.in[3],0,1),2)}
      if(Zfrac.Beta.in[1]==30){Zfrac.draw<-round(rtlnorm(1,Zfrac.Beta.in[2],Zfrac.Beta.in[3],0,1),2)}
      if(Zfrac.Beta.in[1]==4){Zfrac.draw<-round(runif(1,Zfrac.Beta.in[2],Zfrac.Beta.in[3]),2)}
      Input.draws[i,1]<-Zfrac.draw
    }
    if(Zfrac.Beta.in[4]>=0)
    {
      if(Zfrac.Beta.in[4]==10){Beta.draw<-round(rtnorm(1,Zfrac.Beta.in[5],Zfrac.Beta.in[6],0,10),2)}
      if(Zfrac.Beta.in[4]==30){Beta.draw<-round(rtlnorm(1,Zfrac.Beta.in[5],Zfrac.Beta.in[6],0,10),2)}
      if(Zfrac.Beta.in[4]==4){Beta.draw<-round(runif(1,Zfrac.Beta.in[5],Zfrac.Beta.in[6]),2)}
      Input.draws[i,7]<-Beta.draw
    }
    
    #Natural mortality
    if(M.in[1]>=0)
    {
      if(M.in[1]==0){M.draw<-round(rnorm(1,M.in[2],M.in[3]),2)}
      if(M.in[1]==3){M.draw<-round(rlnorm(1,log(M.in[2]),M.in[3]),2)}
      Input.draws[i,3]<-M.draw
    }
    #Growth parameters
    if(sum(L1.in[1:2])>0){L1.draw<-round(rnorm(1,L1.in[1],L1.in[1]),2); Input.draws[i,4]<-L1.draw}
    if(sum(Linf.in[1:2])>0){Linf.draw<-round(rnorm(1,Linf.in[1],Linf.in[1]),2); Input.draws[i,5]<-Linf.draw}
    if(sum(k.in[1:2])>0){k.draw<-round(rnorm(1,k.in[1],k.in[1]),2); Input.draws[i,6]<-k.draw}
    
    #Male draws
    if(genders==T)
    {
      if(M.in[4]>=0)
      {
        if(M.in[4]==0){M.draw.M<-round(rnorm(1,M.in[5],M.in[6]),2)}
        if(M.in[4]==3){M.draw.M<-round(rlnorm(1,log(M.in[5]),M.in[6]),2)}
        Input.draws.M[i,1]<-M.draw.M
      }
      if(sum(L1.in[3:4])>0){L1.draw.M<-round(rnorm(1,L1.in[3],L1.in[4]),2); Input.draws.M[i,2]<-L1.draw.M}
      if(sum(Linf.in[3:4])>0){Linf.draw.M<-round(rnorm(1,Linf.in[3],Linf.in[4]),2); Input.draws.M[i,3]<-Linf.draw.M}
      if(sum(k.in[3:4])>0){k.draw.M<-round(rnorm(1,k.in[3],k.in[4]),2); Input.draws.M[i,4]<-k.draw.M}
    }
    
    ### Change DAT and CTL inputs ###
    #change Depletion
    if(Dep.in[1]>=0)
    {
      dat.new<-readLines(paste(filepath,"/",file.name[1],".ss",sep=""))
      Dep.line<-strsplit(dat.new[grep("DEPLETION",dat.new)], " ")[[1]]
      fleet_num<-Dep.line[3]
      Dep.line[4]<-Dep.draw
      dat.new[grep("DEPLETION",dat.new)]<-paste(Dep.line,collapse=" ")
      write(dat.new,paste(filepath,"/",file.name[1],".ss",sep=""))
    }
    #change M
    if(M.in[1]>=0)
    {
      ctl.new<-readLines(paste(filepath,"/",file.name[2],".ss",sep=""))
      M.line<-strsplit(ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)], " ")[[1]]
      M.line[c(3,4)]<-M.draw
      ctl.new[grep("NatM_p_1_Fem_GP_1",ctl.new)]<-paste(M.line,collapse=" ")      
    }
    #change growth parameters
    if(sum(L1.in[1:2])>0)
    {
      L1.line<-strsplit(ctl.new[grep("L_at_Amin_Fem_GP_1",ctl.new)], " ")[[1]]
      L1.line[c(3,4)]<-L1.draw
      ctl.new[grep("L_at_Amin_Fem_GP_1",ctl.new)]<-paste(L1.line,collapse=" ")
    }
    if(sum(Linf.in[1:2])>0)
    {
      Linf.line<-strsplit(ctl.new[grep("L_at_Amax_Fem_GP_1",ctl.new)], " ")[[1]]
      Linf.line[c(3,4)]<-Linf.draw
      ctl.new[grep("L_at_Amax_Fem_GP_1",ctl.new)]<-paste(Linf.line,collapse=" ")
    }
    if(sum(k.in[1:2])>0)
    {
      k.line<-strsplit(ctl.new[grep("VonBert_K_Fem_GP_1",ctl.new)], " ")[[1]]
      k.line[c(3,4)]<-k.draw
      ctl.new[grep("VonBert_K_Fem_GP_1",ctl.new)]<-paste(k.line,collapse=" ")
    }
    #change male pararmeters
    if(genders==T)
    {
      #change M
      if(M.in[4]>=0)
      {
        M.line.M<-strsplit(ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)], " ")[[1]]
        M.line.M[c(3,4)]<-M.draw.M
        ctl.new[grep("NatM_p_1_Mal_GP_1",ctl.new)]<-paste(M.line.M,collapse=" ")        
      }
      #change growth parameters
      if(sum(L1.in[3:4])>0)
      {
        L1.line.M<-strsplit(ctl.new[grep("L_at_Amin_Mal_GP_1",ctl.new)], " ")[[1]]
        L1.line.M[c(3,4)]<-L1.draw.M
        ctl.new[grep("L_at_Amin_Mal_GP_1",ctl.new)]<-paste(L1.line.M,collapse=" ")
      }
      if(sum(Linf.in[3:4])>0)
      {
        Linf.line.M<-strsplit(ctl.new[grep("L_at_Amax_Mal_GP_1",ctl.new)], " ")[[1]]
        Linf.line.M[c(3,4)]<-Linf.draw.M
        ctl.new[grep("L_at_Amax_Mal_GP_1",ctl.new)]<-paste(Linf.line.M,collapse=" ")
      }
      if(sum(k.in[3:4])>0)
      {
        k.line.M<-strsplit(ctl.new[grep("VonBert_K_Mal_GP_1",ctl.new)], " ")[[1]]
        k.line.M[c(3,4)]<-k.draw.M
        ctl.new[grep("VonBert_K_Mal_GP_1",ctl.new)]<-paste(k.line.M,collapse=" ")
      }
    }
    #change R0
    if(i==1 & R_start[1]==1)
    {
      R0.draw<-round(rtlnorm(1,log(R_start[2]),0.5,3,15),2)
      R0.line<-strsplit(ctl.new[grep("SR_R0",ctl.new)], " ")[[1]]
      R0.line[c(3,4)]<-R0.draw
      ctl.new[grep("SR_R0",ctl.new)]<-paste(R0.line,collapse=" ")
    }
    if(i==1 & R_start[1]==0)
    {
      R0.draw<-R_start[2]
      R0.line<-strsplit(ctl.new[grep("SR_R0",ctl.new)], " ")[[1]]
      R0.line[c(3,4)]<-R0.draw
      ctl.new[grep("SR_R0",ctl.new)]<-paste(R0.line,collapse=" ")
    }
    #change h
    if(h.in[1]>=0)
    {
      h.line<-strsplit(ctl.new[grep("SR_steep",ctl.new)], " ")[[1]]
      h.line[c(3,4)]<-h.draw
      ctl.new[grep("SR_steep",ctl.new)]<-paste(h.line,collapse=" ")
    }
    #change Sfrac and Beta
    if(Zfrac.Beta.in[1]>=0)
    {
      Zfrac.line<-strsplit(ctl.new[grep("Zfrac",ctl.new)], " ")[[1]]
      Zfrac.line[c(3,4)]<-Zfrac.draw
      ctl.new[grep("Zfrac",ctl.new)]<-paste(Zfrac.line,collapse=" ")
    }
    if(Zfrac.Beta.in[4]>=0)
    {
      Beta.line<-strsplit(ctl.new[grep("Beta",ctl.new)], " ")[[1]]
      Beta.line[c(3,4)]<-Beta.draw
      ctl.new[grep("Beta",ctl.new)]<-paste(Beta.line,collapse=" ")
    }
    
    write(ctl.new,paste(filepath,"/",file.name[2],".ss",sep=""))
    
    #Run model
    RUN.SS(paste(filepath,"/",sep=""), ss.exe="SS3",ss.cmd=" -nohess -nox")
    
    #Evaluate convergence and record values
    rep.new<-readLines(paste(filepath,"/Report.sso",sep=""))
    while(abs(as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8]))==0)
    {
      R0.draw<-round(rtlnorm(1,log(R_start),0.5,1,30),2)
      ctl.new<-readLines(paste(filepath,"/",file.name[2],".ctl",sep=""))
      R0.line<-strsplit(ctl.new[grep("SR_R0",ctl.new)], " ")[[1]]
      R0.line[c(3,4)]<-R0.draw
      ctl.new[grep("SR_R0",ctl.new)]<-paste(R0.line,collapse=" ")
      write(ctl.new,paste(filepath,"/",file.name[2],".ctl",sep=""))
      RUN.SS(paste(filepath,"/",sep=""), ss.exe="SS3",ss.cmd=" -nohess -nox")
    }
    forecast.file<-readLines(paste(filepath,"/Forecast-report.sso",sep=""))
    for(iii in 1:length(year0:sb_ofl_yrs[1])){SB.out[i,iii]<-as.numeric(strsplit(rep.new[grep(paste("SPB_",sb.years[iii],sep=""),rep.new)], " ")[[1]][3])}
    Dep.series.out<-SB.out/SB.out[,1]
    Quant.out[i,1]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,2]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,3]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,4]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,5]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,6]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Mal_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,7]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Mal_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,8]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Mal_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,9]<-as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][3])
    Quant.out[i,10]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])
    Quant.out[i,11]<-as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[i,12]<-as.numeric(strsplit(rep.new[grep(paste("SPB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,13]<-as.numeric(strsplit(rep.new[grep(paste("SPB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[i,14]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,15]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,16]<-as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[i,17]<-as.numeric(strsplit(forecast.file[grep("calculate_FMSY",forecast.file)+13],split="[[:blank:]]+")[[1]][2])/as.numeric(strsplit(forecast.file[grep("BIO_Smry_unfished",forecast.file)],split="[[:blank:]]+")[[1]][2])
    Quant.out[i,18]<-as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
    Quant.out[i,19]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
    Quant.out[i,20]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)+2], " ")[[1]][2])
    Quant.out[i,21]<-as.numeric(strsplit(rep.new[grep("Convergence",rep.new)], " ")[[1]][2])
    if(Dep.in[1]>=0){Quant.out[i,22]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num," SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][7])}
    if(Dep.in[1]>=0){Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num," SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][8])}
    Quant.out[i,24]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8])
    Quant.out[i,25]<-as.numeric(strsplit(rep.new[grep(paste("F_",f_yr,sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
    Quant.out[i,26]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,27]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
    warnings.out<-readLines(paste(filepath,"/warning.sso",sep=""))
    num.warning<-as.numeric(strsplit(warnings.out[grep("N warnings",warnings.out)[1]], " ")[[1]][4])
    if(Dep.in[1]>=0){if(abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,21]<0.005&Quant.out[i,14]<100000){i<-i+1}}
    if(Dep.in[1]<0){if(Quant.out[i,21]<0.005&Quant.out[i,14]<100000){i<-i+1}}
    
    #Run using .par file if above conditions not met
    else
    {
      print("*** RUNNING WITH .PAR ***")
      starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
      sum_age_line<-strsplit(starter.new[grep("ss3.par",starter.new)], " ")[[1]]
      sum_age_line[1]<-1
      starter.new[grep("ss3.par",starter.new)]<-paste(sum_age_line, collapse=" ")
      write(starter.new,paste(filepath,"/starter.ss",sep=""))
      RUN.SS(paste(filepath,"/",sep=""), ss.exe="SS3",ss.cmd=" -nohess -nox")
      rep.new<-readLines(paste(filepath,"/Report.sso",sep=""))
      forecast.file<-readLines(paste(filepath,"/Forecast-report.sso",sep=""))
      for(iii in 1:length(year0:sb_ofl_yrs[1])){SB.out[i,iii]<-as.numeric(strsplit(rep.new[grep(paste("SPB_",sb.years[iii],sep=""),rep.new)], " ")[[1]][3])}
      Dep.series.out<-SB.out/SB.out[,1]
      Quant.out[i,1]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,2]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Fem_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,3]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Fem_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,4]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Fem_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,5]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,6]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Mal_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,7]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Mal_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,8]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Mal_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,9]<-as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][3])
      Quant.out[i,10]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])
      Quant.out[i,11]<-as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
      Quant.out[i,12]<-as.numeric(strsplit(rep.new[grep(paste("SPB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])
      Quant.out[i,13]<-as.numeric(strsplit(rep.new[grep(paste("SPB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
      Quant.out[i,14]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
      Quant.out[i,15]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
      Quant.out[i,16]<-as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SPB_Initial",rep.new)], " ")[[1]][3])
      Quant.out[i,17]<-as.numeric(strsplit(forecast.file[grep("calculate_FMSY",forecast.file)+13],split="[[:blank:]]+")[[1]][2])/as.numeric(strsplit(forecast.file[grep("BIO_Smry_unfished",forecast.file)],split="[[:blank:]]+")[[1]][2])
      Quant.out[i,18]<-as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
      Quant.out[i,19]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
      Quant.out[i,20]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)+2], " ")[[1]][2])
      Quant.out[i,21]<-as.numeric(strsplit(rep.new[grep("Convergence",rep.new)], " ")[[1]][2])
      if(Dep.in[1]>=0){Quant.out[i,22]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num," SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][7])}
      if(Dep.in[1]>=0){Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num," SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][8])}
      Quant.out[i,24]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8])
      Quant.out[i,25]<-as.numeric(strsplit(rep.new[grep(paste("F_",f_yr,sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
      Quant.out[i,26]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
      Quant.out[i,27]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
      starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
      sum_age_line<-strsplit(starter.new[grep("ss3.par",starter.new)], " ")[[1]]
      sum_age_line[1]<-0
      starter.new[grep("ss3.par",starter.new)]<-paste(sum_age_line, collapse=" ")
      write(starter.new,paste(filepath,"/starter.ss",sep=""))
      if(Dep.in[1]>=0){if(abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,21]<0.005&Quant.out[i,14]<100000){i<-i+1}}
      if(Dep.in[1]<0){if(Quant.out[i,21]<0.005&Quant.out[i,14]<100000){i<-i+1}}
    }
    ii<-ii+1
  }
  colnames(Quant.out)<-c("M_f","L1_f","Linf_f","k_f","M_m","L1_m","Linf_m","k_m","h","R0","SB0",paste("SPB_",sb_ofl_yrs[1],sep=""),"Dep",paste("OFL_",sb_ofl_yrs[2],sep=""),paste("AdjCatch_",sb_ofl_yrs[2],sep=""),"SBMSY/SB0","BMSY/B0","FMSY","-lnL","LL_survey","Gradient","Dep.Obs","Dep.Exp","R0_init",paste("F_",f_yr,"/FMSY",sep=""),paste("OFL_",sb_ofl_yrs[3],sep=""),paste("AdjCatch_",sb_ofl_yrs[3],sep=""))
  #if(all(is.na(Input.draws.M))=="FALSE"){)
  Input.draws<-cbind(Input.draws,Input.draws.M)
  if(ncol(Input.draws)==10){colnames(Input.draws)<-c("h","Dep","M_f","L1_f","Linf_f","k_f","M_m","L1_m","Linf_m","k_m")}
  if(ncol(Input.draws)==11){colnames(Input.draws)<-c("Sfrac","Dep","M_f","L1_f","Linf_f","k_f","Beta","M_m","L1_m","Linf_m","k_m")}
  end.time<-Sys.time()
  Spp.quant.out<-list(Input.draws,Quant.out,SB.out,Dep.series.out,ii-1,(as.numeric(end.time)-as.numeric(start.time))/60)
  names(Spp.quant.out)<-c("Priors","Posteriors","SB_series","Depeltion_series","Total draws","Runtime_minutes")
  return(Spp.quant.out)
}

loadnnames <- function (file, ...) 
{
  ls.ext <- function(file) {
    local({
      base::load(file)
      base::ls()
    })
  }
  base::load(file, .GlobalEnv, ...)
  ls.ext(file)
}

aurora.BC.out<-SS_output(getwd())
SS_plots(aurora.BC.out)

ARRA.SSS<-SSS.jr.exSSS(getwd(),file.name=c("ARRA_dat3","ARRA_ctl5"),reps=1000,M.in=c(3,0.0405,exp(-3.353),3,0.0429,exp(-3.295)),Dep.in=c(10,0.3,0.2),h.in=c(2,0.779,0.152),L1.in=c(0,0,0,0),Linf.in=c(0,0,0,0),k.in=c(0,0,0,0),Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),R_start=6,sum_age=0,sb_ofl_yrs=c(2013,2015,2016),f_yr=2012,year0=1916,genders=T)
save(ARRA.SSS,file=paste(getwd(),"/SSS.DMP",sep=""))

#Get SSS results for table
obj.name<-loadnnames("C:/Users/Jason.Cope/Documents/Stock Assessments/Aurora/SSS/SSS.DMP")  
obj.out<-get(obj.name)
medians.sss<-c(median(obj.out$Post$SB0),median(obj.out$Post$SPB_2013),median(obj.out$Post$Dep),median(obj.out$Post$'F_2012/FMSY'),median(obj.out$Post$OFL_2015),median(obj.out$Post$OFL_2016))
means.sss<-c(mean(obj.out$Post$SB0),mean(obj.out$Post$SPB_2013),mean(obj.out$Post$Dep),mean(obj.out$Post$'F_2012/FMSY'),mean(obj.out$Post$OFL_2015),mean(obj.out$Post$OFL_2016))
CI95hi<-c(quantile(obj.out$Post$SB0,0.95),quantile(obj.out$Post$SPB_2013,0.95),quantile(obj.out$Post$Dep,0.95),quantile(obj.out$Post$'F_2012/FMSY',0.95),quantile(obj.out$Post$OFL_2015,0.95),quantile(obj.out$Post$OFL_2016,0.95))
CI95lo<-c(quantile(obj.out$Post$SB0,0.05),quantile(obj.out$Post$SPB_2013,0.05),quantile(obj.out$Post$Dep,0.05),quantile(obj.out$Post$'F_2012/FMSY',0.05),quantile(obj.out$Post$OFL_2015,0.05),quantile(obj.out$Post$OFL_2016,0.05))
SSS.outputs.datamod<-rbind(medians.sss,neans.sss,CI95hi,CI95lo)
#save(SSS.outputs.datamod,file="C:/Users/Jason.Cope/Documents/Stock Assessments/Aurora/SSS/SSS_out.DMP")
#load(paste(getwd(),"/SSS_out.DMP",sep=""))
aurora.BC.out
ht <- 6.5; wd<- 6.5
png("C:/Users/Jason.Cope/Documents/Stock Assessments/Aurora/SSS//SSSvSS.png",height=ht,width=wd,pointsize=10,units="in",res=300)
par(mfcol=c(2,1))
  #Spawning biomass
  years<-c(1:length(aurora.BC.out$timeseries$Yr))
  SB.SS<-aurora.BC.out$timeseries$Spawn[years[aurora.BC.out$timeseries$Yr==aurora.BC.out$startyr]:years[aurora.BC.out$timeseries$Yr==aurora.BC.out$endyr+1]]
  SB.SS.hi<-SB.SS+1.96*aurora.BC.out$derived_quants[grep(aurora.BC.out$startyr,aurora.BC.out$derived_quants$LABEL)[1]:grep(aurora.BC.out$endyr+1,aurora.BC.out$derived_quants$LABEL)[1],3]
  SB.SS.lo<-SB.SS-1.96*aurora.BC.out$derived_quants[grep(aurora.BC.out$startyr,aurora.BC.out$derived_quants$LABEL)[1]:grep(aurora.BC.out$endyr+1,aurora.BC.out$derived_quants$LABEL)[1],3]
  SB.SSS<-apply(obj.out$SB_series,2,quantile,probs=0.5)
  SB.SSS.hi<-apply(obj.out$SB_series,2,quantile,probs=0.975)
  SB.SSS.lo<-apply(obj.out$SB_series,2,quantile,probs=0.025)
  col.SS<-col2rgb("gray")
  col.SSS<-col2rgb("red")
  xx <- c(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1)), rev(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1))))
  yy.sb.SS <- c(SB.SS.lo, rev(SB.SS.hi))
  yy.sb.SSS <- c(SB.SSS.lo, rev(SB.SSS.hi))
  plot(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1)),SB.SS,type="l",lwd=2,xlab="Year",ylab="Spawning biomass",ylim=c(0,max(SB.SS.hi,SB.SSS.hi)))
  polygon(xx,yy.sb.SS, col=rgb(col.SS[1]/255,col.SS[2]/255,col.SS[3]/255,alpha=1), border = "black")
  polygon(xx,yy.sb.SSS, col=rgb(col.SSS[1]/255,col.SSS[2]/255,col.SSS[3]/255,alpha=0.3), border = "red")
  lines(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1)),SB.SS,lwd=3,lty=1,col="black")
  lines(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1)),SB.SSS,lwd=3,lty=2,col="red")
  
  #Depletion
  Dep.SS<-aurora.BC.out$timeseries$Spawn[years[aurora.BC.out$timeseries$Yr==aurora.BC.out$startyr]:years[aurora.BC.out$timeseries$Yr==aurora.BC.out$endyr+1]]/aurora.BC.out$timeseries$Spawn[1]
  Dep.SS.hi<-Dep.SS+1.96*c(0,aurora.BC.out$derived_quants[grep(aurora.BC.out$startyr+1,aurora.BC.out$derived_quants$LABEL)[5]:grep(aurora.BC.out$endyr+1,aurora.BC.out$derived_quants$LABEL)[5],3])
  Dep.SS.lo<-Dep.SS-1.96*c(0,aurora.BC.out$derived_quants[grep(aurora.BC.out$startyr+1,aurora.BC.out$derived_quants$LABEL)[5]:grep(aurora.BC.out$endyr+1,aurora.BC.out$derived_quants$LABEL)[5],3])
  Dep.SSS<-apply(obj.out$Depeltion_series,2,quantile,probs=0.5)
  Dep.SSS.hi<-apply(obj.out$Depeltion_series,2,quantile,probs=0.975)
  Dep.SSS.lo<-apply(obj.out$Depeltion_series,2,quantile,probs=0.025)
  yy.dep.SS <- c(Dep.SS.lo, rev(Dep.SS.hi))
  yy.dep.SSS <- c(Dep.SSS.lo, rev(Dep.SSS.hi))
  plot(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1)),Dep.SS,type="l",lwd=2,xlab="Year",ylab="Depletion",ylim=c(0,1.2))
  polygon(xx,yy.dep.SS, col=rgb(col.SS[1]/255,col.SS[2]/255,col.SS[3]/255,alpha=1), border = "black")
  polygon(xx,yy.dep.SSS, col=rgb(col.SSS[1]/255,col.SSS[2]/255,col.SSS[3]/255,alpha=0.3), border = "red")
  lines(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1)),Dep.SS,lwd=3,lty=1,col="black")
  lines(c(aurora.BC.out$startyr:(aurora.BC.out$endyr+1)),Dep.SSS,lwd=3,lty=2,col="red")
  abline(h=c(0.4,0.25),lty=c(2,1))
  dev.off()


##########################################
### Plot DB-SRA catches to new catches ###
##########################################
ARRA.tot.catches<-read.table('clipboard',header=T)
plot(ARRA.tot.catches$Year,ARRA.tot.catches$SSS,type="l",col="blue",lwd=2,xlab="Year",ylab="Removals (mt)",ylim=c(0,round(max(ARRA.tot.catches[,2:3],na.rm=T))))
lines(ARRA.tot.catches$Year,ARRA.tot.catches$SS,lwd=2,col="black")
smartlegend("left","top",c("DB-SRA 2010","SS 2013"),col=c("blue","black"),lwd=2,bty="n")
