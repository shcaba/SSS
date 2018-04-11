##' This file contails all the SSS code and associated functions
##' @author Jason Cope and Chantel Wetzel
##' @export
##' @import msm
##' @import EnvStats
##' @import r4ss




SSS<-function(filepath,file.name,reps=1000,seed.in=19,Dep.in=c(1,0.4,0.1),M.in=c(3,0.1,0.4,3,0.1,0.4),SR_type=3,h.in=c(1,0.6,0.2),FMSY_M.in=c(1,0.5,0.1),
  BMSY_B0.in=c(1,0.5,0.1),L1.in=c(0,0,0,0),Linf.in=c(0,0,0,0),k.in=c(0,0,0,0),Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),R_start=c(0,8),doR0.loop=c(1,4.1,12.1,0.5),
  sum_age=0,sb_ofl_yrs=c(2010,2011,2012),f_yr=2012,year0=1916,genders=F,BH_FMSY_comp=F)
{

  require(msm)
  require(EnvStats)

  set.seed(seed.in)
  start.time<-Sys.time()
  Spp.quant.out<-list()
  Input.draws<-as.data.frame(matrix(NA,nrow=reps,ncol=10))
  if(SR_type==7){Input.draws<-as.data.frame(matrix(NA,nrow=reps,ncol=11))}
  Input.draws.M<-as.data.frame(matrix(NA,nrow=reps,ncol=4))
  if(SR_type>=8 | h.in[1]==99 | BH_FMSY_comp==T)
  {
    SR_expo.out<-SRMQ_par<-rep(NA,reps)
    names(SR_expo.out)<-"gRicker_Beta"
    Input.draws.MQs<-as.data.frame(matrix(NA,nrow=reps,ncol=3))
    colnames(Input.draws.MQs)<-c("FMSY/M","BMSY/B0","FMSY")
  }
  sb.years<-c(year0:sb_ofl_yrs[1])
  Quant.out<-as.data.frame(matrix(NA,nrow=reps,ncol=28))
  Quant.out.bad<-as.data.frame(matrix(NA,1,ncol=28))
  SB.out<-TB.out<-SumAge.out<-SPR.out<-as.data.frame(matrix(NA,nrow=reps,ncol=length(sb.years)))
  colnames(SB.out)<-colnames(TB.out)<-colnames(SumAge.out)<-colnames(SPR.out)<-sb.years
  
  starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
  sum_age_line<-strsplit(starter.new[grep("summary biomass",starter.new)], " ")[[1]]
  sum_age_line[1]<-sum_age
  starter.new[grep("summary biomass",starter.new)]<-paste(sum_age_line, collapse=" ")
  write(starter.new,paste(filepath,"/starter.ss",sep=""))
  i<-1
  ii<-1
  xxx<-1
  n.bad<-0
  while(i<reps+1)
  {
    print(paste0("Attempt ",ii))
    if(ii%%10==0){print(file.name[1])}
    if(file.exists(paste(filepath,"/Report.sso",sep=""))){file.remove(paste(filepath,"/Report.sso",sep=""))}
    ### Draw random values ###
    #- values = no draws taken; one value fixed in the model
    #0 = normal
    #10 = truncated normal
    #1 = symmetric beta
    #2 = full beta (rbeta)
    #3 = lognormal
    #30 = truncated lognormal
    #4 = uniform
    #99 = used only for the steepness parameter. Indicates h will come from FMSY/M prior
    dat.new<-readLines(paste(filepath,"/",file.name[1],sep=""))
    ctl.new<-readLines(paste(filepath,"/",file.name[2],sep=""))
    
    if(length(Dep.in)==3)
    {
      if(Dep.in[1]==2){Dep.draw<-round(1-rbeta.ab(1,1-Dep.in[2],Dep.in[3],0.05,0.95),2)}
      if(Dep.in[1]==4){Dep.draw<-round(runif(1,Dep.in[2],Dep.in[3]),2)}
      if(Dep.in[1]==10){Dep.draw<-round(rtnorm(1,Dep.in[2],Dep.in[3],0.01,1),2)}
      Input.draws[i,2]<-Dep.draw
    }
    if(length(Dep.in)>3)
    {
      Input.draws[i,2]<-Dep.draw<-Dep.in[i]    
    }
    
    
    #Natural mortality
    if(FMSY_M.in[1]<0 & BMSY_B0.in[1]<0)
    {
      if(M.in[1]>=0 & length(M.in)==6)
        {
          if(M.in[1]==0){M.draw<-round(rnorm(1,M.in[2],M.in[3]),2)}
          if(M.in[1]==3){M.draw<-round(rlnorm(1,log(M.in[2]),M.in[3]),2)}
          if(M.in[1]==4){M.draw<-round(runif(1,M.in[2],M.in[3]),2)}
          Input.draws[i,3]<-M.draw
        }
      else{Input.draws[i,3]<-M.draw<-M.in[i,1]}

      if(genders==T)
      {
        if(length(M.in)==6)
        {
          if(M.in[4]==0){M.draw.M<-round(rnorm(1,M.in[5],M.in[6]),2)}
          if(M.in[4]==3){M.draw.M<-round(rlnorm(1,log(M.in[5]),M.in[6]),2)}
          if(M.in[4]==4){M.draw.M<-round(runif(1,M.in[5],M.in[6]),2)}
          Input.draws.M[i,1]<-M.draw.M
          Input.draws[i,7]<-M.draw.M
        }
        if(length(M.in)>6){Input.draws[i,7]<-M.draw.M<-M.in[i,2]}
      }
    }

    #Growth parameters
    if(sum(L1.in[1:2])>0){L1.draw<-round(rnorm(1,L1.in[1],L1.in[1]),2); Input.draws[i,4]<-L1.draw}
    if(sum(Linf.in[1:2])>0){Linf.draw<-round(rnorm(1,Linf.in[1],Linf.in[1]),2); Input.draws[i,5]<-Linf.draw}
    if(sum(k.in[1:2])>0){k.draw<-round(rnorm(1,k.in[1],k.in[1]),2); Input.draws[i,6]<-k.draw}
    #Male draws
    if(genders==T)
    {
      if(sum(L1.in[3:4])>0){L1.draw.M<-round(rnorm(1,L1.in[3],L1.in[4]),2); Input.draws.M[i,2]<-L1.draw.M}
      if(sum(Linf.in[3:4])>0){Linf.draw.M<-round(rnorm(1,Linf.in[3],Linf.in[4]),2); Input.draws.M[i,3]<-Linf.draw.M}
      if(sum(k.in[3:4])>0){k.draw.M<-round(rnorm(1,k.in[3],k.in[4]),2); Input.draws.M[i,4]<-k.draw.M}
    }
    
    #Steenpess
    if(h.in[1]>=0 & length(h.in)==3 & FMSY_M.in[1]<0 & BMSY_B0.in[1]<0)
    {
      if(h.in[1]==2){h.draw<-round(rbeta.ab(1,h.in[2],h.in[3],0.21,0.99),2)}
      if(h.in[1]==10){h.draw<-round(rtnorm(1,h.in[2],h.in[3],0.21,0.99),2)}
      if(h.in[1]==30){h.draw<-round(rlnormTrunc(1,log(h.in[2]),h.in[3],0.21,0.99),2)}
      if(h.in[1]==4){h.draw<-round(runif(1,h.in[2],h.in[3]),2)}
      Input.draws[i,1]<-h.draw
    }
    if(length(h.in)>3 & FMSY_M.in[1]<0 & BMSY_B0.in[1]<0){Input.draws[i,1]<-h.draw<-h.in[i]}
    #Get steepness for BH from FMSY
    if(h.in[1]==99 | BH_FMSY_comp==T & FMSY_M.in[1]>=0)
        {
      SRdat.new<-readLines(paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
      x<-1
      while(x==1)
      {
      #Draw Ms
        if(M.in[1]>=0 & length(M.in)==6)
        {
          if(M.in[1]==0){M.draw<-round(rnorm(1,M.in[2],M.in[3]),2)}
          if(M.in[1]==3){M.draw<-round(rlnorm(1,log(M.in[2]),M.in[3]),2)}
          if(M.in[1]==4){M.draw<-round(rlnorm(1,log(M.in[2]),M.in[3]),2)}
          Input.draws[i,3]<-M.draw
        }
        else{Input.draws[i,3]<-M.draw<-M.in[i,1]}
        #Male draws
        if(genders==T)
          {
            if(length(M.in)==6)
            {
              if(M.in[4]==0){M.draw.M<-round(rnorm(1,M.in[5],M.in[6]),2)}
              if(M.in[4]==3){M.draw.M<-round(rlnorm(1,log(M.in[5]),M.in[6]),2)}
              if(M.in[4]==4){M.draw.M<-round(runif(1,M.in[5],M.in[6]),2)}
              Input.draws.M[i,1]<-M.draw.M
              Input.draws[i,7]<-M.draw.M
            }
            if(length(M.in)>6){Input.draws[i,7]<-M.draw.M<-M.in[i,2]}
           }
        
        if(FMSY_M.in[1]==2){FMSY_M.draw<-round(rbeta.ab(1,FMSY_M.in[2],FMSY_M.in[3],0.0001,10),2)}
        if(FMSY_M.in[1]==10){FMSY_M.draw<-round(rtnorm(1,FMSY_M.in[2],FMSY_M.in[3],0.0001,10),2)}
        if(FMSY_M.in[1]==30){FMSY_M.draw<-round(rlnormTrunc(1,log(FMSY_M.in[2]),FMSY_M.in[3],0.0001,10),2)}
        if(FMSY_M.in[1]==4){FMSY_M.draw<-round(runif(1,FMSY_M.in[2],FMSY_M.in[3]),2)}
        if(length(FMSY_M.in)>3){FMSY_M.draw<-FMSY_M.in[i]}
        Input.draws.MQs[i,1]<-FMSY_M.draw

        if(genders==F)
          {
            SRdat.new[5]<-paste(as.character(M.draw),as.character(M.draw),sep=" ")
            M_FMSY<-M.draw
          }
        if(genders==T)
          {
            SRdat.new[5]<-paste(as.character(M.draw),as.character(M.draw.M),sep=" ")
            M_FMSY<-mean(c(M.draw,M.draw.M))
          }
        SRdat.new[17]<-as.character(FMSY_M.draw*M_FMSY)
        Input.draws.MQs[i,3]<-FMSY_M.draw*M_FMSY
        SRdat.new[21]<-"-0.51"
        write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
        RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
        if(file.exists(paste(filepath,"/h_BMSY_FMSY/bmsyb0.rep",sep=""))==TRUE)
        {
        
        SRMQ_rep.new<-readLines(paste(filepath,"/h_BMSY_FMSY/bmsyb0.rep",sep=""))
        SRMQ_par[i]<-as.numeric(strsplit(readLines(paste(filepath,"/h_BMSY_FMSY/bmsyb0.par",sep="")),split=" ")[[1]][11])
        steep_expo.out<-as.numeric(strsplit(SRMQ_rep.new[1]," ")[[1]])
        if(steep_expo.out[1]<0.25|steep_expo.out[1]>=1|steep_expo.out[2]>10|SRMQ_par[i]>0.1)
          {
            SRdat.new[21]<-"-0.1"
            write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
            RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
            if(steep_expo.out[1]<0.25|steep_expo.out[1]>=1|steep_expo.out[2]>10|SRMQ_par[i]>0.1)
            {
            SRdat.new[21]<-"-1.51"
            write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
            RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
            }
            if(steep_expo.out[1]<0.25|steep_expo.out[1]>=1|steep_expo.out[2]>10|SRMQ_par[i]>0.1)
            {
            SRdat.new[21]<-"-0.51"
            SRdat.new[23]<-"-0.31"
            write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
            RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
            }
          }
        if(steep_expo.out[1]>=0.25&steep_expo.out[1]<1&steep_expo.out[2]<=10&SRMQ_par[i]<0.1){x<-2}          
        }
      }
      Input.draws[i,1]<-h.draw<-steep_expo.out[1]
      SR_expo.out[i]<-steep_expo.out[2]
      file.remove(paste(filepath,"/h_BMSY_FMSY/bmsyb0.rep",sep=""))
      file.remove(paste(filepath,"/h_BMSY_FMSY/bmsyb0.par",sep=""))
    }
    

    if(FMSY_M.in[1]>=0 & BMSY_B0.in[1]>=0)
    {
      SRdat.new<-readLines(paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
      #SRdat.new.mod<-strsplit(SRdat.new[grep("Amax",SRdat.new)+1][1]," ")[[1]]
      #Amax.line<-strsplit(dat.new[grep("Number of ages",dat.new)], " ")[[1]]
      #k.line<-strsplit(ctl.new[grep("VonBert_K_Fem_GP_1",ctl.new)]," ")[[1]]
      #ltwt.exp.line<-strsplit(ctl.new[grep("Wtlen_2_Fem",ctl.new)]," ")[[1]]
      x<-1
      while(x==1)
      {
      #Draw Ms
        if(M.in[1]>=0 & length(M.in)==6)
        {
          if(M.in[1]==0){M.draw<-round(rnorm(1,M.in[2],M.in[3]),2)}
          if(M.in[1]==3){M.draw<-round(rlnorm(1,log(M.in[2]),M.in[3]),2)}
          if(M.in[1]==4){M.draw<-round(runif(1,M.in[2],M.in[3]),2)}
          Input.draws[i,3]<-M.draw
        }
        else{Input.draws[i,3]<-M.draw<-M.in[i,1]}
        #Male draws
        if(genders==T)
          {
            if(length(M.in)==6)
            {
              if(M.in[4]==0){M.draw.M<-round(rnorm(1,M.in[5],M.in[6]),2)}
              if(M.in[4]==3){M.draw.M<-round(rlnorm(1,log(M.in[5]),M.in[6]),2)}
              if(M.in[4]==4){M.draw.M<-round(runif(1,M.in[5],M.in[6]),2)}
              Input.draws.M[i,1]<-M.draw.M
              Input.draws[i,7]<-M.draw.M
            }
            if(length(M.in)>6){Input.draws[i,7]<-M.draw.M<-M.in[i,2]}
           }
        
        if(FMSY_M.in[1]==2){FMSY_M.draw<-round(rbeta.ab(1,FMSY_M.in[2],FMSY_M.in[3],0.0001,10),2)}
        if(FMSY_M.in[1]==10){FMSY_M.draw<-round(rtnorm(1,FMSY_M.in[2],FMSY_M.in[3],0.0001,10),2)}
        if(FMSY_M.in[1]==30){FMSY_M.draw<-round(rlnormTrunc(1,log(FMSY_M.in[2]),FMSY_M.in[3],0.0001,10),2)}
        if(FMSY_M.in[1]==4){FMSY_M.draw<-round(runif(1,FMSY_M.in[2],FMSY_M.in[3]),2)}
        if(length(FMSY_M.in)>3){FMSY_M.draw<-FMSY_M.in[i]}
        if(BMSY_B0.in[1]==2){BMSY_B0.draw<-round(rbeta.ab(1,BMSY_B0.in[2],BMSY_B0.in[3],0.01,0.99),2)*100}
        if(BMSY_B0.in[1]==10){BMSY_B0.draw<-round(rtnorm(1,BMSY_B0.in[2],BMSY_B0.in[3],0.01,0.99),2)*100}
        if(BMSY_B0.in[1]==30){BMSY_B0.draw<-round(rlnormTrunc(1,log(BMSY_B0.in[2]),BMSY_B0.in[3],0.01,0.99),2)*100}
        if(BMSY_B0.in[1]==4){BMSY_B0.draw<-round(runif(1,BMSY_B0.in[2],BMSY_B0.in[3]),2)*100}
        if(length(BMSY_B0.in)>3){BMSY_B0.draw<-BMSY_B0.in[i]}
        Input.draws.MQs[i,1]<-FMSY_M.draw
        Input.draws.MQs[i,2]<-BMSY_B0.draw
        #SRdat.new.mod[1]<-as.numeric(Amax.line[1])
        if(genders==F)
          {
            SRdat.new[5]<-paste(as.character(M.draw),as.character(M.draw),sep=" ")
            M_FMSY<-M.draw
          }
        if(genders==T)
          {
            SRdat.new[5]<-paste(as.character(M.draw),as.character(M.draw.M),sep=" ")
            M_FMSY<-mean(c(M.draw,M.draw.M))
          }
        #print(c(M.draw,M.draw.M,M_FMSY,FMSY_M.draw,BMSY_B0.draw))
        #SRdat.new.mod[3]<-Winf
        #SRdat.new.mod[4]<-as.numeric(k.line[3])
        #SRdat.new.mod[5]<-as.numeric(ltwt.exp.line[3])
        #SRdat.new.mod[6]<-SRdat.new.mod[7]<-Amat
        #SRdat.new.mod[8]<-SR_type-7 #-7 converts from SS SR_type to MQ SR_type
        SRdat.new[17]<-as.character(FMSY_M.draw*M_FMSY)
        SRdat.new[19]<-as.character(BMSY_B0.draw)
        Input.draws.MQs[i,3]<-FMSY_M.draw*M_FMSY
        #SRdat.new[grep("Amax",SRdat.new)+1][1]<-paste(SRdat.new.mod,collapse=" ")
        SRdat.new[21]<-"-0.51"
        write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
        RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
        if(file.exists(paste(filepath,"/h_BMSY_FMSY/bmsyb0.rep",sep=""))==TRUE)
        {
        
        SRMQ_rep.new<-readLines(paste(filepath,"/h_BMSY_FMSY/bmsyb0.rep",sep=""))
        SRMQ_par[i]<-as.numeric(strsplit(readLines(paste(filepath,"/h_BMSY_FMSY/bmsyb0.par",sep="")),split=" ")[[1]][11])
        steep_expo.out<-as.numeric(strsplit(SRMQ_rep.new[1]," ")[[1]])
        if(steep_expo.out[1]<0.25|steep_expo.out[1]>=1|steep_expo.out[2]>10|SRMQ_par[i]>0.1)
          {
            #print("retry 1")
            SRdat.new[21]<-"-0.1"
            write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
            RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
            if(steep_expo.out[1]<0.25|steep_expo.out[1]>=1|steep_expo.out[2]>10|SRMQ_par[i]>0.1)
            {
            #print("retry 2")
            SRdat.new[21]<-"-1.51"
            write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
            RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
            }
            if(steep_expo.out[1]<0.25|steep_expo.out[1]>=1|steep_expo.out[2]>10|SRMQ_par[i]>0.1)
            {
            #print("retry 3")
            SRdat.new[21]<-"-0.51"
            SRdat.new[23]<-"-0.31"
            write(SRdat.new,paste(filepath,"/h_BMSY_FMSY/BMSYB0.dat",sep=""))
            RUN.SS(paste(filepath,"/h_BMSY_FMSY/",sep=""), ss.exe="BMSYB0",ss.cmd=" -nox > out.txt")
            }
          }
        #print(c(steep_expo.out[1],steep_expo.out[2],SRMQ_par[i]>0.1))
        #if(steep_expo.out[1]<0.25|steep_expo.out[1]>=1|steep_expo.out[2]>10|SRMQ_par[i]>0.1){x<-1}
        if(steep_expo.out[1]>=0.25&steep_expo.out[1]<1&steep_expo.out[2]<=10&SRMQ_par[i]<0.1){x<-2}          
        }
      }
      #print(SRMQ_par[i])
      Input.draws[i,1]<-h.draw<-steep_expo.out[1]
      SR_expo.out[i]<-steep_expo.out[2]
      file.remove(paste(filepath,"/h_BMSY_FMSY/bmsyb0.rep",sep=""))
      file.remove(paste(filepath,"/h_BMSY_FMSY/bmsyb0.par",sep=""))
    }
    
    #Three parameter S-R model
    if(Zfrac.Beta.in[1]>=0)
    {
      if(Zfrac.Beta.in[1]==10){Zfrac.draw<-round(rtnorm(1,Zfrac.Beta.in[2],Zfrac.Beta.in[3],0,1),2)}
      if(Zfrac.Beta.in[1]==30){Zfrac.draw<-round(rlnormTrunc(1,log(Zfrac.Beta.in[2]),Zfrac.Beta.in[3],0,1),2)}
      if(Zfrac.Beta.in[1]==4){Zfrac.draw<-round(runif(1,Zfrac.Beta.in[2],Zfrac.Beta.in[3]),2)}
      Input.draws[i,1]<-Zfrac.draw
    }
    if(Zfrac.Beta.in[4]>=0)
    {
      if(Zfrac.Beta.in[4]==10){Beta.draw<-round(rtnorm(1,Zfrac.Beta.in[5],Zfrac.Beta.in[6],0,10),2)}
      if(Zfrac.Beta.in[4]==30){Beta.draw<-round(rlnormTrunc(1,log(Zfrac.Beta.in[5]),Zfrac.Beta.in[6],0,10),2)}
      if(Zfrac.Beta.in[4]==4){Beta.draw<-round(runif(1,Zfrac.Beta.in[5],Zfrac.Beta.in[6]),2)}
      Input.draws[i,7]<-Beta.draw
    }
    
    ### Change DAT and CTL inputs ###
    #change Depletion
    if(Dep.in[1]>=0)
    {
      Dep.line<-strsplit(dat.new[grep("DEPLETION",dat.new)], " ")[[1]]
      fleet_num<-Dep.line[3]
      Dep.line[4]<-Dep.draw
      dat.new[grep("DEPLETION",dat.new)]<-paste(Dep.line,collapse=" ")
      write(dat.new,paste(filepath,"/",file.name[1],sep=""))
    }
    #change M
    Sys.sleep(1)
    if(M.in[1]>=0)
    {
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
    
    SRtype.line<-strsplit(ctl.new[grep("SR_function",ctl.new)], " ")[[1]]
    SRtype.line[1]<-SR_type
    ctl.new[grep("SR_function",ctl.new)]<-paste(SRtype.line,collapse=" ")
    #change R0
    if(R_start[1]==1)
    {
      R0.draw<-round(rlnormTrunc(1,log(R_start[2]),0.5,3,15),2)
      R0.line<-strsplit(ctl.new[grep("R0",ctl.new)], " ")[[1]]
      R0.line[c(3,4)]<-R0.draw
      ctl.new[grep("R0",ctl.new)]<-paste(R0.line,collapse=" ")
    }
    if(R_start[1]==0)
    {
      R0.draw<-R_start[2]
      R0.line<-strsplit(ctl.new[grep("R0",ctl.new)], " ")[[1]]
      R0.line[c(3,4)]<-R0.draw
      ctl.new[grep("R0",ctl.new)]<-paste(R0.line,collapse=" ")
    }
    

    #change h
    if(h.in[1]>=0)
    {
      h.line<-strsplit(ctl.new[grep("SR_steep",ctl.new)], " ")[[1]]
      h.line[c(3,4)]<-h.draw
      ctl.new[grep("SR_steep",ctl.new)]<-paste(h.line,collapse=" ")
    }
    
    if(SR_type>=8)
    {
      h.line<-strsplit(ctl.new[grep("SR_steep",ctl.new)], " ")[[1]]
      h.line[c(3,4)]<-h.draw
      ctl.new[grep("SR_steep",ctl.new)]<-paste(h.line,collapse=" ")
      SRexpo.line<-strsplit(ctl.new[grep("SR_Beta",ctl.new)], " ")[[1]]
      SRexpo.line[c(3,4)]<-SR_expo.out[i]
      ctl.new[grep("SR_Beta",ctl.new)]<-paste(SRexpo.line,collapse=" ")
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
    
    write(ctl.new,paste(filepath,"/",file.name[2],sep=""))
    
    #Run model
    RUN.SS(paste(filepath,"/",sep=""), ss.exe="ss",ss.cmd=" -nohess -nox > out.txt 2>&1")
    
    #Evaluate convergence and record values
    rep.new<-readLines(paste(filepath,"/Report.sso",sep=""))
    Sys.sleep(0.5)
    #The R0 loop. Makes sure R0 is changing from input and that the depletion match is happening. Ensures good models are kept.
    #Profiles over R0 when intial run fails.
    if(doR0.loop[1]>0)
    {
        xx<-1
        #if(doR0.loop[1]==1){R0.explore<-seq(doR0.loop[2],doR0.loop[3],0.5)}
        #if(doR0.loop[1]==2){R0.explore<-seq(doR0.loop[2],doR0.loop[3],-0.5)}
        R0.explore<-seq(doR0.loop[2],doR0.loop[3],doR0.loop[4])
        if(is.na(as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3]))==TRUE){Dep.out.testR0<-10}
        if(is.na(as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3]))==FALSE){Dep.out.testR0<-as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3])}
      #  while(abs(as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8]))==0)
      while(abs(Dep.out.testR0-Dep.draw)>0.01)
      {
        print(paste0("IN THE LOOP: R0= ",R0.explore[xx]))
        ctl.new<-readLines(paste(filepath,"/",file.name[2],sep=""))
        R0.line<-strsplit(ctl.new[grep("R0",ctl.new)], " ")[[1]]
        #R0.draw<-round(rlnormTrunc(1,log(R_start[2]),0.5,4,12),2)
        #R0.explore<-seq(5,11,0.5)
        R0.line[c(3,4)]<-R0.explore[xx]
        #print(R0.explore[xx])
        ctl.new[grep("R0",ctl.new)]<-paste(R0.line,collapse=" ")
        write(ctl.new,paste(filepath,"/",file.name[2],sep=""))
        RUN.SS(paste(filepath,"/",sep=""), ss.exe="ss",ss.cmd=" -nohess -nox > out.txt 2>&1")
        rep.new<-readLines(paste(filepath,"/Report.sso",sep=""))
        if(is.na(as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3]))==TRUE){Dep.out.testR0<-10}
        if(is.na(as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3]))==FALSE){Dep.out.testR0<-as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3])}
        #print(abs(as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8])))
        print(paste0("Is depletion difference < 0.01? ",(abs(as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3])-Dep.draw)<0.01)))
        xx<-xx+1
        if(xx==length(R0.explore)){break}
      }
    }
    #print(FMSY_M.in[i])
    
    #Begin extracting info from run
    forecast.file<-readLines(paste(filepath,"/Forecast-report.sso",sep=""))
    for(iii in 1:length(year0:sb_ofl_yrs[1]))
    	{
    		SB.out[i,iii]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb.years[iii],sep=""),rep.new)], " ")[[1]][3])
    		TB.out[i,iii]<-as.numeric(strsplit(rep.new[grep("TIME_SERIES",rep.new)+3+iii], " ")[[2]][5])
    		SumAge.out[i,iii]<-as.numeric(strsplit(rep.new[grep("TIME_SERIES",rep.new)+3+iii], " ")[[2]][6])
    		SPR.out[i,iii]<-as.numeric(strsplit(rep.new[grep("SPR_series",rep.new)+5+iii], " ")[[2]][8])
     	}
    		
    Dep.series.out<-SB.out/SB.out[,1]
    Quant.out[i,1]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,2]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,3]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Fem_GP_1",rep.new)], " ")[[1]][3])
    Quant.out[i,4]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Fem_GP_1",rep.new)], " ")[[1]][3])
    if(genders==T)
    {
      Quant.out[i,5]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,6]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Mal_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,7]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Mal_GP_1",rep.new)], " ")[[1]][3])
      Quant.out[i,8]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Mal_GP_1",rep.new)], " ")[[1]][3])
    }
    Quant.out[i,9]<-as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][3])
    Quant.out[i,10]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])
    Quant.out[i,11]<-as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[i,12]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,13]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[i,14]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,15]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,16]<-as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
    Quant.out[i,17]<-as.numeric(strsplit(forecast.file[grep("calculate_FMSY",forecast.file)+13],split="[[:blank:]]+")[[1]][2])/as.numeric(strsplit(forecast.file[grep("BIO_Smry_unfished",forecast.file)],split="[[:blank:]]+")[[1]][2])
    Quant.out[i,18]<-as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
    Quant.out[i,19]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
    Quant.out[i,20]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)+2], " ")[[1]][2])
    Quant.out[i,21]<-as.numeric(strsplit(rep.new[grep("Convergence",rep.new)], " ")[[1]][2])
    if(Dep.in[1]>=0){Quant.out[i,22]<-as.numeric(Dep.line[4])}
    if(Dep.in[1]>=0){Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3])}
    #    if(Dep.in[1]>=0){Quant.out[i,22]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num,"_SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][5])}
    #    if(Dep.in[1]>=0){Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num,"_SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][6])}
    Quant.out[i,24]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8])
    Quant.out[i,25]<-as.numeric(strsplit(rep.new[grep(paste("F_",f_yr,sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
    Quant.out[i,26]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,27]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
    Quant.out[i,28]<-as.numeric(strsplit(rep.new[grep("Crash_Pen",rep.new)], " ")[[1]][2])
    warnings.out<-readLines(paste(filepath,"/warning.sso",sep=""))
    num.warning<-as.numeric(strsplit(warnings.out[grep("N warnings",warnings.out)[1]], " ")[[1]][4])
    
    if(Dep.in[1]>=0)
    {
      # if(FMSY_M.in[1]<0 & BMSY_B0.in[1]<0&abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,28]==0&Quant.out[i,14]<100000){i<-i+1}
      # if(FMSY_M.in[1]>=0 & BMSY_B0.in[1]>=0 | h.in[1]==99 | BH_FMSY_comp==T)
      #     {
      #       if(h.in[1]==99){if(abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,28]==0&Quant.out[i,14]<100000&(abs(Quant.out[i,18]-Input.draws.MQs[i,3])/Input.draws.MQs[i,3])<0.6){i<-i+1}}
      #       else{if(abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,28]==0&Quant.out[i,14]<100000&(abs(Quant.out[i,16]-(Input.draws.MQs[i,2])/100)/(Input.draws.MQs[i,2]/100))<0.5){i<-i+1}}
      #     }
      #    else
      # {
        #Run using .par file if above conditions not met
        # print("*** RUNNING WITH .PAR ***")
        # starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
        # sum_age_line<-strsplit(starter.new[grep("ss.par",starter.new)], " ")[[1]]
        # sum_age_line[1]<-1
        # starter.new[grep("ss.par",starter.new)]<-paste(sum_age_line, collapse=" ")
        # write(starter.new,paste(filepath,"/starter.ss",sep=""))
        # RUN.SS(paste(filepath,"/",sep=""), ss.exe="ss",ss.cmd=" -nohess -nox > out.txt 2>&1")
        # rep.new<-readLines(paste(filepath,"/Report.sso",sep=""))
        # forecast.file<-readLines(paste(filepath,"/Forecast-report.sso",sep=""))
        # for(iii in 1:length(year0:sb_ofl_yrs[1])){SB.out[i,iii]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb.years[iii],sep=""),rep.new)], " ")[[1]][3])}
        # Dep.series.out<-SB.out/SB.out[,1]
        # Quant.out[i,1]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,2]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Fem_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,3]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Fem_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,4]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Fem_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,5]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,6]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Mal_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,7]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Mal_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,8]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Mal_GP_1",rep.new)], " ")[[1]][3])
        # Quant.out[i,9]<-as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][3])
        # Quant.out[i,10]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])
        # Quant.out[i,11]<-as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
        # Quant.out[i,12]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])
        # Quant.out[i,13]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
        # Quant.out[i,14]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
        # Quant.out[i,15]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
        # Quant.out[i,16]<-as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
        # Quant.out[i,17]<-as.numeric(strsplit(forecast.file[grep("calculate_FMSY",forecast.file)+13],split="[[:blank:]]+")[[1]][2])/as.numeric(strsplit(forecast.file[grep("BIO_Smry_unfished",forecast.file)],split="[[:blank:]]+")[[1]][2])
        # Quant.out[i,18]<-as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
        # Quant.out[i,19]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
        # Quant.out[i,20]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)+2], " ")[[1]][2])
        # Quant.out[i,21]<-as.numeric(strsplit(rep.new[grep("Convergence",rep.new)], " ")[[1]][2])
        # Quant.out[i,22]<-as.numeric(Dep.line[4])
        # Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3])
        # #        if(Dep.in[1]>=0){Quant.out[i,22]<-as.numeric(Dep.line[4])}
        # #        if(Dep.in[1]>=0){Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3])}
        # #        if(Dep.in[1]>=0){Quant.out[i,22]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num,"_SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][5])}
        # #        if(Dep.in[1]>=0){Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste(fleet_num,"_SURVEY1 ",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][6])}
        # Quant.out[i,24]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8])
        # Quant.out[i,25]<-as.numeric(strsplit(rep.new[grep(paste("F_",f_yr,sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
        # Quant.out[i,26]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
        # Quant.out[i,27]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
        # Quant.out[i,28]<-as.numeric(strsplit(rep.new[grep("Crash_Pen",rep.new)], " ")[[1]][2])
        # starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
        # sum_age_line<-strsplit(starter.new[grep("ss.par",starter.new)], " ")[[1]]
        # sum_age_line[1]<-0
        # starter.new[grep("ss.par",starter.new)]<-paste(sum_age_line, collapse=" ")
        # write(starter.new,paste(filepath,"/starter.ss",sep=""))

      #if(FMSY_M.in[1]<0 & BMSY_B0.in[1]<0)
        if(FMSY_M.in[1]<0 & BMSY_B0.in[1]<0&abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,28]<0.0001&Quant.out[i,14]<100000){i<-i+1}
      #if(FMSY_M.in[1]>=0 & BMSY_B0.in[1]>=0 | h.in[1]==99 | BH_FMSY_comp==T)
        if(FMSY_M.in[1]>=0 & BMSY_B0.in[1]>=0 | h.in[1]==99 | BH_FMSY_comp==T)
          {
            if(h.in[1]==99)
              {
                if(abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,28]<0.0001&Quant.out[i,14]<100000&(abs(Quant.out[i,18]-Input.draws.MQs[i,3])/Input.draws.MQs[i,3])<0.6){i<-i+1}
                #if(abs(Quant.out[i,22]-Quant.out[i,23])>0.01&Quant.out[i,28]!=0&Quant.out[i,14]>100000&(abs(Quant.out[i,18]-Input.draws.MQs[i,3])/Input.draws.MQs[i,3])>0.6)
                else
                  {
                    if(n.bad==0){Quant.out.bad<-Quant.out[i,]}
                    if(n.bad>0){Quant.out.bad<-rbind(Quant.out.bad,Quant.out[i,])}
                    n.bad<-n.bad+1
                    if(length(M.in)>6|length(Dep.in)>3|length(h.in)>3|length(FMSY_M.in)>3|length(BMSY_B0.in)>3)
                      {
                        Quant.out[i,]<-NA
                        i<-i+1
                      }
                  }    
              }

            else
            {
              if(abs(Quant.out[i,22]-Quant.out[i,23])<0.01&Quant.out[i,28]<0.0001&Quant.out[i,14]<100000&(abs(Quant.out[i,16]-(Input.draws.MQs[i,2])/100)/(Input.draws.MQs[i,2]/100))<0.5){i<-i+1}
              #if(abs(Quant.out[i,22]-Quant.out[i,23])>0.01&Quant.out[i,28]!=0&Quant.out[i,14]>100000&(abs(Quant.out[i,16]-(Input.draws.MQs[i,2])/100)/(Input.draws.MQs[i,2]/100))>0.5)
              else
                {
                   if(n.bad==0){Quant.out.bad<-Quant.out[i,]}
                   if(n.bad>0){Quant.out.bad<-rbind(Quant.out.bad,Quant.out[i,])}
                   n.bad<-n.bad+1
                   if(length(M.in)>6|length(Dep.in)>3|length(h.in)>3|length(FMSY_M.in)>3|length(BMSY_B0.in)>3)
                     {
                       Quant.out[i,]<-NA
                       i<-i+1
                     }
                }
            }
          }
                
        else
        {
          if(n.bad==0){Quant.out.bad<-Quant.out[i,]}
          if(n.bad>0){Quant.out.bad<-rbind(Quant.out.bad,Quant.out[i,])}
          n.bad<-n.bad+1
          if(length(M.in)>6|length(Dep.in)>3|length(h.in)>3|length(FMSY_M.in)>3|length(BMSY_B0.in)>3)
            {
              Quant.out[i,]<-NA
              i<-i+1
            }
        }
      # }
    }
    
    if(Dep.in[1]<0)
    {
      if(Quant.out[i,28]<0.0001&Quant.out[i,14]<100000){i<-i+1}
      else
      {
        #Run using .par file if above conditions not met
        print("*** RUNNING WITH .PAR ***")
        starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
        par_line<-strsplit(starter.new[grep("ss.par",starter.new)], " ")[[1]]
        par_line[1]<-1
        starter.new[grep("ss.par",starter.new)]<-paste(par_line, collapse=" ")
        write(starter.new,paste(filepath,"/starter.ss",sep=""))
        RUN.SS(paste(filepath,"/",sep=""), ss.exe="ss",ss.cmd=" -nohess -nox > out.txt 2>&1")
        rep.new<-readLines(paste(filepath,"/Report.sso",sep=""))
        forecast.file<-readLines(paste(filepath,"/Forecast-report.sso",sep=""))
      for(iii in 1:length(year0:sb_ofl_yrs[1]))
    	{
    		SB.out[i,iii]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb.years[iii],sep=""),rep.new)], " ")[[1]][3])
    		TB.out[i,iii]<-as.numeric(strsplit(rep.new[grep("TIME_SERIES",rep.new)+3+iii], " ")[[2]][5])
    		SumAge.out[i,iii]<-as.numeric(strsplit(rep.new[grep("TIME_SERIES",rep.new)+3+iii], " ")[[2]][6])
    		SPR.out[i,iii]<-as.numeric(strsplit(rep.new[grep("SPR_series",rep.new)+5+iii], " ")[[2]][8])
    	}
     Dep.series.out<-SB.out/SB.out[,1]
        Quant.out[i,1]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Fem_GP_1",rep.new)], " ")[[1]][3])
        Quant.out[i,2]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Fem_GP_1",rep.new)], " ")[[1]][3])
        Quant.out[i,3]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Fem_GP_1",rep.new)], " ")[[1]][3])
        Quant.out[i,4]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Fem_GP_1",rep.new)], " ")[[1]][3])
    if(genders==T)
    {
        Quant.out[i,5]<-as.numeric(strsplit(rep.new[grep("NatM_p_1_Mal_GP_1",rep.new)], " ")[[1]][3])
        Quant.out[i,6]<-as.numeric(strsplit(rep.new[grep("L_at_Amin_Mal_GP_1",rep.new)], " ")[[1]][3])
        Quant.out[i,7]<-as.numeric(strsplit(rep.new[grep("L_at_Amax_Mal_GP_1",rep.new)], " ")[[1]][3])
        Quant.out[i,8]<-as.numeric(strsplit(rep.new[grep("VonBert_K_Mal_GP_1",rep.new)], " ")[[1]][3])
    }
        Quant.out[i,9]<-as.numeric(strsplit(rep.new[grep("steep",rep.new)], " ")[[1]][3])
        Quant.out[i,10]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][3])
        Quant.out[i,11]<-as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
        Quant.out[i,12]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])
        Quant.out[i,13]<-as.numeric(strsplit(rep.new[grep(paste("SSB_",sb_ofl_yrs[1],sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
        Quant.out[i,14]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
        Quant.out[i,15]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[2],sep=""),rep.new)], " ")[[1]][3])
        Quant.out[i,16]<-as.numeric(strsplit(rep.new[grep("SSB_MSY",rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("SSB_Initial",rep.new)], " ")[[1]][3])
        Quant.out[i,17]<-as.numeric(strsplit(forecast.file[grep("calculate_FMSY",forecast.file)+13],split="[[:blank:]]+")[[1]][2])/as.numeric(strsplit(forecast.file[grep("BIO_Smry_unfished",forecast.file)],split="[[:blank:]]+")[[1]][2])
        Quant.out[i,18]<-as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
        Quant.out[i,19]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)], " ")[[1]][2])
        Quant.out[i,20]<-as.numeric(strsplit(rep.new[grep("TOTAL",rep.new)+2], " ")[[1]][2])
        Quant.out[i,21]<-as.numeric(strsplit(rep.new[grep("Convergence",rep.new)], " ")[[1]][2])
        Quant.out[i,22]<-as.numeric(Dep.line[4])
        Quant.out[i,23]<-as.numeric(strsplit(rep.new[grep(paste("Bratio_",as.numeric(Dep.line[1]),sep=""),rep.new)], " ")[[1]][3])
        Quant.out[i,24]<-as.numeric(strsplit(rep.new[grep("R0",rep.new)], " ")[[1]][8])
        Quant.out[i,25]<-as.numeric(strsplit(rep.new[grep(paste("F_",f_yr,sep=""),rep.new)], " ")[[1]][3])/as.numeric(strsplit(rep.new[grep("Fstd_MSY",rep.new)], " ")[[1]][3])
        Quant.out[i,26]<-as.numeric(strsplit(rep.new[grep(paste("OFLCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
        Quant.out[i,27]<-as.numeric(strsplit(rep.new[grep(paste("ForeCatch_",sb_ofl_yrs[3],sep=""),rep.new)], " ")[[1]][3])
        Quant.out[i,28]<-as.numeric(strsplit(rep.new[grep("Crash_Pen",rep.new)], " ")[[1]][2])
        starter.new<-readLines(paste(filepath,"/starter.ss",sep=""))
        par_line<-strsplit(starter.new[grep("ss.par",starter.new)], " ")[[1]]
        par_line[1]<-0
        starter.new[grep("ss.par",starter.new)]<-paste(par_line, collapse=" ")
        write(starter.new,paste(filepath,"/starter.ss",sep=""))
        if(Quant.out[i,28]<0.0001&Quant.out[i,14]<100000){i<-i+1}
        else
        {
          if(n.bad==0){Quant.out.bad<-Quant.out[i,]}
          if(n.bad>0){Quant.out.bad<-rbind(Quant.out.bad,Quant.out[i,])}
          n.bad<-n.bad+1
        }
      }
    }
    print(paste0("Kept ",i-1))
    ii<-ii+1
    if(file.exists(paste(filepath,"/Report.sso",sep=""))){file.remove(paste(filepath,"/Forecast-report.sso",sep=""))}
    if(file.exists(paste(filepath,"/Report.sso",sep=""))){file.remove(paste(filepath,"/Report.sso",sep=""))}
  }
  colnames(Quant.out)<-colnames(Quant.out.bad)<-c("M_f","L1_f","Linf_f","k_f","M_m","L1_m","Linf_m","k_m","h","R0","SB0",paste("SSB_",sb_ofl_yrs[1],sep=""),"Term_Yr_Dep",paste("OFL_",sb_ofl_yrs[2],sep=""),paste("AdjCatch_",sb_ofl_yrs[2],sep=""),"SBMSY/SB0","BMSY/B0","FMSY","-lnL","LL_survey","Gradient","Dep.Obs","Dep.Exp","R0_init",paste("F_",f_yr,"/FMSY",sep=""),paste("OFL_",sb_ofl_yrs[3],sep=""),paste("AdjCatch_",sb_ofl_yrs[3],sep=""),"Crash_penalty")
  #if(all(is.na(Input.draws.M))=="FALSE"){)
  if(ncol(Input.draws)==10){colnames(Input.draws)<-c("h","Dep","M_f","L1_f","Linf_f","k_f","M_m","L1_m","Linf_m","k_m")}
  if(ncol(Input.draws)==11){colnames(Input.draws)<-c("Sfrac","Dep","M_f","L1_f","Linf_f","k_f","Beta","M_m","L1_m","Linf_m","k_m")}
  if(SR_type>=8 | h.in[1]==99 | BH_FMSY_comp==T)
    {
      Input.draws<-cbind(Input.draws,Input.draws.MQs,SR_expo.out,SRMQ_par)
      ltcolnames<-length(colnames(Input.draws))
      colnames(Input.draws)[c(ltcolnames-1,ltcolnames)]<-c("Beta","Obj_fxn")
    }
  end.time<-Sys.time()
  Spp.quant.out<-list(Input.draws,Quant.out,SB.out,Dep.series.out,TB.out,SumAge.out,SPR.out,Quant.out.bad,ii-1,(as.numeric(end.time)-as.numeric(start.time))/60)
  names(Spp.quant.out)<-c("Priors","Posteriors","SB_series","Rel_Stock_status_series","Total_Biomass","Summary_Biomass","SPR","Rejected_draws","Total draws","Runtime_minutes")
  SSS.out<-Spp.quant.out
  save(SSS.out,file=paste(filepath,"/SSS_out.DMP",sep=""))
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

