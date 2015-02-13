#RUN SSS
#Dep_V
Dir.in<-"C:/Users/copeja/Desktop/SSS_4_Anna"
source(paste(Dir.in,"/SSS_code.r",sep=""))
ARRA.SSS<-SSS(paste(Dir.in,"/Aurora/SSS",sep=""),
	file.name=c("ARRA_dat3.ss","ARRA_ctl5.ss"),
	reps=1000,
	M.in=c(3,0.0405,exp(-3.353),3,0.0429,exp(-3.295)),
	Dep.in=c(2,0.3089422,0.1726904),
	h.in=c(2,0.779,0.152),
	L1.in=c(0,0,0,0),
	Linf.in=c(0,0,0,0),
	k.in=c(0,0,0,0),
	Zfrac.Beta.in=c(-99,0.2,0.6,-99,0.5,2),
	R_start=6,
	sum_age=0,
	sb_ofl_yrs=c(2013,2013,2014),
	f_yr=2012,year0=1916,genders=T)


#load results from file
load(paste(Dir.in,"/Aurora/SSS/SSS.DMP",sep="")) 