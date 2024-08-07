
rm(list=ls())
n<-function(x){sum(!is.na(x))}
md<-function(x){c(m=mean(x,na.rm=T),sd=sd(x,na.rm=T))}
me<-function(x){c(m=mean(x,na.rm=T),se=sd(x,na.rm=T)/sqrt(sum(!is.na(x))))}
sd.sum<-function(sd){
	sqrt(sum(sd^2,na.rm=TRUE))}

setwd('/R_code')
gc<-read.table('GlobalCarbonPartitioning_Luyssaert_2007.txt',sep='\t',head=TRUE)
gc$NPP.area<-gc$NPP*gc$Area
gc$NPPw.area<-gc$NPPw*gc$Area
gc$NPPf.area<-gc$NPPf*gc$Area

gc$NPP.sd.area<-gc$NPP.sd*gc$Area
gc$NPPw.sd.area<-gc$NPPw.sd*gc$Area
gc$NPPf.sd.area<-gc$NPPf.sd*gc$Area

gc$rNPPf.sd<-
gc$NPPf/gc$NPP*(gc$NPPf.sd/gc$NPPf)

sum(gc$NPP.area)/sum(gc$Area)
sd.sum(gc$NPP.sd.area)/sum(gc$Area)
sum(gc$NPPf.area)/sum(gc$Area)
sd.sum(gc$NPPf.sd.area)/sum(gc$Area)
sum(gc$NPPw.area)/sum(gc$Area)
sd.sum(gc$NPPw.sd.area)/sum(gc$Area)

sum(gc$NPPf.area)/sum(gc$NPP.area)
(sum(gc$NPPf.area)/sum(gc$NPP.area))*(sd.sum(gc$NPPf.sd.area)/sum(gc$NPPf.area))

sum(gc$NPPw.area)/sum(gc$NPP.area)
(sum(gc$NPPw.area)/sum(gc$NPP.area))*(sd.sum(gc$NPPw.sd.area)/sum(gc$NPPw.area))


smg<-read.table('BT_Stand_scaled.txt',sep='\t',head=TRUE)
smg$std.sp<-smg$std
smg$std<-ifelse(!is.na(smg$std2),smg$std2,smg$std.sp)
smg$ttrt<-paste(smg$trt1,smg$trt2,sep='')
smg$std.mx<-exp(log((smg$Bbr*10/smg$std)/(2.528*1.1))*(1/-1.515))*10000
smg$sdi<-(smg$std/2.4711)*(smg$d*100/25.4)^ifelse(smg$spe%in%c('P. abies','P. sylvestris','P. tremula × tremuloides','B. pendula'),1.66,ifelse(smg$spe=='B. gymnorrhiza',1.6175,1.605))
smg$rsd1<-ifelse(smg$spe=='B. gymnorrhiza',smg$std,smg$sdi)/ifelse(smg$spe=='B. gymnorrhiza',smg$std.mx,ifelse(smg$spe=='P. taeda',450,400))
smg$rsd<-ifelse(smg$spe=='B. gymnorrhiza',smg$std,smg$sdi)/ifelse(smg$spe=='B. gymnorrhiza',smg$std.mx,400)

smg$NPPbr<-smg$NPPbr0+smg$brt
smg$NPPaw<-smg$NPPaw0+smg$brt
smg$brtr1<-smg$brt/smg$Bbr.y
smg$brtr.nl<-smg$brt/smg$Bbr.y
smg$h.d.lc<-smg$h.d/smg$lc
smg$rc<-smg$lc/smg$h
smg$bc.d.lc<-smg$bc.d/smg$lc
smg$dh<-smg$d/smg$h
smg$rh<-smg$h.d/smg$h
smg$rNPPbr<-smg$brt/smg$NPPbr
smg$rNPPaw<-smg$brt/smg$NPPaw
smg$rNPPaw0<-smg$NPPaw0/smg$NPPaw
smg$rNPPaw01<-smg$NPPaw/smg$NPPaw0
smg$rev.ord<-(10)-smg$ord.spe
smg$Bbr.0<-smg$Bbr.y-smg$brt
smg$brt.within.crown<-(smg$brt*0.16/0.84)

## Estimates (mean, standard deviation, range of data)
sd(smg$brtr)/mean(smg$brtr)
md(smg$brtr)
range(smg$brtr)
range(smg$rNPPbr)
range(smg$rNPPaw)
md(smg$rNPPbr)
md(smg$rNPPaw)
md(smg$rNPPaw01)

# For each biome type
md(smg$rNPPaw01[smg$site%in%c('TECHS22','TECHS20','TECHS30','ok')])
md(smg$rNPPaw01[smg$site%in%c('dk','fv','ml','mt','sl','nc')])
md(smg$rNPPaw01[smg$site%in%c('hb','ro','Bräc','Gävl','Grän','Möln','Ebbe','ja','fl','fr','rd')])

# For each of shade-tolerant and intolerant species types
range(smg$brtr[smg$sp==1])
range(smg$brtr[smg$sp==2])

######################################################
## Logistic prediction ##############################
h.d<-function(h.d,sp){
	mapply(function(h.d,sp){
		if(sp==1){
			brtr.h.d<-nls(brtr~c+((d-c)/(1+exp(-a*(h.d-b)))),smg[smg$sp==1,],start=list(a=1.74,b=2.55,c=.058,d=.73))
				a<-summary(brtr.h.d)$coe[1]
				b<-summary(brtr.h.d)$coe[2]
				c<-summary(brtr.h.d)$coe[3]
				d<-summary(brtr.h.d)$coe[4]
			return(c+((d-c)/(1+exp(-a*(h.d-b)))))
				}
		if(sp==2){
		brtr.h.d<-lm(brtr~h.d,smg[smg$sp==2,])
				coe<-summary(brtr.h.d)$coe
			return(coe[1]+coe[2]*h.d)
				}
							},h.d,sp)
					}


## Prediction and residuals
smg$h.d.pre<-h.d(smg$h.d,smg$sp)
smg$h.d.res<-smg$brtr-smg$h.d.pre
smg$h.d.res.r<-smg$brtr/smg$h.d.pre
summary(lm(h.d.pre~brtr,smg[smg$sp==1,]))
## Relative residulas vs. std
h.d.std.r<-function(std,rsd,sp,nl){
	mapply(function(std,rsd,sp,nl){
	if(sp==1){
		if(nl==1){
		coe<-
		nls(h.d.res.r~c+((d-c)/(1+exp(-a*(rsd-b)))),smg[smg$sp==1,],start=list(a=6.06,b=.7212,c=.5796,d=1.6384))
			a<-summary(coe)$coe[1]
			b<-summary(coe)$coe[2]
			c<-summary(coe)$coe[3]
			d<-summary(coe)$coe[4]
			return(c+((d-c)/(1+exp(-a*(rsd-b)))))
				}
		if(nl==2){
		coe<-summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))$coe
		return(coe[1]+(rsd)*coe[2])
				}
				}
	if(sp==2){
	coe<-summary(lm(h.d.res.r~std,smg[smg$sp==2,]))$coe
	return(coe[1]+(std)*coe[2])}	
	},std,rsd,sp,nl)
	}

smg$rsd.pre<-h.d.std.r(smg$std,smg$rsd,smg$sp,2)
summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))
smg$rsd.res<-smg$h.d.res.r-smg$rsd.pre


# Residual distribution

########## Whole model ################
h.d.std<-function(h.d,sp,std,rsd,nl){
	mapply(function(h.d,sp,std,rsd,nl){
		if(sp==1){
			brtr.h.d<-nls(brtr~c+((d-c)/(1+exp(-a*(h.d-b)))),smg[smg$sp==1,],start=list(a=2.22,b=2.55,c=.09,d=.81))
				a<-summary(brtr.h.d)$coe[1]
				b<-summary(brtr.h.d)$coe[2]
				c<-summary(brtr.h.d)$coe[3]
				d<-summary(brtr.h.d)$coe[4]
		if(nl==1){
		coe<-
		nls(h.d.res.r~c+((d-c)/(1+exp(-a*(rsd-b)))),smg[smg$sp==1,],start=list(a=6.06,b=.7212,c=.5796,d=1.6384))
			a1<-summary(coe)$coe[1]
			b1<-summary(coe)$coe[2]
			c1<-summary(coe)$coe[3]
			d1<-summary(coe)$coe[4]
		return((c+((d-c)/(1+exp(-a*(h.d-b)))))*(c1+((d1-c1)/(1+exp(-a1*(rsd-b1))))))
				}
		if(nl==2){
		coe<-summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))$coe
		return((c+((d-c)/(1+exp(-a*(h.d-b)))))*(coe[1]+(rsd)*coe[2]))
				}
				}

		if(sp==2){
		brtr.h.d<-lm(brtr~h.d,smg[smg$sp==2,])
				coe1<-summary(brtr.h.d)$coe
		brtr.res.r.std2<-lm(h.d.res.r~std,smg[smg$sp==2,])
				coe2<-summary(brtr.res.r.std2)$coe
			return((coe1[1]+coe1[2]*h.d)*(coe2[1]+coe2[2]*std))
				}
							},h.d,sp,std,rsd,nl)
					}

smg$pre.brtr<-h.d.std(smg$h.d,smg$sp,smg$std,smg$rsd,2)
smg$res.brtr<-smg$brtr-smg$pre.brtr
smg$res.brtr.r<-smg$brtr/smg$pre.brtr


## Figure S1 residual distribution of the developed model in Figure 3 
quartz(w=3.42,h=1.8)
par(mfrow=c(1,2))
par(lwd=.3)
par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,1.1),ylim=c(-0.3,0.3))
lines(c(-1,2),c(0,0),col=8,lwd=7/12,lty=3)
points(res.brtr~pre.brtr,smg,col=col,pch=pch,bg=bg,cex=.5,lwd=.2)
points(res.brtr~pre.brtr,smg[smg$sp==2&smg$site=='dk',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)
points(res.brtr~pre.brtr,smg[smg$sp==2&smg$site=='fr',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)	
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
axis (1,seq(-2,2,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,2,by=.2),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Residuals')),2,line=.7,font=1,cex=7/12)
mtext(expression(paste('Predicted ‡'[B])),1,line=.5,font=1,cex=7/12)

par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,1.1),ylim=c(0,3.5))
lines(c(-1,2),c(1,1),col=8,lwd=7/12,lty=3)
points(res.brtr.r~pre.brtr,smg,col=col,pch=pch,bg=bg,cex=.5,lwd=.2)
points(res.brtr.r~pre.brtr,smg[smg$sp==2&smg$site=='dk',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)
points(res.brtr.r~pre.brtr,smg[smg$sp==2&smg$site=='fr',],pch=pch,bg=bg,col=col,lwd=1,cex=.75)	
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
axis (1,seq(-2,2,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,5,by=1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Normalized residuals')),2,line=.7,font=1,cex=7/12)

## R-square
summary(lm(h.d.res.r~rsd.pre*as.factor(sp),smg))$r.s+(1-summary(lm(h.d.res.r~rsd.pre*as.factor(sp),smg))$r.s)*summary(lm(brtr~h.d.pre*as.factor(sp),smg))$r.s


################################################################################################
## Impact at stand scale
## models for NPPaw with/without branch turnover
################################################################################################
smg1<-smg[smg$sp==1,]
smg2<-smg[smg$sp==2,]

smg1$res.rc<-residuals(lm(rNPPbr~rc,smg1))
summary(lm(rNPPbr~rc*as.factor(sp),smg))
summary(lm(rNPPaw~rc*as.factor(sp),smg))

rNPPbr.rc<-lm(rNPPbr~rc+as.factor(sp),smg)
summary(rNPPbr.rc)
rNPPaw.rc<-lm(rNPPaw~rc+as.factor(sp),smg)
summary(rNPPaw.rc)

md(smg$brtr)
range(smg$brtr)
md(smg$rNPPbr)
range(smg$rNPPbr)

f.rNPPbr.rc<-function(rc,sp){
	x<-rc
	coe<-summary(rNPPbr.rc)$coe
	coe[1]+coe[3]*(sp-1)+coe[2]*(x)}
f.rNPPaw.rc<-function(rc,sp){
	x<-rc
	coe<-summary(rNPPaw.rc)$coe
	coe[1]+coe[3]*(sp-1)+(coe[2])*(x)}


Brtr.mod1<-data.frame(summary(nls(brtr~c+((d-c)/(1+exp(-a*(h.d-b)))),smg[smg$sp==1,],start=list(a=2.22,b=2.55,c=.09,d=.81)))$coe)
Brtr.mod1$mod<-'Brtr.mod1'
Brtr.res.mod11<-data.frame(summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))$coe)
Brtr.res.mod11$mod<-'Brtr.res.mod11'
Brtr.res.mod12<-data.frame(summary(nls(h.d.res.r~c+((d-c)/(1+exp(-a*(rsd-b)))),smg[smg$sp==1,],start=list(a=6.06,b=.7212,c=.5796,d=1.6384)))$coe)
Brtr.res.mod12$mod<-'Brtr.res.mod12'



Brtr.mod2<-data.frame(summary(lm(brtr~h.d,smg[smg$sp==2,]))$coe)
Brtr.mod2$mod<-'Brtr.mod2'
Brtr.res.mod2<-data.frame(summary(lm(h.d.res.r~std,smg[smg$sp==2,]))$coe)
Brtr.res.mod2$mod<-'Brtr.res.mod2'

rNPPbr.mod<-data.frame(summary(rNPPbr.rc)$coe)
rNPPbr.mod$mod<-'rNPPbr.mod'

rNPPaw.mod<-data.frame(summary(rNPPaw.rc)$coe)
rNPPaw.mod$mod<-'rNPPaw.mod'

parameters.scaled<-rbind(Brtr.mod1, Brtr.res.mod11, Brtr.res.mod12, Brtr.mod2, Brtr.res.mod2, rNPPbr.mod, rNPPaw.mod)

write.table(parameters.scaled,'BT_Stand_scaled_parameters.txt',sep='\t',quote=FALSE,row.names=TRUE)

### FIGURE 3

quartz(w=3.42,h=3.6)
par(mfrow=c(2,2))
par(lwd=.3)
par(mai=c(.2,.4,.4,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,5),ylim=c(0,1))
points(brtr~h.d,smg[smg$sp==1,],col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
mtext('Shade-intolerant species',line=-.1,adj=1,font=1,cex=7/12)
summary(lm(h.d.pre~brtr,smg[smg$sp==1,]))
mtext(expression(paste('r'^2,'=0.91; p<.001',sep='')),1,line=-1,adj=1,font=1,cex=7/12)

axis (1,seq(-2,9,by=2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-1,2,by=.5),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Annual turnover rate of branch biomass (y'^-1,')')),2,line=.6,font=1,cex=7/12,at=-.2)
mtext(expression(paste('Height increment (m y'^-1,')')),1,line=.5,font=1,cex=7/12)

for(i in 1){
	df<-smg[smg$sp==i,]
	site<-unique(df$site[df$site!='hb'])
	site<-site[13:1]
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(brtr+h.d~bg+pch+col,df2,FUN=me)
			lines(rep(df3$h.d.m,2),c(-1,1)*df3$brtr.se+df3$brtr.m)
			lines(c(-1,1)*df3$h.d.se+df3$h.d.m,rep(df3$brtr.m,2))
			points(brtr.m~h.d.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
			}}}
		points(brtr~h.d,smg[smg$site=='hb',],bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
	curve(h.d(x,1),add=TRUE,xlim=range(smg$h.d[smg$sp==1],na.rm=TRUE),lwd=0.75)

par(mai=c(.2,.35,.4,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(-0.1,2.5),xlim=c(0,1.45))
points(h.d.res.r~rsd,smg[smg$sp==1,],col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
summary(lm(h.d.res.r~rsd,smg[smg$sp==1,]))
mtext(expression(paste('r'^2,'=0.51; p<.001')),1,line=-1,adj=1,font=1,cex=7/12)
axis (1,seq(-2,8,by=.5),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-2,8,by=1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Relative residuals (observed / predicted)')),2,line=.6,font=1,cex=7/12,at=-.5)
mtext(expression(paste('Relative stand density')),1,line=.5,font=1,cex=7/12)
for(i in 1){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df2$ln.std<-log(df2$std)
			df3<-summaryBy(h.d.res.r+ln.std+std+rsd~bg+pch+col,df2,FUN=me)
			lines(rep(df3$rsd.m,2),c(-1,1)*df3$h.d.res.r.se+df3$h.d.res.r.m)
			lines(c(-1,1)*df3$rsd.se+df3$rsd.m,rep(df3$h.d.res.r.m,2))
			points(h.d.res.r.m~rsd.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=0.2)
			}}}
curve(h.d.std.r(NA,x,1,2),xlim=range(smg$rsd[smg$sp==1],na.rm=TRUE),add=TRUE,lwd=0.75)

par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,0.9),ylim=c(-.01,.09))
mtext('c',line=-.8,adj=0.02,font=2,cex=8/12)
mtext('Shade-tolerant species',line=-.1,adj=1,font=1,cex=7/12)
summary(lm(brtr~h.d,smg[smg$sp==2,]))
mtext(expression(paste('r'^2,'=0.43; p<.001',sep='')),1,line=-1,adj=1,font=1,cex=7/12)
points(brtr~h.d,smg[smg$sp==2,],col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
axis (1,seq(-2,9,by=0.2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-0.4,1,by=0.04),,tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Height increment (m y'^-1,')')),1,line=.5,font=1,cex=7/12)

for(i in 2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(brtr+h.d~bg+pch+col,df2,FUN=me)
			lines(rep(df3$h.d.m,2),c(-1,1)*df3$brtr.se+df3$brtr.m)
			lines(c(-1,1)*df3$h.d.se+df3$h.d.m,rep(df3$brtr.m,2))
			points(brtr.m~h.d.m,df3,bg=bg,pch=pch,col=col,cex=5/12,lwd=0.3)
			points(brtr.m~h.d.m,df3,bg=0,pch=pch,col=1,cex=7/12,lwd=0.3)
			}}}
curve(h.d(x,2),add=TRUE,xlim=range(smg$h.d[smg$sp==2],na.rm=TRUE),lwd=0.65)


par(mai=c(.4,.35,.2,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",ylim=c(0,2.5),xlim=c(250,3400))
mtext('d',line=-.8,adj=0.02,font=2,cex=8/12)
summary(lm(h.d.res.r~std,smg[smg$sp==2,]))
mtext(expression(paste('r'^2,'=0.41; p<.001',sep='')),1,line=-1,adj=1,font=1,cex=7/12)
points(h.d.res.r~std,smg[smg$sp==2,],col=col,pch=pch,bg='white',cex=0.4,lwd=.05)
axis (1,seq(0,6000,by=1000),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=8/12,lwd=0.3)
axis (2,seq(-1,5,by=1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=8/12,lwd=0.3)
mtext(expression(paste('Stand density (stems ha'^-1,')')),1,line=.5,font=1,cex=7/12)

for(i in 2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(h.d.res.r+sdi+std~bg+pch+col,df2,FUN=me)
			lines(rep(df3$std.m,2),c(-1,1)*df3$h.d.res.r.se+df3$h.d.res.r.m)
			lines(c(-1,1)*df3$std.se+df3$std.m,rep(df3$h.d.res.r.m,2))
			points(h.d.res.r.m~std.m,df3,bg=bg,pch=pch,col=col,cex=5/12,lwd=0.3)
			points(h.d.res.r.m~std.m,df3,bg=0,pch=pch,col=1,cex=7/12,lwd=0.3)
			}}}
curve(h.d.std.r(x,NA,2,2),xlim=range(smg$std[smg$sp==2],na.rm=TRUE),add=TRUE,lwd=0.65)


## Figure S2
quartz(w=3.42,h=1.8)
par(mfrow=c(1,2))
par(lwd=.3)
par(mai=c(.4,.4,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(.15,1),ylim=c(0,1.25))
points(rNPPbr~rc,smg,col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('a',line=-.8,adj=0.02,font=2,cex=8/12)
summary(rNPPbr.rc)
mtext(expression(paste('r'^2,'=0.56; p<.001')),1,line=-.9,adj=.05,font=1,cex=7/12)
mtext(expression(paste('')),line=-1.3,adj=.95,font=1,cex=5/12)
axis (1,seq(-2,2,by=.2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,2,by=.5),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[B])),2,line=.6,font=1,cex=7/12)
mtext(expression(paste('Live crown ratio')),1,line=.5,font=1,at=1.1,cex=7/12)

for(i in 1:2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(rNPPbr+rc~bg+pch+col,df2,FUN=me)
			lines(rep(df3$rc.m,2),c(-1,1)*df3$rNPPbr.se+df3$rNPPbr.m)
			lines(c(-1,1)*df3$rc.se+df3$rc.m,rep(df3$rNPPbr.m,2))
			points(rNPPbr.m~rc.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=i*0.1+0.1)
			if(i==2){
			points(rNPPbr.m~rc.m,df3,bg=0,pch=pch,col=1,cex=7.5/12,lwd=0.1)}
			}}
curve(f.rNPPbr.rc(x,i),add=TRUE,xlim=range(smg$rc[smg$sp==i],na.rm=TRUE),lwd=0.75)}

par(mai=c(.4,.35,.2,.07))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(.1,1),ylim=c(0,.3))
points(rNPPaw~rc,smg,col=ifelse(sp==1,bg,col),pch=pch,bg='white',cex=0.4,lwd=.05)
mtext('b',line=-.8,adj=0.02,font=2,cex=8/12)
summary(rNPPaw.rc)
mtext(expression(paste('r'^2,'=0.40; p<.001')),1,line=-.9,adj=.05,font=1,cex=7/12)
mtext(expression(paste('')),line=-1.3,adj=.95,font=1,cex=5/12)
axis (1,seq(-2,2,by=.2),tck=.02,label=TRUE,mgp=c(0,-.3,0),cex.axis=7/12,lwd=0.3)
axis (2,seq(-2,2,by=.1),tck=.02,label=TRUE,mgp=c(0,0,0),cex.axis=7/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[Wa])),2,line=.6,font=1,cex=7/12)

for(i in 1:2){
	df<-smg[smg$sp==i,]
	site<-unique(df$site)
	for (j in site){
		df1<-df[df$site==j,]
		ttrt<-unique(df1$ttrt)
		for (k in ttrt){
			df2<-df1[df1$ttrt==k,]
			df3<-summaryBy(rNPPaw+rc~bg+pch+col,df2,FUN=me)
			lines(rep(df3$rc.m,2),c(-1,1)*df3$rNPPaw.se+df3$rNPPaw.m)
			lines(c(-1,1)*df3$rc.se+df3$rc.m,rep(df3$rNPPaw.m,2))
			points(rNPPaw.m~rc.m,df3,bg=bg,pch=pch,col=col,cex=7/12,lwd=i*0.1+0.1)
			if(i==2){
			points(rNPPaw.m~rc.m,df3,bg=0,pch=pch,col=1,cex=7.5/12,lwd=0.1)}
			}}
curve(f.rNPPaw.rc(x,i),add=TRUE,xlim=range(smg$rc[smg$sp==i],na.rm=TRUE),lty=1,lwd=0.75)}


## Figure 2d-f
quartz(w=4.5,h=2.55)
par(mfrow=c(1,4))
par(lwd=.1,col=0)
par(mai=c(.35,0,.2,.0))
plot(NA,xlab="",ylab="",yaxt="n",xaxt="n",xlim=c(0,10),ylim=c(-.2,9.15))
x.label<-unique(smg[,c('spe','rev.ord')])
x.label<-x.label[order(-x.label$rev.ord),]
x.label$ful.nm<-c('Pinus sylvestris','Pinus taeda','Eucalyptus grandis','Pseudotsuga menziesii','Populus tremula × tremuloides','Betula pendula','Bruguiera gymnorrhiza','Picea abies','Broad-leaved deciduous')

par(lwd=.1,col=1)
mtext(x.label$ful.nm[1:8],las=2,2,at=(x.label$rev.ord[1:8])-0.1,line=-8.5,font=3,cex=6/12)
mtext(x.label$ful.nm[9],las=2,2,at=(x.label$rev.ord[9])-0.1,line=-8.5,font=c(1),cex=6/12)
mtext('Average',las=2,2,at=-0.25,line=-8,font=c(1),cex=6/12)


par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,0.9),ylim=c(-.2,9.25))
points((jitter(rev.ord,.5)-0.5)~brtr,smg,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(smg$sp==1,bg,col))
mtext('a',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-4,4,by=0.5),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('∆'[B],' (y'^-1,')')),1,line=.9,font=1,cex=7/12)

for (i in 1:2){
	df<-smg[smg$sp==i,]
	for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$brtr[dff$brtr>=0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.85
	y.mn<-unique(dff$rev.ord)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
}}


xv<-smg$brtr[smg$brtr>=0]
shape<-mean(xv)^2/var(xv)
rate<-mean(xv)/var(xv)
t.max<-max(dnorm(seq(0.01,1,0.001),mean(xv),var(xv)))*0.85
points(mean(xv),-.3,pch=21,bg=1,cex=8/12)
lines(c(-1,1)*sqrt(var(xv))+mean(xv),rep(-.3,2))

par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,1.1),ylim=c(-.2,9.15))
points((jitter(rev.ord,.5)-0.5)~rNPPbr,smg,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(smg$sp==1,bg,col))
mtext('b',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-4,4,by=0.5),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[B])),1,line=.9,font=1,cex=7/12)

for (i in 1:2){
df<-smg[smg$sp==i,]
for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$rNPPbr[dff$rNPPbr>=0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.85
	y.mn<-unique(dff$rev.ord)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
	# print(unique(dff[,c('site','sp')]))
	# print(gamma_test(x))
	}}

xv<-smg$rNPPbr[dff$rNPPbr>=0]
shape<-mean(xv)^2/var(xv)
rate<-mean(xv)/var(xv)
t.max<-max(dnorm(seq(0.01,1,0.001),mean(xv),var(xv)))*0.85
points(mean(xv),-.3,pch=21,bg=1,cex=8/12)
lines(c(-1,1)*sqrt(var(xv))+mean(xv),rep(-.3,2))

par(lwd=.3)
par(mai=c(.35,.02,.2,.02))
plot(NA,xlab="",ylab="",xaxt="n",yaxt="n",xlim=c(0,.31),ylim=c(-.2,9.15))
points(jitter(rev.ord,.5)-.5~rNPPaw,smg,lwd=0.04,cex=2/12,bg='white',pch=21,col=ifelse(smg$sp==1,bg,col))
mtext('c',line=-1,adj=0.9,font=2,cex=8/12)
axis (1,seq(-1,1,by=.1),tck=.02,label=TRUE,mgp=c(0,-.2,0),cex.axis=9.5/12,lwd=0.3)
mtext(expression(paste('Turnover'[B],' : NPP'[Wa])),1,line=.9,font=1,cex=7/12)
mtext('Scaled plots',line=-.1,adj=1,font=1,cex=7/12)


for (i in 1:2){
df<-smg[smg$sp==i,]
for (j in 1:length(unique(df$site))){
	dff<-df[df$site==unique(df$site)[j],]
	x<-dff$rNPPaw[dff$brtr>0]
	x.mean<-mean(x,na.rm=TRUE)
	x.var<-var(x,na.rm=TRUE)
	xn<-min(x,na.rm=TRUE)
	xm<-max(x,na.rm=TRUE)
	xl<-length(!is.na(x))
	x.01<-range(x,na.rm=TRUE)
	shape<-x.mean^2/x.var
	rate<-x.mean/x.var
	col<-ifelse(i==1,dff$bg,dff$col)
	y.mx<-max(dgamma(seq(0.01,2,0.001),shape,rate))*1.85
	y.mn<-unique(dff$rev.ord)
curve(dgamma(x,shape,rate)/y.mx+(y.mn-.4),col=col,lwd=dff$lwd,lty=dff$lty,xlim=x.01,add=TRUE)
}}
xv<-smg$rNPPaw[dff$rNPPaw>=0]
shape<-mean(xv)^2/var(xv)
rate<-mean(xv)/var(xv)
t.max<-max(dnorm(seq(0.01,1,0.001),mean(xv),var(xv)))*0.85
points(mean(xv),-.3,pch=21,bg=1,cex=8/12)
lines(c(-1,1)*sqrt(var(xv))+mean(xv),rep(-.3,2))


