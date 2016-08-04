
utilglm<-function(x , beta, family, criterion = "D"){

eta<-beta%*%t(x)

if(!is.function(family)){
	if(is.list(family)){
	stuff<-family} else{
		family2 <- get(family, mode = "function", envir = parent.frame())
		stuff<-family2()}
		} else{
	stuff<-family()}
mu<-stuff$linkinv(eta)
w<-(stuff$mu.eta(eta)^2)/stuff$variance(mu)

if(criterion=="D"){
eval<-.Call("Dcpp", x, w, PACKAGE = "acebayes")}
if(criterion=="A"){
eval<-.Call("Acpp", x, w, PACKAGE = "acebayes")}
if(criterion=="E"){
eval<-.Call("Ecpp", x, w, PACKAGE = "acebayes")}

as.vector(eval)}

pval <- function(oldeval, neweval, binary){
if(binary){
old_n<-length(oldeval)
new_n<-length(neweval)
old_sum<-sum(oldeval)
new_sum<-sum(neweval)
new_beta_sam<-rbeta(n = 10000, shape1 = 1 + new_sum, shape2 = 1 + new_n - new_sum)
out<-mean(pbeta(q = new_beta_sam, shape1 = 1 + old_sum, shape2 = 1 + old_n - old_sum))} else{
details<-as.vector(.Call( "pvalcpp", oldeval, neweval, PACKAGE = "acebayes" ))
out<-1-pt(details[1],df=details[2])}
out}

noisyexpectutil<-function(utility, B, d, i, j, Dij){
temp<-d
z<-c()
for(L in 1:length(Dij)){
temp[i,j]<-Dij[L]
z[L]<-mean(utility(d=temp, B=B))}
list(z=z,Dij=Dij)}

distmat <- function(Dij){
	.Call( "distcpp",Dij,PACKAGE = "acebayes" )
}

GPpred <- function(paras, dist, z, newDij, Dij){
as.vector(.Call( "GPpredcpp",paras, dist, z, newDij, Dij, PACKAGE = "acebayes" ))
}

FisherScoring2D<-function(par, Dij, z, dist, tol = 1e-10, max.iter = 25){

singular<-0
QQ<-length(Dij)

theta<-par

e1<-exp(theta[1])
e2<-exp(theta[2])

G<-exp(-e2*dist)
R<-dist*G
B<-e1*diag(QQ)+G
iB<-c()
try(iB<-solve(B),silent=TRUE)
if(!is.null(iB)){
M1<-iB%*%iB
M2<-iB%*%R
M3<-M2%*%iB
M1z<-as.vector(M1%*%matrix(z,ncol=1))
M3z<-as.vector(M3%*%matrix(z,ncol=1))

D1<-0.5*e1*(sum(M1z*z) - sum(diag(iB)))
D2<-0.5*e2*(sum(diag(M2)) - sum(M3z*z))
A11<--0.5*e1*e1*sum(diag(M1))
A12<-0.5*e1*e2*sum(diag(M3))
A22<--0.5*e2*e2*sum(M2*t(M2))
F<-A11*A22-A12*A12
counter<-1

	while(max(abs(c(D1,D2)))>tol & counter<max.iter & abs(F)<Inf & abs(F)>0 & !any(is.na(c(D1,D2,F)))){
	last<-theta
	theta<-theta-c(A22*D1 - A12*D2,A11*D2 - A12*D1)/F
	e1<-exp(theta[1])
	e2<-exp(theta[2])

	G<-exp(-e2*dist)
	R<-dist*G
	B<-e1*diag(QQ)+G
	iB<-c()
	try(iB<-solve(B),silent=TRUE)
	if(!is.null(iB)){
		M1<-iB%*%iB
		M2<-iB%*%R
		M3<-M2%*%iB
		M1z<-as.vector(M1%*%matrix(z,ncol=1))
		M3z<-as.vector(M3%*%matrix(z,ncol=1))

		D1<-0.5*e1*(sum(M1z*z) - sum(diag(iB)))
		D2<-0.5*e2*(sum(diag(M2)) - sum(M3z*z))
		A11<--0.5*e1*e1*sum(diag(M1))
		A12<-0.5*e1*e2*sum(diag(M3))
		A22<--0.5*e2*e2*sum(M2*t(M2))
		F<-A11*A22-A12*A12

		counter<-counter+1}   else{
		
		counter<-max.iter+2}

	}

if(counter==(max.iter+2)){
theta<-last}

} else{

singular<-1}

list(par=theta, singular=singular)}

acephase1<-function(utility, start.d, B = c(20000,1000), Q = 20, N1 = 20, lower, upper, limits = NULL, progress = FALSE, binary = FALSE){

ptm<-proc.time()[3]

if(!is.matrix(upper)){
UPPER<-upper+0*start.d} else{
UPPER<-upper}
if(!is.matrix(lower)){
LOWER<-lower+0*start.d} else{
LOWER<-lower}
DIFF<-UPPER-LOWER

if(length(limits)==0){
limits2<-function(i, j, d){
c(LOWER[i,j],sort(runif(9998))*DIFF[i,j]+LOWER[i,j],UPPER[i,j])}} else{
limits2<-limits}

n<-dim(start.d)[1]
k<-dim(start.d)[2]

DESIGN<-start.d
eval<-utility(d = start.d, B = B[1])
curr<-mean(eval)
curr2<-curr

counter<-1

best<-DESIGN
inner<-DESIGN
inner_eval<-curr

while(counter<=N1){

for(i in 1:n){
for(j in 1:k){

xxx<-as.vector(optimumLHS(n=Q,k=1))*DIFF[i,j]+LOWER[i,j]
#####xxx<-c(0,as.vector(optimumLHS(n=Q-2,k=1)),1)*DIFF[i,j]+LOWER[i,j]
yyy<-noisyexpectutil(utility=utility, B=B[2], d=DESIGN, i=i, j=j, Dij=xxx)$z

A.array<-distmat(xxx)

meany<-mean(yyy)
sdy<-sd(yyy)
zzz<-(yyy-meany)/sdy

optgp<-FisherScoring2D(par=c(0,0), Dij=xxx, z=zzz, dist=A.array)
opt<-NULL
if(optgp$singular!=1){
xxxz<-limits2(i = i, j = j, d = DESIGN)
yyyz<-meany+sdy*GPpred(paras=optgp$par,dist=A.array,z=matrix(zzz,ncol=1),newDij=xxxz,Dij=xxx)

opt<-NULL
start<-xxxz[max(yyyz)==yyyz]
if(length(start)==1){
opt<-list(best.x=start)}}

if(!is.null(opt)){
	old_DESIGN<-DESIGN
	new_DESIGN<-DESIGN
	new_DESIGN[i,j]<-opt$best.x


	old_eval<-utility(d=old_DESIGN,B=B[1])
	new_eval<-utility(d=new_DESIGN,B=B[1])

	the.p.val<-pval(old_eval,new_eval,binary)
	the.p.val<-ifelse(is.na(the.p.val),0,the.p.val)

	if(the.p.val>=runif(1)){
		DESIGN<-new_DESIGN
		curr<-c(curr,mean(new_eval))

			if(curr[length(curr)]>inner_eval){
			inner_eval<-curr[length(curr)]
			inner<-new_DESIGN}

		} else{
		DESIGN<-old_DESIGN
		curr<-c(curr,curr[length(curr)])}

	} else{
	curr<-c(curr,curr[length(curr)])}

}}     						#### End of scan through all elements

old_DESIGN<-best
new_DESIGN<-inner

old_eval<-utility(d=old_DESIGN,B=B[1])
new_eval<-utility(d=new_DESIGN,B=B[1])

the.p.val<-pval(old_eval,new_eval,binary)
the.p.val<-ifelse(is.na(the.p.val),0,the.p.val)

if(the.p.val>=runif(1)){
best_eval<-mean(new_eval)
best<-inner} else{
inner<-best
best_eval<-mean(old_eval)
inner_eval<-best_eval}

curr2<-c(curr2,best_eval)

if(progress){
cat("Phase I iteration ", counter, " out of ",N1," (Current value = ",best_eval,") \n",sep="")}
counter<-counter+1
}
					
ptm<-proc.time()[3]-ptm

output<-list(start.d = start.d, phase1.d = best, phase2.d = best, phase1.trace = curr2, phase2.trace = NULL, Q = Q, N1 = N1, N2 = 0, glm = FALSE, criterion = "NA", family = "NA", prior = "NA", time = ptm, binary = binary)

class(output)<-"ace"

output}

acephase2<-function(utility, start.d, B = c(20000,1000), N2 = 100, progress = FALSE, binary = FALSE){

ptm<-proc.time()[3]

n<-dim(start.d)[1]
k<-dim(start.d)[2]

DESIGN<-start.d
CAND<-DESIGN

eval<-utility(d=DESIGN,B=B[1])

best<-DESIGN
best_ob<-eval

curr2<-mean(eval)
counter2<-1
if(progress){
cat("Phase II iteration ", counter2, " out of ",N2," (Current value = ",curr2,") \n",sep="")}
while(counter2<N2){

crt<-c()
for(j in 1:n){
INTER<-rbind(DESIGN,CAND[j,])
crt[j]<-mean(utility(d=INTER,B=B[2]))}

#INTER<-rbind(DESIGN,CAND[(1:n)[max(crt)==crt],])
potpts<-(1:n)[max(crt)==crt]
if(length(potpts)>1){
potpts<-sample(x=potpts,size=1)}   ### To get around situations where augmenting two different design pts leads to same expected utility (quite common with binary = TRUE)
INTER<-rbind(DESIGN,CAND[potpts,])

crt<-c()
for(j in 1:(n+1)){
INTER2<-as.matrix(INTER[-j,],nrow=n)
crt[j]<-mean(utility(d=INTER2,B=B[2]))}

#new_DESIGN<-matrix(INTER[-((1:(n+1))[max(crt)==crt]),],nrow=n,dimnames=dimnames(CAND))
potpts<-(1:(n+1))[max(crt)==crt]
if(length(potpts)>1){
potpts<-sample(x=potpts,size=1)}   
new_DESIGN<-matrix(INTER[-potpts,],nrow=n,dimnames=dimnames(CAND))

old_DESIGN<-DESIGN

old_eval<-utility(d=old_DESIGN,B=B[1])
new_eval<-utility(d=new_DESIGN,B=B[1])

the.p.val<-pval(old_eval,new_eval,binary)
the.p.val<-ifelse(is.na(the.p.val),0,the.p.val)

if(the.p.val>runif(1)){
	curr2<-c(curr2,mean(new_eval))
	DESIGN<-new_DESIGN
	
	the.p.val<-pval(best_ob,new_eval,binary)
	the.p.val<-ifelse(is.na(the.p.val),0,the.p.val)

	if(the.p.val>=runif(1)){
		best<-DESIGN
		best_ob<-new_eval}
	} else{
	curr2<-c(curr2,mean(old_eval))}

counter2<-counter2+1
if(progress){
cat("Phase II iteration ", counter2, " out of ",N2," (Current value = ",curr2[length(curr2)],") \n",sep="")}
}

ptm<-proc.time()[3]-ptm

output<-list(start.d = start.d, phase1.d = start.d, phase2.d = best, phase1.trace = NULL, phase2.trace = curr2, Q = NULL, N1 = 0, N2 = N2, glm = FALSE, criterion = "NA", family = "NA", prior = "NA", time = ptm, binary = binary)

class(output)<-"ace"

output}


ace<-function(utility, start.d, B = c(20000,1000), Q = 20, N1 = 20, N2 = 100, lower = -1, upper = 1, limits = NULL, progress = FALSE, binary = FALSE){

ptm<-proc.time()[3]

if(N1>0){
interim<-acephase1(utility = utility, start.d = start.d, B = B, Q = Q, N1 = N1, lower = lower, upper = upper, limits = limits, progress = progress, binary = binary)
interim.d<-interim$phase1.d
interim.trace<-interim$phase1.trace} else{
interim.d<-start.d
interim.trace<-NULL}

if(N2>0){
last<-acephase2(utility = utility, start.d = interim.d, B = B, N2 = N2, progress = progress, binary = binary)
last.d<-last$phase2.d
last.trace<-last$phase2.trace} else{
last.d<-interim.d
last.trace<-NULL}

ptm<-proc.time()[3]-ptm

output<-list(start.d = start.d, phase1.d = interim.d, phase2.d = last.d, phase1.trace = interim.trace, phase2.trace = last.trace, Q = Q, N1 = N1, N2 = N2, glm = FALSE, criterion = "NA", prior = "NA", time = ptm, binary = binary)

class(output)<-"ace"

output}

aceglm<-function(formula, start.d, family, prior, criterion = "D", B = c(20000,1000),  Q = 20, N1 = 20, N2 = 100, 
lower = -1, upper = 1, progress = FALSE){

inte<-function(d, B){
x<-model.matrix(object = formula, data = data.frame(d))
beta<-prior(B)
utilglm(x = x, beta = beta, family = family, criterion = criterion)}

output<-ace(utility = inte, start.d = start.d, B = B, Q = Q, N1 = N1, N2 = N2, lower = lower, upper = upper, progress = progress)
output$glm<-TRUE
output$criterion<-criterion
#output$family<-family
output$prior<-prior
output$phase1.d<-output$phase1.d
output$phase2.d<-output$phase2.d

output}

plot.ace<-function(x,...){

if(length(x$phase1.trace)>0 & length(x$phase2.trace)>0){
ulim<-max(c(x$phase1.trace,x$phase2.trace))
llim<-min(c(x$phase1.trace,x$phase2.trace))
plot(1:length(x$phase1.trace),x$phase1.trace,xlab="Phase I iteration",ylab="Observation of expected utility",ylim=c(llim,ulim),type="l",xlim=c(0,length(x$phase1.trace)))
new_z<-axTicks(side=1)
new_x<-(1:length(x$phase2.trace))*(length(x$phase1.trace)/length(x$phase2.trace))
lines(new_x,x$phase2.trace,col=2)
legend(x="bottomright",legend=c("Phase I","Phase II"),col=c(1,2),lty=c(1,1),bty="n")
axis(side=3,labels=new_z/(length(x$phase1.trace)/length(x$phase2.trace)),at=new_z)
mtext("Phase II iteration", side=3, line = par("mgp")[1]) }

if(length(x$phase1.trace)>0 & length(x$phase2.trace)==0){
plot(1:length(x$phase1.trace),x$phase1.trace,type="l",xlab="Phase I iteration",ylab="Observation of expected utility")
legend(x="bottomright",legend=c("Phase I"),col=1,lty=1,bty="n") }

if(length(x$phase1.trace)==0 & length(x$phase2.trace)>0){
plot(1:length(x$phase2.trace),x$phase1.trace,type="l",xlab="Phase II iteration",ylab="Observation of expected utility",col=2)
legend(x="bottomright",legend=c("Phase II"),col=2,lty=1,bty="n") }

}

print.ace<-function(x,...){

hrs<-round(x$time%/%3600,0)
mins<-round((x$time%%3600)%/%60,0)
secs<-round((x$time%%3600)%%60,0)
hrs<-ifelse(hrs<10,paste("0",hrs,sep=""),hrs)
mins<-ifelse(mins<10,paste("0",mins,sep=""),mins)
secs<-ifelse(secs<10,paste("0",secs,sep=""),secs)

if(x$glm==TRUE){
cat("Generalised Linear Model \n")
cat("Criterion = Bayesian ",x$criterion,"-optimality \n",sep="")} else{
cat("User-defined utility \n")}
cat("\n")
cat("Number of runs = ",dim(x$phase2.d)[1],"\n",sep="")
cat("\n")
cat("Number of factors = ",dim(x$phase2.d)[2],"\n",sep="")
cat("\n")
cat("Number of Phase I iterations = ",x$N1,"\n",sep="")
cat("\n")
cat("Number of Phase II iterations = ",x$N2,"\n",sep="")
cat("\n")
cat("Computer time = ",paste(hrs,":",mins,":",secs,sep=""),"\n",sep="")

}

summary.ace<-function(object,...){

print.ace(x=object)}

################################ Utility Functions ###########################################################

################################################################################################################################
### 4.1 Linear model
################################################################################################################################

utillinmod<-function(d, B){
x<-cbind(1,d,d^2,d[,1]*d[,2])
log_deter<-as.vector(.Call("LMcpp", x, PACKAGE = "acebayes"))
log_deter+rnorm(B)}

optdeslinmod<-function(n, type = "ACE"){
if(type=="ACE"){
des<-linmodoptdesigns[linmodoptdesigns$n==n,]} else{
des<-truelinmodoptdesigns[truelinmodoptdesigns$n==n,]}
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

################################################################################################################################
### 4.2 Compartmental model
################################################################################################################################

### n=18 Gotwalt et al model Bayesian D-optimality

utilcomp18bad<-function(d, B){
low<-c(0.01884,0.298)
upp<-c(0.09884,8.298)
sam<-cbind(runif(n=B,min=low[1],max=upp[1]),runif(n=B,min=low[2],max=upp[2]))
as.vector(.Call("utilcomp18badcpp", d , sam, PACKAGE = "acebayes"))}

optdescomp18bad<-function(type = "ACE"){
if(type=="ACE"){
des<-comp18badoptdesign}
if(type=="Gotwalt"){
des<-comp18badoptdesign_got}
if(type=="Atkinson"){
des<-comp18badoptdesign_atk}
des<-as.matrix(des)
rownames(des)<-NULL
des}

### n=15 Ryan et al model SIG

utilcomp15bad<-function(d, B){
d2<-(d+1)*12
theta<-exp(cbind(rnorm(n=B,mean=log(0.1),sd=sqrt(0.05)),rnorm(n=B,mean=log(1),sd=sqrt(0.05)),rnorm(n=B,mean=log(20),sd=sqrt(0.05))))
as.vector(.Call("utilcomp15badcpp", d2 , theta, PACKAGE = "acebayes"))}

optdescomp15bad<-function(){
des<-comp15badoptdesign
des<-as.matrix(des)
rownames(des)<-NULL
des}

### n=15 Ryan et al model SIG

#inidescomp15sig<-function(rep){
#des<-comp15badinidesign$time[comp15badinidesign$rep==rep]
#des<-as.matrix(des)
#rownames(des)<-NULL
#des}

utilcomp15sig<-function(d, B){
D<-400
sigadd<-0.1
sigpro<-0.01
d2<-12*(as.vector(d)+1)
n1<-length(d2)
sam<-cbind(rnorm(n=2*B,mean=log(0.1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(20),sd=sqrt(0.05)))
sam<-exp(sam)
mu<-(D/matrix(rep(sam[,3],n1),ncol=n1))*(matrix(rep(sam[,2],n1),ncol=n1)/(matrix(rep(sam[,2],n1),ncol=n1)-matrix(rep(sam[,1],n1),ncol=n1)))*(exp(-matrix(rep(sam[,1],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE))-exp(-matrix(rep(sam[,2],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE)))
vv<-sigadd+sigpro*(mu^2)
y<-matrix(rnorm(n=n1*B,mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,])),ncol=n1)

frho<-as.vector(.Call("rowSumscpp", log(vv[-(1:B),]), PACKAGE = "acebayes"))
loglik<-as.vector(.Call("rowSumscpp", matrix(dnorm(x=as.vector(y),mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,]),log=TRUE),ncol=n1), PACKAGE = "acebayes"))

rsll4<-as.vector(.Call("utilcomp15sigcpp", y, mu[-(1:B),], vv[-(1:B),], frho, PACKAGE = "acebayes"))
MY3<-log(rsll4/B)

eval<-loglik-MY3

eval}

optdescomp15sig<-function(){
des<-comp15sigoptdesign
des<-as.matrix(des)
rownames(des)<-NULL
des}

### n=15 Ryan et al model SIG DRS

utilcomp15sigDRS<-function(d, B){
D<-400
sigadd<-0.1
sigpro<-0.01
d2<-24*qbeta(p=seq(from=0,to=1,length=17)[-c(1,17)],shape1=d[1,1],shape2=d[2,1])
n1<-length(d2)
sam<-cbind(rnorm(n=2*B,mean=log(0.1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(1),sd=sqrt(0.05)),rnorm(n=2*B,mean=log(20),sd=sqrt(0.05)))
sam<-exp(sam)
mu<-(D/matrix(rep(sam[,3],n1),ncol=n1))*(matrix(rep(sam[,2],n1),ncol=n1)/(matrix(rep(sam[,2],n1),ncol=n1)-matrix(rep(sam[,1],n1),ncol=n1)))*(exp(-matrix(rep(sam[,1],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE))-exp(-matrix(rep(sam[,2],n1),ncol=n1)*matrix(rep(d2,2*B),ncol=n1,byrow=TRUE)))
vv<-sigadd+sigpro*(mu^2)
y<-matrix(rnorm(n=n1*B,mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,])),ncol=n1)

frho<-as.vector(.Call("rowSumscpp", log(vv[-(1:B),]), PACKAGE = "acebayes"))
loglik<-as.vector(.Call("rowSumscpp", matrix(dnorm(x=as.vector(y),mean=as.vector(mu[1:B,]),sd=sqrt(vv[1:B,]),log=TRUE),ncol=n1), PACKAGE = "acebayes"))

rsll4<-as.vector(.Call("utilcomp15sigcpp", y, mu[-(1:B),], vv[-(1:B),], frho, PACKAGE = "acebayes"))
MY3<-log(rsll4/B)

eval<-loglik-MY3

eval}

optdescomp15sigDRS<-function(){
des<-comp15sigDRSoptdesign
des<-as.matrix(des)
rownames(des)<-NULL
des}

################################################################################################################################
### 4.3 Logistic Regression
################################################################################################################################

### Standard Logistic Regression - Bayesian D-optimality

utillrbad<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
x<-cbind(1,d)
eval<-as.vector(.Call("LRDcpp", x, beta, PACKAGE = "acebayes"))
eval}

optdeslrbad<-function(n, type = "ACE"){
if(type=="ACE"){
des<-LRBADoptdesign[LRBADoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]}
if(type=="Gotwalt1"){
des<-LRBADoptdesign2[LRBADoptdesign2$n==n & LRBADoptdesign2$type==0,]}
if(type=="Gotwalt2"){
des<-LRBADoptdesign2[LRBADoptdesign2$n==n & LRBADoptdesign2$type==1,]}
if(type=="Woods"){
des<-LRBADoptdesign2[LRBADoptdesign2$n==n & LRBADoptdesign2$type==2,]}
if(type!="ACE"){
des<-as.matrix(des)
des<-des[,-c(1,2)]}
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - Bayesian D-optimality

utilhlrbad<-function(d, B){
x<-cbind(1,d)
m<-6
n<-dim(x)[1]
G<-n/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
marker<-rep(1:G,each=m)
S<-t(bet*t(matrix(1-sqrt(runif(5*B)),ncol=5)))
gam<-matrix(0,ncol=G*5,nrow=B)
z<-matrix(0,nrow=n,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
gam[,5*(i-1)+j]<-runif(n=B,min=-S[,j],max=S[,j])}}

#eval<-as.vector(test.cpp(x=x,z=z,beta=beta,gam=gam,S=S))
eval<-as.vector(.Call("HLRDcpp", x, z, beta, gam, S, PACKAGE = "acebayes"))
eval}

optdeshlrbad<-function(n){
des<-HLRBADoptdesign[HLRBADoptdesign$n==n,] 
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Standard Logistic Regression - SIG

inideslrsig<-function(n, rep){
des<-LRSIGinidesign[LRSIGinidesign$n==n & LRSIGinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utillrsig<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)

sam<-matrix(runif(n=2*B*5),ncol=5)
for(jj in 1:5){
sam[,jj]<-sam[,jj]*(upp[jj]-low[jj])+low[jj]}

x<-cbind(1,d)
n1<-dim(x)[1]

rho<-1/(1+exp(-sam%*%t(x)))
Z<-log(1-rho[-(1:B),])
frho<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))

y<-matrix(rbinom(n=n1*B,size=1,prob=as.vector(rho[1:B,])),ncol=n1)
Z<-dbinom(x=y,size=1,prob=rho[1:B,],log=TRUE)                        ## loglik
rsll<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))
sam<-sam[-(1:B),]

rsll4<-as.vector(.Call("siglrcpp", y, x, sam, frho, PACKAGE = "acebayes"))

MY3<-log(rsll4/B)
eval<-rsll-MY3

eval}

optdeslrsig<-function(n){
des<-LRSIGoptdesign[LRSIGoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - SIG

inideshlrsig<-function(n, rep){
des<-LRSIGinidesign[LRSIGinidesign$n==n & LRSIGinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utilhlrsig<-function(d, B){
B2d<-B
B4d<-B
x<-cbind(1,d)
m<-6
n1<-dim(x)[1]
G<-n1/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
sam<-t(t(6*matrix(runif(n=5*(B2d+2*B4d)),ncol=5))+low)
marker<-rep(1:G,each=m)
zam<-matrix(0,ncol=G*5,nrow=B2d+2*B4d)
z<-matrix(0,nrow=n1,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
zam[,5*(i-1)+j]<-bet[j]*(1-sqrt(runif(B2d+2*B4d)))*(2*runif(B2d+2*B4d)-1)}}
etax<-sam%*%t(x)
etaz<-zam%*%t(z)
rho<-1/(1+exp(-etax-etaz))
y<-matrix(rbinom(n=B2d*n1,size=1,prob=as.vector(rho[1:B2d,])),ncol=n1)

frho<-as.vector(.Call("rowSumscpp", log(1-rho[((B2d+1):(B2d+B4d)),]), package="acebayes"))

rsll4<-.Call("sighlrcpp", cbind(x,z), y, cbind(sam[-(1:B2d),],zam[-(1:B2d),]), frho, sam[1:B2d,], exp(etax[1:B2d,]), exp(etaz[-(1:(B2d+B4d)),]), package="acebayes")

rsll4<-log(rsll4)

eval<-rsll4[,2]-rsll4[,1]

eval}

optdeshlrsig<-function(n){
des<-HLRSIGoptdesign[HLRSIGoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Standard Logistic Regression - Bayesian A-optimality

utillrbaa<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
x<-cbind(1,d)
eval<-as.vector(.Call("LRAcpp", x, beta, PACKAGE = "acebayes"))
eval}

optdeslrbaa<-function(n){
des<-LRBAAoptdesign[LRBAAoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - Bayesian A-optimality

utilhlrbaa<-function(d, B){
x<-cbind(1,d)
m<-6
n<-dim(x)[1]
G<-n/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
beta<-t(t(6*matrix(runif(n=5*B),ncol=5))+low)
marker<-rep(1:G,each=m)
S<-t(bet*t(matrix(1-sqrt(runif(5*B)),ncol=5)))
gam<-matrix(0,ncol=G*5,nrow=B)
z<-matrix(0,nrow=n,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
gam[,5*(i-1)+j]<-runif(n=B,min=-S[,j],max=S[,j])}}

#eval<-as.vector(test.cpp(x=x,z=z,beta=beta,gam=gam,S=S))
eval<-as.vector(.Call("HLRAcpp", x, z, beta, gam, S, PACKAGE = "acebayes"))
eval}

optdeshlrbaa<-function(n){
des<-HLRBAAoptdesign[HLRBAAoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Standard Logistic Regression - NSEL

inideslrnsel<-function(n, rep){
des<-LRNSELinidesign[LRNSELinidesign$n==n & LRNSELinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utillrnsel<-function(d, B){
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)

sam<-matrix(runif(n=2*B*5),ncol=5)
for(jj in 1:5){
sam[,jj]<-sam[,jj]*(upp[jj]-low[jj])+low[jj]}

x<-cbind(1,d)
n1<-dim(x)[1]

rho<-1/(1+exp(-sam%*%t(x)))
Z<-log(1-rho[-(1:B),])
frho<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))

y<-matrix(rbinom(n=n1*B,size=1,prob=as.vector(rho[1:B,])),ncol=n1)
fsam<-sam[1:B,]
sam<-sam[-(1:B),]

rsll4<-.Call("nsellrcpp", y, x, sam, frho, PACKAGE = "acebayes")

Z<-(rsll4-fsam)^2
eval<-as.vector(.Call("rowSumscpp", Z, PACKAGE = "acebayes"))

-eval}

optdeslrnsel<-function(n){
des<-LRNSELoptdesign[LRNSELoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

### Hierarchical Logistic Regression - NSEL

inideshlrnsel<-function(n, rep){
des<-LRSIGinidesign[LRSIGinidesign$n==n & LRSIGinidesign$rep==rep,-c(1,2)]
des<-as.matrix(des)
rownames(des)<-NULL
des}

utilhlrnsel<-function(d, B){
x<-cbind(1,d)
m<-6
n1<-dim(x)[1]
G<-n1/m
low<-c(-3,4,5,-6,-2.5)
upp<-c(3,10,11,0,3.5)
bet<-c(3,3,3,1,1)
sam<-t(t(6*matrix(runif(n=5*2*B),ncol=5))+low)
marker<-rep(1:G,each=m)
zam<-matrix(0,ncol=G*5,nrow=2*B)
z<-matrix(0,nrow=n1,ncol=G*5)
for(i in 1:G){
for(j in 1:5){
z[marker==i,5*(i-1)+j]<-x[marker==i,j]
zam[,5*(i-1)+j]<-bet[j]*(1-sqrt(runif(2*B)))*(2*runif(2*B)-1)}}
etax<-sam%*%t(x)
etaz<-zam%*%t(z)
rho<-1/(1+exp(-etax-etaz))
y<-matrix(rbinom(n=B*n1,size=1,prob=as.vector(rho[1:B,])),ncol=n1)
frho<-as.vector(.Call("rowSumscpp", log(1-rho[-(1:B),]), PACKAGE = "acebayes"))
rsll4<-.Call("nselhlrcpp", cbind(x,z), y, cbind(sam[-(1:B),],zam[-(1:B),]), frho, PACKAGE = "acebayes")

MY3<-(rsll4-sam[1:B,])^2
eval<-as.vector(.Call("rowSumscpp", MY3, PACKAGE = "acebayes"))

-eval
}

optdeshlrnsel<-function(n){
des<-HLRNSELoptdesign[HLRNSELoptdesign$n==n,]
des<-as.matrix(des)
des<-des[,-1]
rownames(des)<-NULL
des}

################################################################################################################################
### 4.4 Beetles
################################################################################################################################

utilbeetle<-function(d, B){
nd<-dim(d)[1]
L<-1.6907
U<-1.8839
R<-U-L
ck<-log(-log(0.5));lambda<-60
Xda<-cbind(1,(d+1)/2)
Xda[,2]<-Xda[,2]*R+L
Xda2<-cbind(Xda,Xda[,2]^2)
msam<-sample(x=1:6,size=B,prob=as.vector(probs),replace=TRUE)
tmsam<-c(length(msam[msam==1]),length(msam[msam==2]),length(msam[msam==3]),length(msam[msam==4]),length(msam[msam==5]),length(msam[msam==6]))
wr1<-sample(x=1:10000,size=tmsam[1],replace=TRUE)
wr2<-sample(x=1:10000,size=tmsam[2],replace=TRUE)
wr3<-sample(x=1:10000,size=tmsam[3],replace=TRUE)
wr4<-sample(x=1:10000,size=tmsam[4],replace=TRUE)
wr5<-sample(x=1:10000,size=tmsam[5],replace=TRUE)
wr6<-sample(x=1:10000,size=tmsam[6],replace=TRUE)
sam1<-as.matrix(lrlin[wr1,],ncol=2)
sam2<-as.matrix(lrquad[wr2,],ncol=3)
sam3<-as.matrix(colin[wr3,],ncol=2)
sam4<-as.matrix(coquad[wr4,],ncol=3)
sam5<-as.matrix(prlin[wr5,],ncol=2)
sam6<-as.matrix(prquad[wr6,],ncol=3)
eta1<-sam1%*%t(Xda)
eta2<-sam2%*%t(Xda2)
eta3<-sam3%*%t(Xda)
eta4<-sam4%*%t(Xda2)
eta5<-sam5%*%t(Xda)
eta6<-sam6%*%t(Xda2)
phi0<-c(-sam1[,1]/sam1[,2],
0.5*(-sam2[,2]+sqrt(sam2[,2]^2 - 4*sam2[,1]*sam2[,3]))/sam2[,3],
(ck-sam3[,1])/sam3[,2],
0.5*(-sam4[,2]+sqrt(sam4[,2]^2 - 4*(sam4[,1]-ck)*sam4[,3]))/sam4[,3],
-sam5[,1]/sam5[,2],
0.5*(-sam6[,2]+sqrt(sam6[,2]^2 - 4*sam6[,1]*sam6[,3]))/sam6[,3])
mu<-lambda*rbind(1/(1+exp(-eta1)),1/(1+exp(-eta2)),1-exp(-exp(eta3)),1-exp(-exp(eta4)),pnorm(eta5),pnorm(eta6))
y<-matrix(rpois(n=prod(dim(mu)),lambda=as.vector(mu)),ncol=nd)

wr1<-sample(x=1:10000,size=tmsam[1],replace=TRUE)
wr2<-sample(x=1:10000,size=tmsam[2],replace=TRUE)
wr3<-sample(x=1:10000,size=tmsam[3],replace=TRUE)
wr4<-sample(x=1:10000,size=tmsam[4],replace=TRUE)
wr5<-sample(x=1:10000,size=tmsam[5],replace=TRUE)
wr6<-sample(x=1:10000,size=tmsam[6],replace=TRUE)
sam1<-as.matrix(lrlin[wr1,],ncol=2)
sam2<-as.matrix(lrquad[wr2,],ncol=3)
sam3<-as.matrix(colin[wr3,],ncol=2)
sam4<-as.matrix(coquad[wr4,],ncol=3)
sam5<-as.matrix(prlin[wr5,],ncol=2)
sam6<-as.matrix(prquad[wr6,],ncol=3)
eta1<-sam1%*%t(Xda)
eta2<-sam2%*%t(Xda2)
eta3<-sam3%*%t(Xda)
eta4<-sam4%*%t(Xda2)
eta5<-sam5%*%t(Xda)
eta6<-sam6%*%t(Xda2)
mu1<-1/(1+exp(-eta1))
mu2<-1/(1+exp(-eta2))
mu3<-1-exp(-exp(eta3))
mu4<-1-exp(-exp(eta4))
mu5<-pnorm(eta5)
mu6<-pnorm(eta6)

phi<-c(-sam1[,1]/sam1[,2],
0.5*(-sam2[,2]+sqrt(sam2[,2]^2 - 4*sam2[,1]*sam2[,3]))/sam2[,3],
(ck-sam3[,1])/sam3[,2],
0.5*(-sam4[,2]+sqrt(sam4[,2]^2 - 4*(sam4[,1]-ck)*sam4[,3]))/sam4[,3],
-sam5[,1]/sam5[,2],
0.5*(-sam6[,2]+sqrt(sam6[,2]^2 - 4*sam6[,1]*sam6[,3]))/sam6[,3])

mu<-lambda*rbind(mu1,mu2,mu3,mu4,mu5,mu6)
Lmu<-log(mu)
fy<-lfactorial(y)

frho<--as.vector(.Call("rowSumscpp", mu, PACKAGE = "acebayes"))
ncr<--as.vector(.Call("rowSumscpp", fy, PACKAGE = "acebayes"))

rsll<-as.vector(.Call("beetlecpp", phi, y, Lmu, frho, ncr, PACKAGE = "acebayes"))

eval<-(phi0-rsll)^2

names(eval)<-NULL

-eval}

optdesbeetle<-function(n){
des<-beetleoptdesign[beetleoptdesign$n==n,]
des<-as.matrix(des)
des<-matrix(des[,-1],nrow=n)
rownames(des)<-NULL
des}



