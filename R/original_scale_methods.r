
mpfrthr <- 9; mpfrpbn <- 10 # mpfr - high  precision kick in levels for alternating series
mat.sort.c <- function(mat,mat.sz,n) {mat[rank(mat[,n],ties.method="first"),] <- mat[1:mat.sz,]; return(mat)} # Define matrix sort (ranking along column)
mat.sort.r <- function(mat,mat.sz,n) {mat[,rank(mat[n,],ties.method="first")] <- mat[,1:mat.sz]; return(mat)} # Define matrix sort (ranking along a row)
ins.sort.c <- function(vec,sor.mat,sor.mat.sz,n){ # Insert a vector into a sorted matrix (ranking along a column)
    r.min <- 1; r.max <- sor.mat.sz; tau <- vec[n]
    if(tau > sor.mat[r.max,n]){r.min <- r.max <- sor.mat.sz+1}else{
        while(r.max > r.min){
            r.mid <- floor(0.5*(r.max+r.min))
            if(tau <= sor.mat[r.mid,n]){r.max <- r.mid}else{r.min <- r.mid+1}}}
    return(rbind(sor.mat[seq_len(r.min-1),,drop=FALSE],vec,sor.mat[seq_len(sor.mat.sz-r.max+1)+r.max-1,,drop=FALSE]))}

ins.sort.r <- function(vec,sor.mat,sor.mat.sz,n){ # Insert a vector into a sorted matrix (ranking along a row)
    r.min <- 1; r.max <- sor.mat.sz; tau <- vec[n]
    if(tau > sor.mat[n,r.max]){r.min <- r.max <- sor.mat.sz+1}else{
        while(r.max > r.min){
            r.mid <- floor(0.5*(r.max+r.min))
            if(tau <= sor.mat[n,r.mid]){r.max <- r.mid}else{r.min <- r.mid+1}}}
    return(cbind(sor.mat[,seq_len(r.min-1),drop=FALSE],vec,sor.mat[,seq_len(sor.mat.sz-r.max+1)+r.max-1,drop=FALSE]))}



#############################################
##### 2.1.1 - Bessel Functionals - See Pollock et al. 2015
#############################################

eazeta        <- function(n,s,t,x,y,L,U){if(max(x-U,y-U,L-x,L-y)>=0){1}else{j<-1:(ceiling(n/2));P<--2/(t-s);D<-U-L;D1<-D*j+L;D2<-D*j-U;z<-y-x;if(gtools::even(n)){sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y))-exp(P*j^2*D^2-P*j*D*z)-exp(P*j^2*D^2+P*j*D*z))}else{sum(exp(P*(D1-x)*(D1-y))+exp(P*(D2+x)*(D2+y)))-sum(exp(P*j[1:length(j)-1]^2*D^2-P*j[1:length(j)-1]*D*z)+exp(P*j[1:length(j)-1]^2*D^2+P*j[1:length(j)-1]*D*z))}}} # Zeta Functional
eagamma        <- function(n,s,t,x,y,L,U){1-eazeta(n,s,t,x,y,L,U)} # Gamma Functional
eapsi        <- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j-(xoy-m))*exp(P*(abs(u-m)*j-(xoy-m)))} # Psi Functional
eachi        <- function(j,s,t,m,xoy,u){P<--2*abs(u-m)*j/(t-s);(2*abs(u-m)*j+(xoy-m))*exp(P*(abs(u-m)*j+(xoy-m)))} # Chi Functional
eadel2      <- function(n,s,t,m,xoy,u){if(max(xoy-u,m-xoy)>=0){0}else{if(gtools::even(n)){j<-1:(n/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)}else{if(n>1){j<-1:((n-1)/2);1-(sum(eapsi(j,s,t,m,xoy,u)-eachi(j,s,t,m,xoy,u)))/(xoy-m)-eapsi(max(j)+1,s,t,m,xoy,u)/(xoy-m)}else{1-eapsi(1,s,t,m,xoy,u)/(xoy-m)}}}}
eadelR        <- function(n,s,t,x,y,m,u){if(x==m){xI<-1}else{xI<-0};if(y==m){yI<-1}else{yI<-0};if(max(xI,yI)==1){delT<-2}else{delT<-1};if(m>min(x,y)){x<--x;y<--y;m<--m;u<--u};if(max(x-u,y-u,m-x,m-y)>=0){out<-0};if(delT==1){out<-eagamma(n,s,t,x,y,m,u)/(1-exp(-2*(x-m)*(y-m)/(t-s)))};if(delT==2){if(xI*yI==0){xoy<-max(x,y);out<-eadel2(n,s,t,m,xoy,u)}else{out<-0}};if(out<0){out<-0};if(out>1){out<-1};if((t-s)==0){out <- 1}; out}
eadelC         <- function(mt,s,t,x,y,m,u){if(mt>=mpfrthr){pbn<-mt*mpfrpbn;s<-mpfr(s,precBits=pbn);t<-mpfr(t,precBits=pbn);x<-mpfr(x,precBits=pbn);y<-mpfr(y,precBits=pbn);m<-mpfr(m,precBits=pbn);u<-mpfr(u,precBits=pbn);c(s1=eadelR(mt,s,t,x,y,m,u),s2=eadelR(mt+1,s,t,x,y,m,u))}else{eadel_pair_cpp(mt,s,t,x,y,m,u)}}

#############################################
#### 2.1.2 - Simulate Bessel 3 Process Mid Points along with Bessel Exceedance Evaluation - See Pollock (PhD Thesis, 2013)
#############################################

scale.midpt <- function(q, s, tau, x, y, minI, bdry){

    ## Repeat until acceptance
    repeat{
        ### Simulation of Bessel proposal
        c.scaling     <- c(tau-s,minI*(x-y)*(tau-q)/(tau-s)^1.5,sqrt((tau-q)*(q-s))/(tau-s))
        bb.sims          <- rnorm(3,0,sd=c.scaling[3])
        w                 <- y + minI*sqrt(c.scaling[1]*((c.scaling[2]+bb.sims[1])^2+bb.sims[2]^2+bb.sims[3]^2))

        ### Determine acceptance / rejection of proposal
        u               <- runif(1,0,1) # Uniform RV to make accept / reject decision
        counter     <- ceiling(sqrt((tau-s)+(bdry-y)^2)/(2*abs(bdry-y))); if(gtools::even(counter)==1){counter <- counter + 1} # Determining the minimum threshold for computing the boundary
        repeat{
            bounds     <- eadelC(counter,s,q,x,w,y,bdry)*eadelC(counter,q,tau,w,y,y,bdry) # Determine crossing probability
            if(u <= bounds[1]){accI <- 1; break} # Determine whether to accept
            if(u > bounds[2]){accI <- 0; break} # Determine whether to reject
            counter <- counter + 2} # Increase counter

        ### If accept then break loop
        if(accI == 1){break}}
    ### Output midpt
    list(w=w)}

#############################################
##### 2.2.2 - Proposal and rejection sampler
#############################################

dev.pr        <- function(Jst.t=0.64,Jst.rat=0.5776972){ # First Hitting Time (J*) Proposal Function
    U <- runif(1,0,1) # Simulate uniform
    if(U < Jst.rat){ # Case 1, U <- p/(p+q)
        E <- rexp(1,1); X <- Jst.t + 8*E/(pi^2)
    }else{ # Case 2, U >= p/(p+q)
        repeat{E <- rexp(2,1); if(E[1]^2 <= 2*E[2]/Jst.t){X <- Jst.t/(1+Jst.t*E[1])^2; break}}}
    list(X=X,U=U,E=E)}

dev.rej         <- function(X,Jst.t=0.64){ # First Hitting Time (J*) Rejection Sampler given a proposed point (X)
    if(X <= Jst.t){ # Series sampler 1, used around critical point 0.64 (t*<= 0.64, see Appendix C)
        S <- exp(-1/(2*X))/2; n <- -1; Y <- runif(1,0,1)*S # Initialise alternating sequence and simulate appropriate Uniform
        repeat{
            n <- n+2 # Indexation of sequence for subsequent odd / gtools::even terms
            S <- S - (n+0.5)*exp(-2*(n+0.5)^2/X); if(Y <= S){Acc <- 1; n <- n-1; break} # Odd n
            S <- S + (n+1.5)*exp(-2*(n+1.5)^2/X); if(Y >= S){Acc <- 0; break}} # gtools::even n
    }else{ # Series sampler 2, used around critical point 0.64 (t*> 0.64, see Appendix C)
        S <- exp(-pi^2*X/8)/2; n <- -1; Y <- runif(1,0,1)*S # Initialise alternating sequence and simulate appropriate Uniform
        repeat{
            n <- n+2 # Indexation of sequence for subsequent odd / gtools::even terms
            S <- S - (n+0.5)*exp(-(n+0.5)^2*pi^2*X/2); if(Y <= S){Acc <- 1; n <- n-1; break} # Odd n
            S <- S + (n+1.5)*exp(-(n+1.5)^2*pi^2*X/2); if(Y >= S){Acc <- 0; break}}} # gtools::even n
    list(Acc=Acc,S=S,n=n,Y=Y,X=X)} # Output A/R. Note S will not match target density as constants have been removed (see Appendix C)

#############################################
##### 2.2.3 - First Hitting Time (J*) Simulation
#############################################

bm.pass        <- function(s=0,x=0,theta=1,Jst.t=0.64,Jst.rat=0.5776972){ # theta denotes passage level
    repeat{sim <- dev.rej(dev.pr()$X); if(sim$Acc==1){break}}
    tau <- s + theta^2*sim$X; minI <- 2*rbinom(1,1,0.5)-1; y <- x - theta*minI
    list(tau=tau,y=y,minI=minI)}


#############################################
#### 1.2 - Resampling Algorithms
#############################################
##### 1.2.1 - Multinomial Resampling
#############################################

multi.resamp     <- function(p.wei,n=length(p.wei)){ # Multinomial Resampling
  if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
    p.idx        <- sample(1:length(p.wei),n,replace=TRUE,prob=p.wei) # Sampled Index
  }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
  list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.2 - Systematic Resampling
#############################################

system.resamp    <- function(p.wei,n=length(p.wei)){ # Systematic Resampling
  if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
    cum.wei        <- cumsum(p.wei)/sum(p.wei) # Normalise and calculate cumulative weights
    samp.vec    <- seq(runif(1,0,1/n),1,1/n) # Vector of cdf samples
    p.idx        <- numeric(0); for(i in 1:n){p.idx <- c(p.idx,length(cum.wei[samp.vec[i]>=cum.wei])+1)} # Sample from cdf
  }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
  list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.3 - Stratified Resampling
#############################################

strat.resamp    <- function(p.wei,n=length(p.wei)){ # Stratified Resampling
  vec.length  <- length(p.wei) # Calculate the length of the input vector
  cum.wei        <- cumsum(p.wei) # Calculate cumulative weights
  if(cum.wei[vec.length]>0){ # Check whether resampling possible
    cum.wei     <- cum.wei/cum.wei[vec.length]  # Normalise cumulative weights
    samp.vec    <- seq(0,(n-1)/n,1/n) + runif(n,0,1/n) # Vector of cdf samples
    p.idx <- findInterval(samp.vec,cum.wei)+1 # Sample from cdf
  }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
  list(p.idx=p.idx)} # Return particle indices

#############################################
##### 1.2.4 - Residual Resampling
#############################################

resid.resamp    <- function(p.wei,n=length(p.wei),nest.resamp=strat.resamp){ # Residual Resampling
  if((sum(p.wei)>0)&(n>0)){ # Check whether resampling possible
    det.resamp    <- floor(n*p.wei/sum(p.wei))
    p.idx        <- c(rep(1:length(p.wei),det.resamp),nest.resamp(p.wei-det.resamp*(sum(p.wei)/n),n=n-sum(det.resamp))$p.idx)
  }else{p.idx     <- numeric(0)} # If resampling not possible return empty vector
  list(p.idx=p.idx)} # Return particle indices
