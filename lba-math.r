fptcdf=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) return(pnorm(chi/z,mean=driftrate,sd=sddrift,lower.tail=F))
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu ; xx=chiminuszu-x0max
  chizu=chiminuszu/zs ; chizumax=xx/zs
  tmp1=zs*(dnorm(chizumax)-dnorm(chizu))
  tmp2=xx*pnorm(chizumax)-chiminuszu*pnorm(chizu)
  1+(tmp1+tmp2)/x0max
}

fptpdf=function(z,x0max,chi,driftrate,sddrift) {
  if (x0max<1e-10) return( (chi/z^2)*dnorm(chi/z,mean=driftrate,sd=sddrift))
  zs=z*sddrift ; zu=z*driftrate ; chiminuszu=chi-zu
  chizu=chiminuszu/zs ; chizumax=(chiminuszu-x0max)/zs
  (driftrate*(pnorm(chizu)-pnorm(chizumax)) + 
    sddrift*(dnorm(chizumax)-dnorm(chizu)))/x0max
}

rlba=function(n,b,A,vs,s,t0,st0=0,truncdrifts=TRUE){
  n.with.extras=ceiling(n*(1+3*prod(pnorm(-vs))))
  drifts=matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE)
  if (truncdrifts) {
    repeat {
      drifts=rbind(drifts,matrix(rnorm(mean=vs,sd=s,n=n.with.extras*length(vs)),ncol=length(vs),byrow=TRUE))
      tmp=apply(drifts,1,function(x) any(x>0))
      drifts=drifts[tmp,]
      if (nrow(drifts)>=n) break
    }
  }
  drifts=drifts[1:n,]
  drifts[drifts<0]=0
  starts=matrix(runif(min=0,max=A,n=n*length(vs)),ncol=length(vs),byrow=TRUE)
  ttf=t((b-t(starts)))/drifts
  rt=apply(ttf,1,min)+t0+runif(min=-st0/2,max=+st0/2,n=n)
  resp=apply(ttf,1,which.min)
  list(rt=rt,resp=resp)
}


                        

dlba=function(t,t0,x0max,chi,v1,v0,sdI,truncdrifts=TRUE) {
  # Generates defective PDF for responses on node #1.
  # "truncdrifts" sets whether that part of the multi-variate
  # normal distribution on drift rates which would otherwise
  # lead to non-terminating trials is truncated.
  if(t < 0){
    drift = c(v0,v1)
  }else{
    drift = c(v1,v0)
  }
  t = abs(t)
  t = t - t0
  x0max = rep(x0max,2)
  chi = rep(chi,2)
  
  sdI = rep(drift,2)
  
  N=length(drift) # Number of responses.
  if (N>2) {
    tmp=array(dim=c(length(t),N-1))
    for (i in 2:N) tmp[,i-1]=fptcdf(z=t,x0max=x0max[i],chi=chi[i],
      driftrate=drift[i],sddrift=sdI[i])
    G=apply(1-tmp,1,prod)
  } else {
    G=1-fptcdf(z=t,x0max=x0max[2],chi=chi[2],driftrate=drift[2],sddrift=sdI[2])
  }
  out=G*fptpdf(z=t,x0max=x0max[1],chi=chi[1],driftrate=drift[1],sddrift=sdI[1])
  if (truncdrifts) {
    out=out/(1-prod(pnorm(-drift/sdI)))
    out[t<=0]=0
    out = pmax(out,1e-10)
    return(out)
  } else {
    out = pmax(out,1e-10)
    return(out)
  }
}

n1PDF=function(t,x0max,chi,drift,sdI,st0=0,truncdrifts=TRUE) {
  N=length(drift) # Number of responses
  if (length(x0max)<N) x0max=rep(x0max,length.out=N)
  if (length(chi)<N) chi=rep(chi,length.out=N)
  if (length(sdI)<N) sdI=rep(sdI,length.out=N)
  if (length(st0)>1) st0=st0[1] # Only ONE non-decision time.
  if (st0==0) return(n1PDFfixedt0(t,x0max,chi,drift,sdI,truncdrifts=truncdrifts))
  tmpf=function(t,x0max,chi,drift,sdI,st0,truncdrifts)
    n1PDFfixedt0(t,x0max,chi,drift,sdI,truncdrifts=truncdrifts)/st0
  outs=numeric(length(t))
  for (i in 1:length(outs)) outs[i]=integrate(f=tmpf,lower=t[i]-st0/2,upper=t[i]+st0/2,
    x0max=x0max,chi=chi,drift=drift,sdI=sdI,st0=st0,truncdrifts=truncdrifts)$value
  outs
}


n1CDF=function(t,x0max,chi,drift,sdI,st0=0,truncdrifts=TRUE) {
  # Generates defective CDF for responses on node #1. 
  N=length(drift) # Number of responses
  if (length(x0max)<N) x0max=rep(x0max,length.out=N)
  if (length(chi)<N) chi=rep(chi,length.out=N)
  if (length(sdI)<N) sdI=rep(sdI,length.out=N)
  if (length(st0)>1) stop("Only one value of st0 allowed.")
  if (st0<1e-6) st0=0 # Integral can fail for small st0.
  outs=numeric(length(t)) ; bounds=c(-st0/2,t)
  for (i in 1:length(t)) {
    tmp="error"
    repeat {
      if (bounds[i]>=bounds[i+1]) {outs[i]=0;break}
      tmp=try(integrate(f=n1PDF,lower=bounds[i],upper=bounds[i+1],subdivisions=1000,
        x0max=x0max,chi=chi,drift=drift,sdI=sdI,st0=st0,truncdrifts=truncdrifts)$value,silent=T)
      if (is.numeric(tmp)) {outs[i]=tmp;break}
      # Try smart lower bound.
      if (bounds[i]<=0) {
	bounds[i]=max(c((chi-0.98*x0max)/(max(mean(drift),drift[1])+2*sdI)[1],0))
	next
      }
      # Try smart upper bound.
      if (bounds[i+1]==Inf) {
	bounds[i+1]=0.02*max(chi)/(mean(drift)-2*mean(sdI))
	next
      }
      stop("Error in n1CDF that I could not catch.")
    }
  }
  cumsum(outs)
}

n1mean=function(x0max,chi,drift,sdI,truncdrifts=TRUE) {
  # Generates mean RT for responses on node #1. 
   pc=n1CDF(Inf,x0max,chi,drift,sdI,truncdrifts=truncdrifts)
   fn=function(t,x0max,chi,drift,sdI,st0=0,truncdrifts,pc)
     t*n1PDF(t,x0max,chi,drift,sdI,st0,truncdrifts=truncdrifts)/pc
   tmp=integrate(f=fn,lower=0,upper=100*chi,x0max=x0max,chi=chi,pc=pc,
     drift=drift,sdI=sdI,st0=st0,truncdrifts=truncdrifts)$value
   list(mean=tmp,p=pc)
}

