hetsurr.fun = function(y1, y0, s1, s0, w1, w0, wf.grd = NULL, h0= NULL, h1= NULL, h2= NULL, h3= NULL, h4= NULL, var.want =FALSE,  type = "cont", test.want = FALSE, c.adj = 1) {
	n1 = length(y1)
	n0 = length(y0)
	if(test.want) {var.want = TRUE}
	if(type == "cont")	{
		if(is.null(wf.grd)){
			w.all = c(w1, w0)
			wf.grd = seq(quantile(w.all, 0.10), quantile(w.all, 0.90), length = 50)
		}
		if(is.null(h0) | is.null(h1) | is.null(h2) | is.null(h3) | is.null(h4)) {
			hsw1=c.adj*(bw.nrd(w1)*length(w1)^(1/5))*(length(w1)^(-0.4))   #h2
			hyw0 =c.adj*(bw.nrd(w0)*length(w0)^(1/5))*(length(w0)^(-0.4))  #h0
			hyw1 =c.adj*(bw.nrd(w1)*length(w1)^(1/5))*(length(w1)^(-0.4))  #h1  
			hsw1 =hyw0                                                     #h2. setting h2=h0
			hy2s1 =c.adj*(sqrt(4)*bw.nrd(s1)*length(s1)^(1/5))*(length(s1)^(-0.4)) #h3
			hy2w1 =c.adj*(sqrt(4)*bw.nrd(w1)*length(w1)^(1/5))*(length(w1)^(-0.4)) #h4
			hsw1=hyw0 #setting h2=h0
			hyw1=hy2w1 #setting h1=h4
		}
		else{
			hsw1=h2   #h2
			hyw0 =h0 #h0
			hyw1 =h1  #h1  
			hsw1 =hyw0   #h2. setting h2=h0
			hy2s1 =h3 #h3
			hy2w1 =h4 #h4
			hsw1=hyw0 #setting h2=h0
			hyw1=hy2w1 #setting h1=h4			
		}
		
		B.wf = length(wf.grd)
		dsw1=dw1=dw0=rep(0, n1)
		weights11=matrix(0, n1, n1)
		weightw11=matrix(0, n1, n1)
		weights10=matrix(0, n1, n0)
		weightw10=matrix(0, n1, n0)
		weightyw10=matrix(0, n1, n0)
		weightys10=matrix(0, n1, n0)
		
		for(i in 1:n1)
			{weights11[i,]=kf(s1[i]-s1, hy2s1)
  			weightw11[i,]=kf(w1[i]-w1, hy2w1)
  			weights10[i,]=kf(s1[i]-s0, hy2s1)
  			weightw10[i,]=kf(w1[i]-w0, hsw1)
  			weightyw10[i,]=kf(w1[i]-w0, hy2w1)
  			weightys10[i,]=kf(s1[i]-s0, hy2s1)
  			dsw1[i]=y1[i]-mu1swy(s1=s1,w1=w1,y1=y1,s=s1[i], w=w1[i], hy2w1=hy2w1, hy2s1=hy2s1)
  			dw1[i]=y1[i]-mu1wy(w1[i],hyw1=hyw1,y1=y1,w1=w1)}

		weightsw00=matrix(0, n0, n0)
		for(i in 1:n0)
  			{weightsw00[i,]=kf(w0[i]-w0, hsw1)}

		d0=dw0=rep(0, n0)
		for(i in 1:n0)
   			{dw0[i]=y0[i]-mu0wy(w0[i], hyw0=hyw0,y0=y0,w0=w0)
    		wt=weightys10*weightyw10[,i]
    		wt=t(t(wt)/apply(wt, 2, sum))
    		wt=wt%*%weightsw00[,i]/sum(weightsw00[,i])
    		d0[i]=mu1swy(s1=s1,w1=w1,y1=y1,s=s0[i], w=w0[i], hy2w1=hy2w1, hy2s1=hy2s1)-sum(y1*wt) }


		delta.grd=delta.s.grd=R.s.grd=rep(NA, B.wf)
		weightsw0=weightyw0=matrix(NA, B.wf, n0)
		weightyw1=weighty2w1=matrix(NA, B.wf, n1)

		r01=apply(weights10*weightw10, 1, mean)/apply(weights11*weightw11, 1, mean)

		eta20=eta10=matrix(NA, B.wf, n0)
		eta21=eta11=matrix(NA, B.wf, n1)
		mu1=mu0=rep(0, B.wf)
		for(bb in 1:B.wf)
   			{wf=wf.grd[bb]
    		mu1[bb]=mu1wy(wf,hyw1=hyw1,y1=y1,w1=w1)
    		mu0[bb]=mu0wy(wf,hyw0=hyw0,y0=y0,w0=w0)
    		delta.grd[bb]=mu1[bb]-mu0[bb]
    		delta.s.grd[bb]=mu10(s1=s1,w1=w1,y1=y1,w0=w0,s0=s0, w=wf, hy2w1=hy2w1, hy2s1=hy2s1, hsw1=hsw1)-mu0wy(wf,hyw0=hyw0,y0=y0,w0=w0)
    		R.s.grd[bb] = 1-delta.s.grd[bb]/delta.grd[bb]
    
        
    		if(var.want) {
    			weightsw0[bb,]=kf(wf-w0, hsw1)
    			weightyw0[bb,]=kf(wf-w0, hyw0)
    			weightyw1[bb,]=kf(wf-w1, hyw1)
    			weighty2w1[bb,]=kf(wf-w1, hy2w1)        

   				eta20[bb,]=weightsw0[bb,]*d0/mean(weightsw0[bb,])-weightyw0[bb,]*dw0/mean(weightyw0[bb,])
   				eta21[bb,]=weighty2w1[bb,]*dsw1/mean(weighty2w1[bb,])*r01 #10/16/21 changed denominator from mean(weightw11[bb,]) to mean(weighty2w1[bb,])

   				eta10[bb,]=-weightyw0[bb,]*dw0/mean(weightyw0[bb,])
   				eta11[bb,]=weightyw1[bb,]*dw1/mean(weightyw1[bb,])
   			}
   		}

		if(var.want) {
			BBB=500
			rstar=deltastar=deltastar.s=matrix(NA, B.wf, BBB)
			for(bbb in 1:BBB)
  				{xi0=rnorm(n0)
   				xi1=rnorm(n1)
   				for(bb in 1:B.wf)
      				{deltastar[bb,bbb]=mean(eta10[bb,]*xi0)+mean(eta11[bb,]*xi1)
       				deltastar.s[bb,bbb]=mean(eta20[bb,]*xi0)+mean(eta21[bb,]*xi1)
       				rstar[bb, bbb]=(deltastar.s[bb,bbb]-(1-R.s.grd[bb])*deltastar[bb,bbb])/delta.grd[bb]
       			}
   			}
 
			## the SE estimator for Rhat  
			sigma.r.grd=apply(rstar, 1, sd)
			## the cutoff value for the superemum statistic in constructing CB of R
			c0.r=quantile(apply(abs(rstar)/sigma.r.grd, 2, max), 0.95)

			## the SE estimator for delta
			sigma.delta.grd=apply(deltastar, 1, sd)
		
			## the SE estimator for delta.s
			sigma.delta.s.grd=apply(deltastar.s, 1, sd)
		}

		if(test.want) {
			### The omnibus test
			rhat=mean(delta.s.grd)/mean(delta.grd)
			dstar1=deltastar.s-rhat*deltastar
			dstar2=delta.grd%*%t((apply(deltastar.s, 2, mean)-rhat*apply(deltastar, 2, mean))/mean(delta.grd))
			dstar=dstar1-dstar2
			sigma.d.grd=apply(dstar, 1, sd)
			c0.d=quantile(apply(abs(dstar)/sigma.d.grd, 2, max), 0.95)
			t.star = apply(abs(dstar)/sigma.d.grd, 2, max)
			t.obs = max(abs(delta.s.grd-rhat*delta.grd)/sigma.d.grd)
			p.ob = mean(t.star >= t.obs)

			### The trend test
			sigma.cov=sd(apply(rstar*(wf.grd-mean(wf.grd)), 2, mean))
			z.trend=mean((R.s.grd-mean(R.s.grd))*(wf.grd-mean(wf.grd)))/sigma.cov
			p.trend = 2*(1-pnorm(abs(z.trend)))
		}

		if(!var.want & !test.want) {results = list("w.values" = wf.grd,"delta.w" = delta.grd, "delta.w.s" = delta.s.grd, "R.w.s" = R.s.grd)}
		if(var.want & !test.want) {results = list("w.values" = wf.grd,"delta.w" = delta.grd, "delta.w.s" = delta.s.grd, "R.w.s" = R.s.grd, "se.delta.w" = sigma.delta.grd, "se.delta.w.s" = sigma.delta.s.grd, "se.R.w.s" = sigma.r.grd, "conf.delta.w.lower" = delta.grd - 1.96*sigma.delta.grd, "conf.delta.w.upper" = delta.grd + 1.96*sigma.delta.grd, "conf.delta.w.s.lower" = delta.s.grd - 1.96*sigma.delta.s.grd, "conf.delta.w.s.upper" = delta.s.grd + 1.96*sigma.delta.s.grd,  "conf.R.w.s.lower" = R.s.grd - 1.96*sigma.r.grd, "conf.R.w.s.upper" = R.s.grd + 1.96*sigma.r.grd, "band.R.w.s.lower" = R.s.grd - c0.r*sigma.r.grd, "band.R.w.s.upper" = R.s.grd + c0.r*sigma.r.grd)}
		if(var.want & test.want) {results = list("w.values" = wf.grd,"delta.w" = delta.grd, "delta.w.s" = delta.s.grd, "R.w.s" = R.s.grd, "se.delta.w" = sigma.delta.grd, "se.delta.w.s" = sigma.delta.s.grd, "se.R.w.s" = sigma.r.grd, "conf.delta.w.lower" = delta.grd - 1.96*sigma.delta.grd, "conf.delta.w.upper" = delta.grd + 1.96*sigma.delta.grd, "conf.delta.w.s.lower" = delta.s.grd - 1.96*sigma.delta.s.grd, "conf.delta.w.s.upper" = delta.s.grd + 1.96*sigma.delta.s.grd,  "conf.R.w.s.lower" = R.s.grd - 1.96*sigma.r.grd, "conf.R.w.s.upper" = R.s.grd + 1.96*sigma.r.grd, "band.R.w.s.lower" = R.s.grd - c0.r*sigma.r.grd, "band.R.w.s.upper" = R.s.grd + c0.r*sigma.r.grd, "omnibus.test.statistic" = t.obs, "omnibus.p.value" = p.ob,"trend.test.statistic" = z.trend, "trend.p.value" = p.trend)}
		}
		
	
	if(type == "discrete") {
		if(is.null(wf.grd)){
			wf.grd = sort(unique(c(w1, w0)))
		}
		S.diff.w = vector(length = length(wf.grd))
    	for(jj in 1:length(wf.grd)) {
		S.diff.w[jj] = mean(y1[w1 == wf.grd[jj]]) - mean(y0[w0 == wf.grd[jj]])
    	}
		delta.s.w = delta.s.estimate.w(sone = s1, szero=s0, yone=y1, yzero=y0, wone=w1, wzero=w0, w.all = wf.grd) 
		R.s.w = 1-delta.s.w/S.diff.w 	
		if(var.want) {
			ll.1 = length(y1)
			ll.0 = length(y0)
			boot.num = 500
			num.w = length(wf.grd)
    		boot.mat = matrix(nrow = boot.num,ncol = length(wf.grd)*3)
				for(kk in 1:boot.num) {
					index.boot.1 = sample(1:ll.1, ll.1, replace = T) 
					index.boot.0 = sample(1:ll.0, ll.0, replace = T) 
					s1.b = s1[index.boot.1]
					y1.b = y1[index.boot.1]
					s0.b = s0[index.boot.0]
					y0.b = y0[index.boot.0]
					u1.b = w1[index.boot.1]
					u0.b = w0[index.boot.0]

					#			delta as a function of w
					S.diff.w.boot = vector(length = length(wf.grd))
    				for(jj in 1:length(wf.grd)) {
						S.diff.w.boot[jj] = mean(y1.b[u1.b == wf.grd[jj]]) - mean(y0.b[u0.b == wf.grd[jj]])
    				}
					boot.mat[kk,1:num.w] = S.diff.w.boot
					delta.s.w.boot = delta.s.estimate.w(sone = s1.b, szero=s0.b, yone=y1.b, yzero=y0.b, wone=u1.b, wzero=u0.b, w.all = wf.grd) 
 					boot.mat[kk,(num.w+1):(num.w+num.w)] = delta.s.w.boot
					boot.mat[kk,(num.w*2+1):(num.w*3)] = 1-delta.s.w.boot/S.diff.w.boot 	
		
				} #closes for loop
			var.results.mat = diag(var(boot.mat))
			se.delta.w = sqrt(var.results.mat[1:num.w])
			se.delta.w.s = sqrt(var.results.mat[(num.w+1):(num.w+num.w)])
			se.R.w.s = sqrt(var.results.mat[(num.w*2+1):(num.w*3)])
			#testing
			#make contrast matrix
			cont = diag(1,nrow = length(wf.grd)-1, ncol = length(wf.grd)) + cbind(rep(0,length(wf.grd)-1), diag(-1,nrow = length(wf.grd)-1, ncol = length(wf.grd)-1))
			vec.R = R.s.w
			delta.test = cont %*% as.matrix(vec.R)
			sand = solve(cont %*% var(boot.mat[,(length(wf.grd)*2+1):(length(wf.grd)*3)]) %*% t(cont))
			G = t(delta.test) %*% sand %*% delta.test
			test.stat = as.numeric(G)
			discrete.p.value = as.numeric(1-pchisq(G, num.w-1))
		}
		if(!var.want & !test.want) {results = list("w.values" = wf.grd, "delta.w" = S.diff.w, "delta.w.s" = delta.s.w, "R.w.s" = R.s.w)}
		if(var.want & !test.want) {results = list("w.values" = wf.grd, "delta.w" = S.diff.w, "delta.w.s" = delta.s.w, "R.w.s" = R.s.w, "se.delta.w" = se.delta.w, "se.delta.w.s" = se.delta.w.s, "se.R.w.s" = se.R.w.s, "conf.delta.w.lower" = S.diff.w - 1.96*se.delta.w, "conf.delta.w.upper" = S.diff.w + 1.96*se.delta.w, "conf.delta.w.s.lower" = delta.s.w - 1.96*se.delta.w.s, "conf.delta.w.s.upper" = delta.s.w + 1.96*se.delta.w.s,  "conf.R.w.s.lower" = R.s.w - 1.96*se.R.w.s, "conf.R.w.s.upper" = R.s.w + 1.96*se.R.w.s)}
		if(var.want & test.want) {results = list("w.values" = wf.grd, "delta.w" = S.diff.w, "delta.w.s" = delta.s.w, "R.w.s" = R.s.w, "se.delta.w" = se.delta.w, "se.delta.w.s" = se.delta.w.s, "se.R.w.s" = se.R.w.s, "conf.delta.w.lower" = S.diff.w - 1.96*se.delta.w, "conf.delta.w.upper" = S.diff.w + 1.96*se.delta.w, "conf.delta.w.s.lower" = delta.s.w - 1.96*se.delta.w.s, "conf.delta.w.s.upper" = delta.s.w + 1.96*se.delta.w.s,  "conf.R.w.s.lower" = R.s.w - 1.96*se.R.w.s, "conf.R.w.s.upper" = R.s.w + 1.96*se.R.w.s, "test.statistic" = test.stat, "p.value" = discrete.p.value)}
	} #close discrete
	return(results)
}


hetsurr.plot = function(hetsurr.results, xlab.name = "Baseline Covariate", placement = "bottomleft") {
	oldpar <- par(no.readonly = TRUE)
	on.exit(par(oldpar)) 
	wf.grd = hetsurr.results$w.values
	par(mfrow = c(3,1))
	plot(wf.grd, hetsurr.results$delta.w, typ = "l",lwd=2, ylab = expression(Delta), xlab = xlab.name, ylim = c(min(hetsurr.results$conf.delta.w.lower),max(hetsurr.results$conf.delta.w.upper)))
	points(wf.grd, hetsurr.results$conf.delta.w.lower, lwd=2,lty = 2, typ = "l")
	points(wf.grd, hetsurr.results$conf.delta.w.upper, lwd = 2, lty = 2, typ = "l")
	legend(placement,c("Estimate","Pointwise Confidence Interval"), lty = c(1,2), col=c("black", "black"), lwd = 2)

	plot(wf.grd, hetsurr.results$delta.w.s, typ = "l",lwd=2,  ylab = expression(Delta[S]), xlab = xlab.name, ylim = c(min(hetsurr.results$conf.delta.w.s.lower),max(hetsurr.results$conf.delta.w.s.upper)))
	points(wf.grd, hetsurr.results$conf.delta.w.s.lower, lwd=2,lty = 2, typ = "l")
	points(wf.grd, hetsurr.results$conf.delta.w.s.upper, lwd = 2, lty = 2, typ = "l")
	legend(placement,c("Estimate","Pointwise Confidence Interval"), lty = c(1,2), col=c("black", "black"), lwd = 2)

	plot(wf.grd, hetsurr.results$R.w.s, typ = "l",lwd=2, ylab = expression(R[S]), xlab = xlab.name, ylim=c(min(hetsurr.results$band.R.w.s.lower),max(hetsurr.results$band.R.w.s.upper)))
	points(wf.grd, hetsurr.results$conf.R.w.s.lower, lwd=2,lty = 2, typ = "l")
	points(wf.grd, hetsurr.results$conf.R.w.s.upper, lwd = 2, lty = 2, typ = "l")
	points(wf.grd, hetsurr.results$band.R.w.s.lower, lwd=2,lty = 3, typ = "l")
	points(wf.grd, hetsurr.results$band.R.w.s.upper, lwd = 2, lty = 3, typ = "l")
	legend(placement, c("Estimate","Pointwise Confidence Interval","Confidence Band"), lty = c(1,2,3), col=c("black", "black", "black"), lwd = 2)

	
}


delta.s.estimate.w=function(sone, szero, yone, yzero, wone, wzero, w.all, weight.perturb = NULL) 
{
    	delta.s = vector(length = length(w.all))
    	for(jj in 1:length(w.all)) {
			R.s = R.s.estimate(sone=sone[wone == w.all[jj]], szero=szero[wzero == w.all[jj]], yone = yone[wone == w.all[jj]], yzero = yzero[wzero == w.all[jj]], extrapolate = TRUE)
			delta.s[jj] = R.s$delta.s
    	}
    return(delta.s)
}

## kernel function
kf=function(x, h){return(dnorm(x/h)/h)}

## estimate mu1(s, w)=E(Y1|S1=s, W1=w)
mu1swy=function(s1,w1,y1,s, w, hy2w1, hy2s1)
{m=length(s)
 res=rep(NA, m)
 for(i in 1:m)
   {weight=kf(s[i]-s1, hy2s1)*kf(w-w1, hy2w1)
    res[i]=sum(weight*y1)/sum(weight)
    }
  if(sum(is.na(res))!=0) {
  	ind = which(is.na(res))
  	for(i in ind) 
   		{mat = cbind(abs(s-s[i]), res) 
   		 mmm = which(mat[,1] == min(mat[-i,1]))[1]   
   		 res[i] = res[mmm]  
   		}}
return(res)
}


## estimate mu10(w)=\int mu1(s, w)dF0(s|w)
mu10=function(s1,w1,y1,w0,s0, w, hy2w1, hy2s1, hsw1)
    {weight0=kf(w0-w, hsw1)
     return(sum(mu1swy(s1=s1,w1=w1,y1=y1,s=s0, w=w, hy2w1=hy2w1, hy2s1=hy2s1)*weight0)/sum(weight0))
    }

## estimate mu1(w)=E(Y1|W1=w)
mu1wy=function(w,hyw1,y1,w1)
    {weight=kf(w-w1, hyw1)
     return(sum(weight*y1)/sum(weight))
    }


## estimate mu0(w)=E(Y0|W0=w)
mu0wy=function(w,hyw0,y0,w0)
    {weight=kf(w-w0, hyw0)
     return(sum(weight*y0)/sum(weight))
    }
