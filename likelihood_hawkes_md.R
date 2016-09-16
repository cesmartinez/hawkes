# negative lambda0 forbidden
# negative alpha forbidden

likelihood_hawkes_md = function(brtscomplete,alpha,beta,lambda0,TS,cond='stem'){
      
    alpha   = round(alpha,digits=14)
    beta    = round(beta,digits=14)
    lambda0 = round(lambda0,digits=14)
  
    if(all(lambda0>=0)==F){return(-Inf)}
    if(all(alpha>=0)==F){return(-Inf)}
    if(all(beta>=0)==F){return(-Inf)}
    if((cond %in% c('crown','stem'))==F){return( 'cond should be \'stem\' or \'crown\' ')}
    brts    = brtscomplete[2,]
    process = brtscomplete[1,]
    
    di = length(lambda0)
    alpha = matrix(alpha,di,di)
    beta  = matrix(beta,di,di)
    nbtimes = length(brts)
    

    #durations between the branching times and the end of the process
    brtsdiff = TS - brts
    brtsdiff = matrix(rep(brtsdiff,di),nrow=di,ncol=nbtimes,byrow=T)

   
    pr=process[1]
    
    #variable part of the rates
    ratesvar = matrix(rep(0,di^2),di,di)
    ratesvar[,pr] = alpha[,pr]
    
    ##initialization for the log likelihood
    if(cond == 'stem'){ll = log(lambda0[pr])}
    if(cond=='crown'){ll=0}
    
    ############################################################
    
    #loop that calculates ratesvar at each branchingtimes      
    for(i in 2:nbtimes){
    betatime = exp(- beta * (brts[i] - brts[i-1]))
    
    ratesvar = ratesvar * betatime
    pr=process[i]
    
    #we add the log of the rate of the process that has just jumped
    ll = ll + log(lambda0[pr] + sum(ratesvar[pr,])) 
    
    #update of the rates
    ratesvar[,pr] = ratesvar[,pr] + alpha[,pr]
    
    # #check that the rates are positive (we allow negative lambda)
    # if(ratesvar + matrix(lambda0)
    
    
    }
    
    ####################################################################################
    
    #ratesvar integrated intensity between branching times and the end of the process
    intratesvar = matrix(rep(0,nbtimes),nrow=di, ncol=nbtimes)
    
    for(j in 1:di){
      ind1 = which(process==j)
      #check that there is at least one event of this type
      if(length(ind1>0)){
        #case where beta>0
        ind2 = which(beta[,j]==0)
        if(length(ind2)>0){
        intratesvar[ind2,ind1] = alpha[ind2,j] * brtsdiff[ind2,ind1]
        intratesvar[-ind2,ind1] = (alpha[-ind2,j]/beta[-ind2,j]) * (1 - exp(-beta[-ind2,j] * brtsdiff[ind2,ind1]))  
        }else{
          intratesvar[,ind1] = (alpha[,j]/beta[,j]) * (1 - exp(-beta[,j] * brtsdiff[,ind1]))}
      }
    }
    

    return(ll-sum(lambda0 * TS) -sum(intratesvar))
    
}



likoptim_md = function(par,brtscomplete,TC,cond='stem'){
  ll=-likelihood_hawkes_md(brtscomplete=brtscomplete, alpha=par[2], beta=par[3],cond=cond, lambda0=par[1], TS=TC)
  return(ll)
}



