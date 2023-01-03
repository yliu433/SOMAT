

### fixed_method: different statistics for testing the fixed effects
### acr_method: for all cont. traits random effects, use matrix form rather than cross product

SOMAT<-function(Y,X,G,W,bin_flag=NULL,fixed_method="Q",acr_method="MF"){

  n<-dim(Y)[1]
  K<-dim(Y)[2]
  d<-dim(X)[2]
  p<-dim(G)[2]
  qW<-dim(W)[1]
  
  if(is.null(bin_flag)){
    bin_flag<-c()
    for(k in 1:K){
      bin_flag[k]<-as.numeric(all(Y[,k]==0 | Y[,k]==1)) 
    }
  }
  
  #*************************************************************
  #
  #                 fixed effects / theta
  #
  #*************************************************************
  
  Utheta<-matrix(NA,qW,K) # each column for score statistic for trait k
  U_forV<-list() # each entry for n independent terms for trati k
  phi_tilde<-rep(1,K)
  
  Utheta_forV<-NULL
  
  # default values for cont. traits
  if(any(bin_flag==0)){
    F_theta_def <- t(X)%*%X
    F_theta_inv_def <- solve(F_theta_def)
    E_theta_def <- t(X)%*%G%*%t(W)
  }

  
  for (k in 1:K){
    y_work=Y[,k]
    
    if ( bin_flag[k]==1 ){#binary
      
      f1=glm(y_work~X-1, family="binomial")    #f1$fitted = exp(eta)/(1+exp(eta))
      
      alpha_hat=f1$coef
      
      eta1=X %*%alpha_hat   
      
      exp_eta=exp(eta1)
      
      fit_tilde=exp_eta/(1+exp_eta)
      
      RES =y_work - fit_tilde     
      
      wX<-sweep(X,1,as.vector(exp_eta/(1+exp_eta)^2),"*")
      
      F_theta= t(wX)%*%X          #equiv to B*
      
      F_theta_inv=solve(F_theta)   
      
      E_theta= t(wX)%*%(G%*%t(W))  #for theta
      
      phi_tilde[k]=1.0
      
      Utheta[,k]= ( W %*% t(G) %*% RES )/phi_tilde[k] # q by 1
      
      UfV <- W %*% t(G) - t(E_theta)%*%(F_theta_inv)%*%t(X)
      
      U_forV[[k]]<- sweep(UfV,2, RES,"*")/ phi_tilde[k] # q by n
      
    }else{#continuous
      
      f1=glm(y_work~X-1, family="gaussian")  
      
      alpha_hat=f1$coef
      
      eta1=X %*%alpha_hat
      
      RES = y_work - eta1 # residuals
      
      phi_tilde[k]=  sum( RES^2 ) / ( n  - d )
      
      Utheta[,k]= ( W %*% t(G) %*% RES )/phi_tilde[k] # q by 1
      
      F_theta_inv <- F_theta_inv_def  
      
      E_theta <- E_theta_def
      
      UfV <- W %*% t(G) - t(E_theta)%*%(F_theta_inv)%*%t(X)
      
      U_forV[[k]]<- sweep(UfV,2, RES,"*")/ phi_tilde[k] # q by n
      
    }
    
    Utheta_forV=rbind( Utheta_forV, U_forV[[k]] )
    
  }
  
  V_theta= (Utheta_forV)%*%t(Utheta_forV)
  
  Utheta_ts<-matrix(Utheta,ncol=1)  # test statistic
  
  
  if(fixed_method=="Q"){
    T_theta <- c(t(Utheta_ts)%*%solve(V_theta)%*% Utheta_ts)
    
    pval_theta= pchisq( T_theta, df=(K*qW), lower.tail=F)
  }else if(fixed_method=="Q1"){
    T_theta <- c(t(Utheta_ts)%*% Utheta_ts)
    theta_eigen<-Re(eigen(V_theta)$values)
    pval_theta=SKAT_davies(q=T_theta,lambda=theta_eigen)$Qq
  }else if(fixed_method=="T"){
    vars<-diag(V_theta)
    Z<-Utheta_ts/sqrt(vars)
    VZ_theta<-V_theta/outer(sqrt(vars),sqrt(vars))
    #VZ_theta<-cov2cor(V_theta)
    VZ_theta_inv<-solve(VZ_theta)
    T_theta<-(sum(VZ_theta_inv%*%Z)/sqrt(sum(VZ_theta_inv)))^2
    pval_theta <- pchisq( T_theta, df=(1), lower.tail=F)
  }else if(fixed_method=="T1"){
    vars<-diag(V_theta)
    Z<-Utheta_ts/(vars)
    VZ_theta<-V_theta/outer((vars),(vars))
    #VZ_theta<-cov2cor(V_theta)
    VZ_theta_inv<-solve(VZ_theta)
    T_theta<-(sum(VZ_theta_inv%*%Z)/sqrt(sum(VZ_theta_inv)))^2
    pval_theta <- pchisq( T_theta, df=(1), lower.tail=F)
  }else{
    print("please specify the correct methods for the fixed method")
  }

    
  #****************************************************************
  #
  #                random effects / tau
  #
  #****************************************************************
  
  Rtau<-matrix(NA,p,K) # each column for half of score statistic for trait k
  phi_hat<-rep(1,K)
  
  M=cbind(X, G%*%t(W))
  
  if(all(bin_flag==0) & acr_method=="MF"){  ### all continuous
    
    E_tau<-matrix(NA,n,K) # each column for residuals of trait k
    
    for (k in 1:K){
      y_work=Y[,k]
      
      f_tau=glm(y_work~M-1, family="gaussian")  
      
      eff_hat=f_tau$coef
      
      eta_for_tau= M %*%eff_hat 
      
      E_tau[,k]<-y_work - eta_for_tau
      
      phi_hat[k]=  sum( (y_work - eta_for_tau)^2) / ( n  - d - qW)
      
      Rtau[,k]= t(G)%*%(E_tau[,k])/ phi_hat[k]

      
    }
    
    # eigen by kronecker
    PM= diag(1,n) - M%*%solve(t(M)%*%M)%*%t(M)
    lam=eigen( t(G)%*%PM%*%G )$values
    
    ####
    Omega=cov(sweep(E_tau,2,(phi_hat),"/")) #called Omega* in the paper
    
    psai=eigen( Omega )$values
    
    tau_eigen= sort( kronecker(lam, psai),
                     decreasing=T)
    
  }else{ ### there exist binary traits
    
    U_forV<-list() # each entry for n independent terms for trati k
    Utau_forV<-NULL
    
    # default values for cont. traits
    if(any(bin_flag==0)){
      F_tau_def <- t(M)%*%M
      F_tau_inv_def <- solve(F_tau_def)
      E_tau_def <- t(M)%*%G
    }
    
    
    for (k in 1:K){
      y_work=Y[,k]
      
      if ( bin_flag[k]==1 ){#binary
        
        f1=glm(y_work~M-1, family="binomial")    #f1$fitted = exp(eta)/(1+exp(eta))
        
        alpha_hat=f1$coef
        
        eta1=M %*%alpha_hat   
        
        exp_eta=exp(eta1)
        
        fit_tau=exp_eta/(1+exp_eta)
        
        RES =y_work - fit_tau   
        
        wX<-sweep(M,1,as.vector(exp_eta/(1+exp_eta)^2),"*")
        
        F_tau= t(wX)%*%M          
        
        F_tau_inv=solve(F_tau)   
        
        E_tau = t(wX)%*%G  #for tau
        
        phi_hat[k]=1.0
        
        Rtau[,k]= ( t(G) %*% RES )/phi_hat[k] # p by 1
        
        UfV <- t(G) - t(E_tau)%*%(F_tau_inv)%*%t(M)
        
        U_forV[[k]]<- sweep(UfV,2, RES,"*")/ phi_hat[k] # p by n
        
      }else{#continuous
        
        f1=glm(y_work~M-1, family="gaussian")  
        
        alpha_hat=f1$coef
        
        eta1=M %*%alpha_hat
        
        RES = y_work - eta1 # residuals
        
        phi_hat[k]=  sum( RES^2 ) / ( n  - d -qW )
        
        Rtau[,k]= (t(G) %*% RES )/phi_hat[k] # p by 1
        
        F_tau_inv <- F_tau_inv_def  
        
        E_tau <- E_tau_def
        
        UfV <- t(G) - t(E_tau)%*%(F_tau_inv)%*%t(M)
        
        U_forV[[k]]<- sweep(UfV,2, RES,"*")/ phi_hat[k] # p by n
        
      }
      
      Utau_forV=rbind( Utau_forV, U_forV[[k]] )
      
    }
    
    V_tau= (Utau_forV)%*%t(Utau_forV)
    
    tau_eigen=Re(eigen(V_tau)$values)
    
  }
  
  Rtau_ts<-matrix(Rtau,ncol=1)
  
  T_tau<-c(t(Rtau_ts)%*%Rtau_ts)
  
  pval_tau=SKAT_davies(q=T_tau,lambda=tau_eigen)$Qq
    
   
  
  #################################################################
  #         overall pval by fisher
  #################################################################
  fisher=-2*( log(pval_theta) + log(pval_tau) )
  
  pval_fisher= pchisq( fisher, df=4, lower.tail=F)
    
  #################################################################
  #         overall pval by oSOMAT
  #################################################################
  pval_oSOMAT<-NA
  rho_opt1<-t_star<-NA
  
  if(fixed_method!="Q1"){
    ### observed fixed and random effects
    u_theta<-T_theta;u_tau2<-T_tau
    
    ### random effect eigenvalue function
    eigen_k<-function(k){sum(tau_eigen^k)}
    
    ### compute optimal rho for observed t_star
    
    if(fixed_method=="Q"){
      csdf <- qW*K
    }else{
      csdf <- 1
    }
      
    sum1<-(u_theta-csdf)*eigen_k(2)+(u_tau2-eigen_k(1))*csdf
    
    rho_opt10<-(u_tau2-eigen_k(1))*csdf/sum1
    
    c_k<-function(rho,k){
      (1-rho)^k*csdf+(rho)^k*eigen_k(k)
    }
    
    m1<-function(rho){
      ((1-rho)*u_theta+rho*u_tau2-c_k(rho,1))/sqrt(2*c_k(rho,2))
    }
    
    t1<-u_theta-csdf;t2<-u_tau2-eigen_k(1)
    
    if(t1>0 & t2>0){
      rho_opt1<-rho_opt10
      t_star<-m1(rho_opt1)
    }else{
      rho_list<-c(0,1)
      T_list<-m1(rho_list)
      t_star<-max(T_list)
      rho_opt1<-rho_list[which.max(T_list)]
    }
    
    #### end; rho_opt1 & t_star
    
    ### compute the p-value
    t_star1<-t_star*sqrt(2*csdf)+csdf
    t_star2<-t_star*sqrt(2*eigen_k(2))+eigen_k(1)
    if(t_star<=0){
      pval_oSOMAT<-1-(1-SKAT_liu(q=t_star2,lambda=tau_eigen))*pchisq(t_star1,df=csdf)
    }else{
      q1<-1-(1-SKAT_liu(q=t_star2,lambda=tau_eigen))*pchisq(csdf,df=csdf)
      delta<-function(u_theta_var){
        m_value<-sqrt(t_star^2-(u_theta_var-csdf)^2/(2*csdf))*sqrt(2*eigen_k(2))+eigen_k(1)
        Fq<-1-SKAT_liu(q=m_value,lambda=tau_eigen)
        f_inf<-Fq*dchisq(u_theta_var,df=csdf)
        return(f_inf)
      }
      q2<-integrate(delta,csdf,t_star1)$value
      pval_oSOMAT<-q1-q2
      
    }
    
  }
  
  return(list(pval_oSOMAT =pval_oSOMAT,pval_fSOMAT=pval_fisher,pval_theta=pval_theta,
              pval_tau=pval_tau,rho_opt=c(rho_opt1,t_star),bin_flag=bin_flag))
  
}
