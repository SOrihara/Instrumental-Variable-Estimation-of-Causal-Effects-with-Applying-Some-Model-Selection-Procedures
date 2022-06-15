

###Version History########

## 2022-06-15	Ver 1.0 has released (creator: Shunichiro Orihara, Yokohama City Univ.)

##########################


###READ ME################

## This is a function for 2SRI & LIML estimator with model selection (see Orihara et al., 2022).
## Please prepare a matrix of expectation variables for treatment & outcome models.

## TT: treatment variable, XX1: matrix of expectation variables for the treatment model
## YY: outcome variable,   XX2: matrix of expectation variables for the outcome model

##########################


###PROGRAMS###############

### 2SRI estimator with model selection ###
TSRI_func <- function(XX1,XX2,family="probit",model_select=1){   ### family: assumed outcome model, model_select==1: AIC; else: BIC 
TS1 <- lm(TT~-1+XX1); TT_rdl <- TS1$resi   ### 1st stage model
TS2 <- glm(YY~-1+TT_rdl+XX2,family=binomial(family))

if(model_select==1){
  LIC1 <- AIC(TS1); LIC2 <- AIC(TS2)
}
else{
  LIC1 <- BIC(TS1); LIC2 <- BIC(TS2)
}

alpha <- as.numeric(TS1$coef); beta <- as.numeric(TS2$coef)

kekka <- list(c(LIC1,LIC2),alpha,beta)
names(kekka) <- c("Information criterions","Coef. of the treartment model","Coef. of the outcome model")

return(kekka)
}
### End of 2SRI ###

### LIML estimator with model selection ###
LIML_func <- function(XX1,XX2,model_select=1){   ### model_select==1: AIC; else: BIC

ini=rep(0,ncol(XX1)+ncol(XX2))

optim_func <- function(theta){
  alpha <- theta[1:(ncol(XX1))]; beta <- theta[(ncol(XX1)+1):length(ini)]
  sigma <- theta[length(ini)+1]; rho <- theta[length(ini)+2]
  
  rho2 <- 1-rho^2; rho2_ <- replace(rho2,0.01>rho2,0.01)   ### Adjust of correlation of unmeasured covariates
  rho1 <- replace(rho2_,0.99<rho2_,0.99)
  
  nprob1_ <- pnorm((XX2%*%beta+rho*(TT-XX1%*%alpha))/sqrt(rho1))
  nprob0_ <- 1-pnorm((XX2%*%beta+rho*(TT-XX1%*%alpha))/sqrt(rho1))
  
  nprob1 <- replace(nprob1_,nprob1_<=0,2^-1074)
  nprob0 <- replace(nprob0_,nprob0_<=0,2^-1074)
  
  llik <- c(YY*log(nprob1)+(1-YY)*log(nprob0)-(TT-XX1%*%alpha)^2/(2*sigma^2)-log(sqrt(2*pi*sigma^2)))
  
  return(-sum(llik))
}
YY_coef <- optim(par=c(ini,1,0.2),fn=optim_func,method="L-BFGS-B",lower=c(rep(-Inf,length(ini)),0.1,-0.99),upper=c(rep(Inf,length(ini)),5,0.99))$par
YY_para <- length(YY_coef)

alpha <- YY_coef[1:(ncol(XX1))]; beta <- YY_coef[(ncol(XX1)+1):length(ini)]
sigma <- YY_coef[length(ini)+1]; rho <- YY_coef[length(ini)+2]

nprob1 <- pnorm((XX2%*%beta+rho*(TT-XX1%*%alpha))/sqrt(1-rho^2))
nprob0 <- 1-pnorm((XX2%*%beta+rho*(TT-XX1%*%alpha))/sqrt(1-rho^2))

llik <- c(YY*log(nprob1)+(1-YY)*log(nprob0)-(TT-XX1%*%alpha)^2/(2*sigma^2)-log(sqrt(2*pi*sigma^2)))

if(model_select==1) LIC <- -2*sum(llik)+2*YY_para
else LIC <- -2*sum(llik)+log(nn)*YY_para

kekka <- list(LIC,alpha,beta,sigma,rho)
names(kekka) <- c("Information criterion","Coef. of the treartment model","Coef. of the outcome model",
				  "Variance of the unmeasured covariate","Correlation of the unmeasured covariates")
				  
return(kekka)
}

### End of LIML ###

##########################


###EXAMPLES###############

#rm(list=ls(all.names=T))
set.seed(1234)

nn <- 1000

#Covariates
XX1 <- rnorm(nn)
XX2 <- rbinom(nn,1,0.5)
XX3 <- rnorm(nn)

#Unmeasured covariates
rho_1 <- 0.3; sigma_1 <- c(1,rho_1,rho_1,1); Sigma_1 <- matrix(sigma_1,ncol=2)
UU_1 <- mvrnorm(nn,c(0,0),Sigma_1)

#Valid IV
ZZ <-  rbinom(nn,1,0.5)

#Treatment
TT_11 <- 1+0.2*ZZ+XX2+XX3+UU_1[,1]

#Outcome
YY_lm_1111 <- 0.5+0.2*TT_11+0.5*(XX1+XX2)+UU_1[,2]
YY_1111 <- rep(0,nn)
for(ii in 1:nn){
  if(YY_lm_1111[ii]>=0) YY_1111[ii] <- 1
  else YY_1111[ii] <- 0
}

TT <- TT_11; XX_1 <- cbind(1,ZZ,XX2,XX3)
YY <- YY_1111; XX_2 <- cbind(1,TT,XX1,XX2)

TSRI_func(XX_1,XX_2)
LIML_func(XX_1,XX_2)
##########################


