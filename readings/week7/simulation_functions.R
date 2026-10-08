######################################
# (a) the data generating function   #
######################################
generate_data<-function(TYPE, subjects, Timesteps, Replicates, HETEROGENEITY, MEASUREMENT_ERROR, parA, parB, parC,  parD, parF, parG, errorEpsilon, errorLambda, errorKappa, varianceNu, varianceMu, covarMuNu, POP) {
  # see Box 1 in main text for equations
  if (HETEROGENEITY==FALSE)  { varianceNu<-varianceMu<-covarMuNu<-0 }
  ## Generate the data
  varnames<-c("X","Y", "Z", "SubjectID",  "timestep", "Xlagged", "Ylagged", "Zlagged",  "YXlagged", "ZXlagged","Xoriginal", "Yoriginal", "YX", "ZX", "deviationX", "deviationYlagged", "deviationY", "deviationYXlagged", "deviationZXlagged", "meanX", "meanY", "meanZ", "sdX", "sdY", "sdZ","sdXY", "SubjMeanX", "SubjMeanY", "SubjMeanZ")
  datalong<- array(NA,dim = c(subjects*Timesteps,length(varnames), Replicates), dimnames=list(NULL,varnames))
  
  for(r in 1: Replicates) {
    deviation_subject <- rmvnorm(n=subjects, mean=c(0,0), sigma=matrix(c(varianceNu,covarMuNu,covarMuNu,varianceMu), ncol=2), method="eigen")
    errorY<-rnorm(subjects*Timesteps,0,errorEpsilon) 
    errorX<-rnorm(subjects*Timesteps,0,errorKappa) 
    errorZ<-rnorm(subjects*Timesteps,0,errorLambda) 
    counter <- 0
    for(i in 1:subjects) {
      for (j in 1:Timesteps) {
        counter<-counter+1
        datalong[counter,"SubjectID",r]<-i
        datalong[counter,"timestep",r]<-j
        if(TYPE=="TRADEOFF") { 
          ifelse(j==1, datalong[counter,"X",r]<-parC+deviation_subject[i,1]+errorX[counter], datalong[counter,"X",r]<-parC+parD*datalong[(counter-1),"Y",r]+deviation_subject[i,1]+errorX[counter])
          datalong[counter,"Y",r]<-parA+parB*datalong[counter,"X",r]+deviation_subject[i,2]+errorY[counter]
        }       
        if(TYPE=="GROUP") { 
          ifelse(j==1, datalong[counter,"X",r]<-rpois(1,POP),datalong[counter,"X",r]<-parC+parD*datalong[(counter-1),"Y",r]+parF*datalong[(counter-1),"Z",r]*datalong[(counter-1),"X",r]+errorX[counter])
          datalong[counter,"Y",r]<-parA+parB*datalong[counter,"X",r]+deviation_subject[i,2]+errorY[counter]
          datalong[counter,"Z",r]<-parG+deviation_subject[i,1]+errorZ[counter]
        }
        if(TYPE=="DENSDEP")  { 
          ifelse(j==1, datalong[counter,"X",r]<-rpois(1,POP), datalong[counter,"X",r]<-parC+parD*datalong[(counter-1),"Y",r]*datalong[(counter-1),"X",r]+parF*datalong[(counter-1),"Z",r]*datalong[(counter-1),"X",r]+errorX[counter])
          datalong[counter,"Y",r]<-parA+ parB*datalong[counter,"X",r]+deviation_subject[i,2]+errorY[counter]
          datalong[counter,"Z",r]<-parG+ deviation_subject[i,1]+errorZ[counter]
        }
      }
    }
    
    # rescale all variable X, Y and Z to have zero mean and sd=1
    datalong[,"meanX",r]<-mean(datalong[,"X",r])
    datalong[,"meanY",r]<-mean(datalong[,"Y",r])
    datalong[,"meanZ",r]<-mean(datalong[,"Z",r])
    datalong[,"sdX",r]<-sd(datalong[,"X",r])
    datalong[,"sdY",r]<-sd(datalong[,"Y",r])
    datalong[,"sdZ",r]<-sd(datalong[,"Z",r])
    datalong[,"X",r]<-(datalong[,"X",r]-mean(datalong[,"X",r]))/sd(datalong[,"X",r]) # z-scores
    datalong[,"Y",r]<-(datalong[,"Y",r]-mean(datalong[,"Y",r]))/sd(datalong[,"Y",r]) # z-scores
    datalong[,"Z",r]<-(datalong[,"Z",r]-mean(datalong[,"Z",r]))/sd(datalong[,"Z",r]) # z-scores
    
    # add measurement error if needed
    if(MEASUREMENT_ERROR=="X") { 
      errorXobs<-sqrt((1-Reliability)*var(datalong[,"X",r]))     #  measurement error X 
      datalong[,"Xoriginal",r]<-datalong[,"X",r]
      datalong[,"X",r]<-datalong[,"X",r]+rnorm(length(datalong[,"X",r]),0,errorXobs)
    }
    if (MEASUREMENT_ERROR=="Y") {  
      errorYobs<-sqrt((1-Reliability)*var(datalong[,"Y",r]))     #  measurement error Y 
      datalong[,"Yoriginal",r]<-datalong[,"Y",r]
      datalong[,"Y",r]<-datalong[,"Y",r]+rnorm(length(datalong[,"Y",r]),0,errorYobs)
    }
    # calculate derived variables (lagged and centered terms)
    selection<-which(datalong[,"timestep",r]>1) 
    datalong[selection,"Ylagged",r]<-datalong[(selection-1),"Y",r]
    datalong[selection,"Xlagged",r]<-datalong[(selection-1),"X",r]
    if(TYPE!="TRADEOFF") { datalong[selection,"Zlagged",r]<-datalong[(selection-1),"Z",r]}
    datalong[,"YX",r]<-datalong[,"X",r]*datalong[,"Y",r]
    datalong[,"sdXY",r]<-sd(datalong[,"YX",r])
    datalong[,"ZX",r]<-datalong[,"X",r]*datalong[,"Z",r]
    datalong[selection,"YXlagged",r]<-datalong[(selection-1),"YX",r]
    datalong[selection,"ZXlagged",r]<-datalong[(selection-1),"ZX",r]
    datalong[,"deviationX",r]<-wgdev(datalong[,"X",r], datalong[,"SubjectID",r]) # calculates the within subject deviation in X by subtracting the subject mean X of all X values.
    datalong[,"deviationYlagged",r]<-wgdev(datalong[,"Ylagged",r], datalong[,"SubjectID",r])+mean(datalong[,"Ylagged",r], na.rm=T)
    datalong[,"deviationY",r]<-wgdev(datalong[,"Y",r], datalong[,"SubjectID",r])+mean(datalong[,"Y",r], na.rm=T)
    datalong[,"deviationY",r]<-wgdev(datalong[,"Y",r], datalong[,"SubjectID",r])+mean(datalong[,"Y",r], na.rm=T)
    datalong[,"deviationYXlagged",r]<-wgdev(datalong[,"YXlagged",r], datalong[,"SubjectID",r])+mean(datalong[,"YXlagged",r], na.rm=T)
    datalong[,"deviationZXlagged",r]<-wgdev(datalong[,"ZXlagged",r], datalong[,"SubjectID",r])+mean(datalong[,"ZXlagged",r], na.rm=T)
    datalong[,"SubjMeanX",r]<-datalong[,"X",r]-datalong[,"deviationX",r]
    datalong[,"SubjMeanY",r]<-wgmean(datalong[,"Y",r], datalong[,"SubjectID",r])
    datalong[,"SubjMeanZ",r]<-wgmean(datalong[,"Z",r], datalong[,"SubjectID",r])
  }
  return(datalong)
}

######################################
# (b) the lavaan_data  function      #
######################################
lavaan_data<-function(data, Timesteps, subjects, TYPE) {
  lavaanData<-as.data.frame(matrix(data = NA, nrow = subjects, ncol = Timesteps*3, byrow = FALSE, dimnames = NULL))
  ifelse(TYPE=="TRADEOFF", colnames(lavaanData) <-c(sprintf("y%d", 1:Timesteps),sprintf("x%d", 1:Timesteps)), colnames(lavaanData) <-c(sprintf("y%d", 1:Timesteps),sprintf("x%d", 1:Timesteps),sprintf("z%d", 1:Timesteps)))
  for (i in 1:Timesteps) {
    aaa<-subset(data,data$timestep==i)
    lavaanData[,i]<-aaa$Y
    lavaanData[,(i+Timesteps)]<-aaa$X
    if(TYPE!="TRADEOFF") {  
      lavaanData[,(i+Timesteps+Timesteps)]<-aaa$Z   
       }
  }
  return(lavaanData)
}

######################################
# (c) the stan_data  function        #
######################################
stan_data<-function(data, subjects) {
  stanData <- list(subj = data$SubjectID,   ## subject identifier 
                   Y = data$Y,              ## response variable equation 1
                   X = data$X,              ## response variable equation 2 
                   Z = data$Z,              ## response variable equation 3
                   N = nrow(data),          ## number of observations
                   L = max(data$timestep),
                   S = subjects) ## number of subjects
  return(stanData)
}

######################################
# (d) the lavaan_model  function     #
######################################
lavaan_model<-function(TYPE, Timesteps, ERROR) {
  if(TYPE=="TRADEOFF") {
    if(ERROR=="NONE") {
      Tline1<-paste("Mu =~1*y1 ",paste(sprintf("+1*y%d", 2:Timesteps), collapse = ''))
      Tline2<-paste("Nu =~1*x1 ",paste(sprintf("+1*x%d", 2:Timesteps), collapse = ''))
      Tline3<-"Mu~~varMu*Mu" #intercept factor variance of y 
      Tline4<-"Nu~~varNu*Nu" #intercept factor variance of x
      Tline5<-"Mu~~covMuNu*Nu" # covariance intercepts factors x and y
      Tline6<-paste(sprintf("y%d~inty*1+b*x%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline7<-"x1~intx1*1"  #intercept only model for t=1
      Tline8<-paste(sprintf("x%d~intx*1+d*y%d", 2:Timesteps, (2:Timesteps)-1), sep="\n") 
      Tline9<-paste(sprintf("y%d~~vary*y%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline10<-"x1~~varx1*x1"  #let lavaan estimate the variance of the x var, but for varx at t=1 is different than for t>1
      Tline11<-paste(sprintf("x%d~~varx*x%d", 2:Timesteps, 2:Timesteps), sep="\n")  # let lavaan estimate the variance of the x variables
      LavModel<-paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(Tline1,Tline2), sep=""),Tline3), sep=""),Tline4), sep=""),Tline5), sep=""),Tline6), sep=""),Tline7), sep=""),Tline8), sep=""),Tline9), sep=""),Tline10), sep=""),Tline11), sep="")
    }
    if(ERROR=="X") {
      Tline1<-paste(sprintf("xL%d =~1*x%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline2<-paste("Mu =~1*y1",paste(sprintf("+1*y%d", 2:Timesteps), collapse = ''))
      Tline3<-paste("Nu =~1*xL1 ",paste(sprintf("+1*xL%d", 2:Timesteps), collapse = ''))
      Tline4<-"Mu~~varMu*Mu" #intercept factor variance of y 
      Tline5<-"Nu~~varNu*Nu" #intercept factor variance of x
      Tline6<-"Mu~~covMuNu*Nu" # covariance intercepts factors x and y
      Tline7<-paste(sprintf("y%d~inty*1+b*xL%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline8<-"xL1~intx1*1"  #intercept only model for t=1
      Tline9<-paste(sprintf("xL%d~intx*1+d*y%d", 2:Timesteps, (2:Timesteps)-1), sep="\n") 
      Tline10<-"xL1~~varx1*xL1"  #let lavaan estimate the variance of the x var, but for varx at t=1 is different than for t>1
      Tline11<-paste(sprintf("xL%d~~varx*xL%d", 2:Timesteps, 2:Timesteps), sep="\n")  # let lavaan estimate the variance of the x variables
      Tline12<-paste(sprintf("y%d~~vary*y%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline13<-paste(sprintf("x%d~~%f*x%d", 1:Timesteps, (1-Reliability), 1:Timesteps), sep="\n")
      LavModel<-paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(Tline1,Tline2), sep=""),Tline3), sep=""),Tline4), sep=""),Tline5), sep=""),Tline6), sep=""),Tline7), sep=""),Tline8), sep=""),Tline9), sep=""),Tline10), sep=""),Tline11), sep=""),Tline12), sep=""),Tline13), sep="")
    }
    if(ERROR=="Y") {
      Tline1<-paste(sprintf("yL%d =~1*y%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline2<-paste("Mu =~1*yL1",paste(sprintf("+1*yL%d", 2:Timesteps), collapse = ''))
      Tline3<-paste("Nu =~1*x1 ",paste(sprintf("+1*x%d", 2:Timesteps), collapse = ''))
      Tline4<-"Mu~~varMu*Mu" #intercept factor variance of y 
      Tline5<-"Nu~~varNu*Nu" #intercept factor variance of x
      Tline6<-"Mu~~covMuNu*Nu" # covariance intercepts factors x and y
      Tline7<-paste(sprintf("yL%d~inty*1+b*x%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline8<-"x1~intx1*1"  #intercept only model for t=1
      Tline9<-paste(sprintf("x%d~intx*1+d*yL%d", 2:Timesteps, (2:Timesteps)-1), sep="\n") 
      Tline10<-"x1~~varx1*x1"  #let lavaan estimate the variance of the x var, but for varx at t=1 is different than for t>1
      Tline11<-paste(sprintf("x%d~~varx*x%d", 2:Timesteps, 2:Timesteps), sep="\n")  # let lavaan estimate the variance of the x variables
      Tline12<-paste(sprintf("yL%d~~vary*yL%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Tline13<-paste(sprintf("y%d~~%f*y%d", 1:Timesteps, (1-Reliability), 1:Timesteps), sep="\n")
      LavModel<-paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(Tline1,Tline2), sep=""),Tline3), sep=""),Tline4), sep=""),Tline5), sep=""),Tline6), sep=""),Tline7), sep=""),Tline8), sep=""),Tline9), sep=""),Tline10), sep=""),Tline11), sep=""),Tline12), sep=""),Tline13), sep="")
    }
  }
  if(TYPE=="GROUP") {
      Gline1<-paste("Mu =~1*y1 ",paste(sprintf("+1*y%d", 2:Timesteps), collapse = ''))
      Gline2<-paste("Nu =~1*z1 ",paste(sprintf("+1*z%d", 2:Timesteps), collapse = ''))
      Gline3<-paste("Ksi =~1*x1 ",paste(sprintf("+1*x%d", 2:Timesteps), collapse = ''))
      Gline4<-"Mu~~varMu*Mu" #intercept factor variance of y 
      Gline5<-"Nu~~varNu*Nu" #intercept factor variance of z
      Gline6<-"Ksi~~varKsi*Ksi" #intercept factor variance of x
      Gline7<-"Mu~~covMuNu*Nu" # covariacne intercepts factors z and y
      Gline8<-"Mu~~0*Ksi" # we constrain the covariance between y and x to be zero
      Gline9<-"Nu~~0*Ksi"# we constrain the covariance between z and x to be zero
      Gline10<-paste(sprintf("y%d~inty*1+b*x%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Gline11<-"x1~intx1*1"
      Gline12<-paste(sprintf("x%d~intx*1+d*y%d +f*z%d:x%d", 2:Timesteps, (2:Timesteps)-1,(2:Timesteps)-1 ,(2:Timesteps)-1), sep="\n") 
      Gline13<-paste(sprintf("z%d~intz*1", 1:Timesteps, 1:Timesteps), sep="\n")
      Gline14<-paste(sprintf("y%d~~vary*y%d", 1:Timesteps, 1:Timesteps), sep="\n")
      Gline15<-"x1~~varx1*x1"
      Gline16<-paste(sprintf("x%d~~varx*x%d", 2:Timesteps, 2:Timesteps), sep="\n")
      Gline17<-paste(sprintf("z%d~~varz*z%d", 1:Timesteps, 1:Timesteps), sep="\n")
      LavModel<-paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(Gline1,Gline2), sep=""),Gline3), sep=""),Gline4), sep=""),Gline5), sep=""),Gline6), sep=""),Gline7), sep=""),Gline8), sep=""),Gline9), sep=""),Gline10), sep=""),Gline11), sep=""),Gline12), sep=""),Gline13), sep=""),Gline14), sep=""),Gline15), sep=""),Gline16), sep=""),Gline17), sep="")
  }
  if(TYPE=="DENSDEP") {
    Dline1<-paste("Mu =~1*y1 ",paste(sprintf("+1*y%d", 2:Timesteps), collapse = ''))
    Dline2<-paste("Nu =~1*z1 ",paste(sprintf("+1*z%d", 2:Timesteps), collapse = ''))
    Dline3<-paste("Ksi =~1*x1 ",paste(sprintf("+1*x%d", 2:Timesteps), collapse = ''))
    Dline4<-"Mu~~varMu*Mu" #intercept factor variance of y 
    Dline5<-"Nu~~varNu*Nu" #intercept factor variance of z
    Dline6<-"Ksi~~varKsi*Ksi" #intercept factor variance of x
    Dline7<-"Mu~~covMuNu*Nu" # covariacne intercepts factors z and y
    Dline8<-"Mu~~0*Ksi"
    Dline9<-"Nu~~0*Ksi"
    Dline10<-paste(sprintf("y%d~inty*1+b*x%d", 1:Timesteps, 1:Timesteps), sep="\n")
    Dline11<-"x1~intx1*1"
    Dline12<-paste(sprintf("x%d~intx*1+d*y%d:x%d +f*z%d:x%d", 2:Timesteps, (2:Timesteps)-1,(2:Timesteps)-1,(2:Timesteps)-1 ,(2:Timesteps)-1), sep="\n") 
    Dline13<-paste(sprintf("z%d~intz*1", 1:Timesteps, 1:Timesteps), sep="\n")
    Dline14<-paste(sprintf("y%d~~vary*y%d", 1:Timesteps, 1:Timesteps), sep="\n")
    Dline15<-"x1~~varx1*x1"
    Dline16<-paste(sprintf("x%d~~varx*x%d", 2:Timesteps, 2:Timesteps), sep="\n")
    Dline17<-paste(sprintf("z%d~~varz*z%d", 1:Timesteps, 1:Timesteps), sep="\n")
    LavModel<-paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(paste(c(Dline1,Dline2), sep=""),Dline3), sep=""),Dline4), sep=""),Dline5), sep=""),Dline6), sep=""),Dline7), sep=""),Dline8), sep=""),Dline9), sep=""),Dline10), sep=""),Dline11), sep=""),Dline12), sep=""),Dline13), sep=""),Dline14), sep=""),Dline15), sep=""),Dline16), sep=""),Dline17), sep="")
  }
  return(LavModel)
}

######################################
# (e) the stan_model  function       #
######################################
stan_model<-function(TYPE, Timesteps, subjects) {
  if(TYPE=="TRADEOFF") {
    # for multiple subjects
    stanModel<-"data {
  int<lower=1> N; //number of data points  
  int<lower=1> S; //number of subjects  
  int<lower=1> L; //number of timesteps  
  vector [N] Y;   
  vector [N] X;
}
parameters {  
  vector[5] beta; //fixed intercept and slope first equation (Y~X)  
  vector<lower=0>[3] sigma_e; // error sd of response vairables  
  vector<lower=0>[2] sigma_u; //among subject sd  
  cholesky_factor_corr[2] L_u;  
  matrix[2,S] z_u;
}
transformed parameters {  
  matrix[2,S] u;  
  u = diag_pre_multiply(sigma_u, L_u) * z_u; //subject random effects  
}
model {  
  //priors  
  sigma_e ~ normal(0,1);  
  beta ~ normal(0,2);  
  L_u ~ lkj_corr_cholesky(2.0);  
  to_vector(z_u) ~ normal(0,1);  
  //likelihood  
  for (s in 1:S) {    
    Y[(1+(s-1)*L):(L+(s-1)*L)]~ normal(beta[1]+ beta[2] * X[(1+(s-1)*L):(L+(s-1)*L)] + u[1,s], sigma_e[1]);  
    X[(2+(s-1)*L):(L+(s-1)*L)] ~ normal(beta[3] + beta[4]* Y[(1+(s-1)*L):((L-1)+(s-1)*L)]+ u[2,s], sigma_e[2]);      // equation 2 for t>1 
    X[(1+(s-1)*L)] ~ normal(beta[5]+ u[2,s], sigma_e[3]);   // equation 2 for t=1
  }
}
"

# for single subject time series
stanModel1<-"data {
  int<lower=1> N; //number of data points
  int<lower=1> L; //number of timesteps
  vector [N] Y; 
  vector [N] X;
}
parameters {
  vector[5] beta; //fixed intercept and slope first equation (Y~X)
  vector<lower=0>[3] sigma_e; // error sd of response vairables  
}
model {
  //priors
  sigma_e ~ normal(0,1);
  beta ~ normal(0,2);
  //likelihood
  Y[1:L] ~ normal(beta[1]+ beta[2] * X[1:L], sigma_e[1]);     
  X[2:L] ~ normal(beta[3] + beta[4]* Y[1:(L-1)], sigma_e[2]);
  X[1] ~ normal(beta[5], sigma_e[3]);
}
"
  }

if(TYPE=="GROUP") {
  # for multiple subjects
  stanModel<-"data {
  int<lower=1> N; //number of data points
  int<lower=1> S; //number of subjects
  int<lower=1> L; //number of timesteps
  vector [N] Y; 
  vector [N] X;
  vector [N] Z; 
}
parameters {
  vector[7] beta; //fixed intercept and slope first equation (Y~X)
  vector<lower=0>[4] sigma_e; // error sd of response vairables  
  real<lower=0> sigma_v; //sd of random intercept effect for X
  vector[S] eta;
  vector<lower=0>[2] sigma_u; //among subject sd
  cholesky_factor_corr[2] L_u;
  matrix[2,S] z_u;
}
transformed parameters {
  matrix[2,S] u;
  vector[S] v;
  u = diag_pre_multiply(sigma_u, L_u) * z_u; //subj random effects
  v = sigma_v * eta;
}
model {
  //priors
  sigma_e ~ normal(0,1);
  beta ~ normal(0,2);
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  eta ~ normal(0,1);
  //likelihood
  for (s in 1:S) {
    Y[(1+(s-1)*L):(L+(s-1)*L)] ~ normal(beta[1] + beta[2] * X[(1+(s-1)*L):(L+(s-1)*L)] + u[1,s], sigma_e[1]);     // equation 1
    X[(2+(s-1)*L):(L+(s-1)*L)] ~ normal(beta[3] + beta[4] * Y[(1+(s-1)*L):((L-1)+(s-1)*L)]+ beta[6]* X[(1+(s-1)*L):((L-1)+(s-1)*L)].*Z[(1+(s-1)*L):((L-1)+(s-1)*L)] + v[s], sigma_e[2]); // equation 2f or t>1
    X[(1+(s-1)*L)] ~ normal(beta[5] + v[s], sigma_e[3]);  // equation for t=1
    Z[(1+(s-1)*L):(L+(s-1)*L)] ~ normal(beta[7] + u[2,s], sigma_e[4]);      
  }
}
"
# for single subjects
stanModel1<-"data {
  int<lower=1> N; //number of data points
  int<lower=1> L; //number of timesteps
  vector [N] Y; 
  vector [N] X;
  vector [N] Z; 
}
parameters {
  vector[6] beta; //fixed intercept and slope first equation (Y~X)
  vector<lower=0>[3] sigma_e; // error sd of response vairables  
}
model {
  //priors
  sigma_e ~ normal(0,2);
  beta ~ normal(0,1);
  //likelihood
  Y[1:L] ~ normal(beta[1] + beta[2]*X[1:L], sigma_e[1]);     
  X[2:L] ~ normal(beta[3] + beta[4]*Y[1:(L-1)] + beta[6]*X[1:(L-1)].*Z[1:(L-1)], sigma_e[2]);
  Z[1:L] ~ normal(beta[5], sigma_e[3]);     
}
"
}  

if(TYPE=="DENSDEP") {
  # for multiple subjects
  stanModel<-"data {
  int<lower=1> N; //number of data points
  int<lower=1> S; //number of subjects
  int<lower=1> L; //number of timesteps
  vector [N] Y; 
  vector [N] X;
  vector [N] Z; 
}
parameters {
  vector[7] beta; //fixed intercept and slope first equation (Y~X)
  vector<lower=0>[4] sigma_e; // error sd of response vairables  
  real<lower=0> sigma_v; //sd of random intercept effect for X
  vector[S] eta;
  vector<lower=0>[2] sigma_u; //among subject sd
  cholesky_factor_corr[2] L_u;
  matrix[2,S] z_u;
}
transformed parameters {
  matrix[2,S] u;
  vector[S] v;
  u = diag_pre_multiply(sigma_u, L_u) * z_u; //subj random effects
  v = sigma_v * eta;
}
model {
  //priors
  sigma_e ~ normal(0,1);
  beta ~ normal(0,2);
  L_u ~ lkj_corr_cholesky(2.0);
  to_vector(z_u) ~ normal(0,1);
  eta ~ normal(0,1);
  //likelihood
  for (s in 1:S) {
    Y[(1+(s-1)*L):(L+(s-1)*L)] ~ normal(beta[1] + beta[2] * X[(1+(s-1)*L):(L+(s-1)*L)] + u[1,s], sigma_e[1]);     // equation 1
    X[(2+(s-1)*L):(L+(s-1)*L)] ~ normal(beta[3] + beta[4] * Y[(1+(s-1)*L):((L-1)+(s-1)*L)].*X[(1+(s-1)*L):((L-1)+(s-1)*L)]+ beta[6]* X[(1+(s-1)*L):((L-1)+(s-1)*L)].*Z[(1+(s-1)*L):((L-1)+(s-1)*L)] + v[s], sigma_e[2]); // equation 2f or t>1
    X[(1+(s-1)*L)] ~ normal(beta[5] + v[s], sigma_e[3]);  // equation for t=1
    Z[(1+(s-1)*L):(L+(s-1)*L)] ~ normal(beta[7] + u[2,s], sigma_e[4]);      
  }
}
"
# for single subjects
stanModel1<-"data {
  int<lower=1> N; //number of data points
  int<lower=1> L; //number of timesteps
  vector [N] Y; 
  vector [N] X;
  vector [N] Z; 
}
parameters {
  vector[7] beta; //fixed intercept and slope first equation (Y~X)
  vector<lower=0>[4] sigma_e; // error sd of response vairables  
}
model {
  //priors
  sigma_e ~ normal(0,1);
  beta ~ normal(0,2);
  //likelihood
  Y[1:L] ~ normal(beta[1] + beta[2]*X[1:L], sigma_e[1]);     
  X[2:L] ~ normal(beta[3] + beta[4]*Y[1:(L-1)].*X[1:(L-1)] + beta[6]*X[1:(L-1)].*Z[1:(L-1)], sigma_e[2]);
  X[1] ~ normal(beta[5] , sigma_e[3]);
  Z[1:L] ~ normal(beta[7], sigma_e[4]);     
}
"
}
ifelse(subjects==1, StaModel<-stanModel1, StaModel<-stanModel)
return(StaModel)
}

######################################
# (f) the summarize_output function  #
######################################
summarize_output<-function(estimatesB, significanceB, methodnames, TYPES, TYPE,  subjects, Timesteps, Replicates, HETEROGENEITY, MEASUREMENT_ERROR, parB, parD, covarMuNu, POP) {
  comparisonnames <- c("method","median", "%bias", "mean", "mean_bias", "2.5%", "25%", "75%", "97.5%", "SD", "SignCorr", "NonSign","SignOppo", "TYPE", "Subjects", "Timesteps", "Replicates", "HETEROGENEITY", "MEASERROR", "parB", "parD", "Covar")
  comparisonrows<-methodnames[1:(which(methodnames=="sdX")-1)]
  comparison<-as.data.frame(matrix(data = NA, nrow = length(comparisonrows), ncol = length(comparisonnames), byrow = FALSE, dimnames = list(NULL, comparisonnames)))
  comparison[,"method"]<-comparisonrows
  for (i in 1:length(comparisonrows)) { 
    comparison[i,"median"]<-median(estimatesB[,i]*(estimatesB[,"sdY"]/estimatesB[,"sdX"])) 
    comparison[i,"mean"]<-mean(estimatesB[,i]*(estimatesB[,"sdY"]/estimatesB[,"sdX"]))
    comparison[i,"2.5%"]<-quantile(estimatesB[,i]*(estimatesB[,"sdY"]/estimatesB[,"sdX"]),0.025, na.rm=T)
    comparison[i,"25%"]<-quantile(estimatesB[,i]*(estimatesB[,"sdY"]/estimatesB[,"sdX"]),0.25, na.rm=T)
    comparison[i,"75%"]<-quantile(estimatesB[,i]*(estimatesB[,"sdY"]/estimatesB[,"sdX"]),0.75, na.rm=T)
    comparison[i,"97.5%"]<-quantile(estimatesB[,i]*(estimatesB[,"sdY"]/estimatesB[,"sdX"]),0.975, na.rm=T)
    comparison[i,"SD"]<-sd(estimatesB[,i]*(estimatesB[,"sdY"]/estimatesB[,"sdX"]))
    if(is.na(significanceB[1,i])==FALSE) {
      comparison[i,"SignCorr"]<-length(which(significanceB[,i]==3))/length(significanceB[,i])*100
      comparison[i,"NonSign"]<-length(which(significanceB[,i]==1))/length(significanceB[,i])*100
    } else {
      comparison[i,"SignCorr"]<-comparison[i,"NonSign"]<-NA
    }
  }
  comparison[,"%bias"]<-round(100*(comparison[,"median"]-ParB[which(TYPES==TYPE)])/abs(ParB[which(TYPES==TYPE)]),1)
  comparison[,"mean_bias"]<-round(100*(comparison[,"mean"]-ParB[which(TYPES==TYPE)])/abs(ParB[which(TYPES==TYPE)]),1)
  comparison[,"TYPE"]<-TYPE
  comparison[,"Subjects"]<-subjects
  comparison[,"Timesteps"]<-Timesteps
  comparison[,"Replicates"]<-Replicates
  comparison[,"HETEROGENEITY"]<-HETEROGENEITY
  comparison[,"MEASERROR"]<-MEASUREMENT_ERROR
  comparison[,"parB"]<-ParB[which(TYPES==TYPE)]
  comparison[,"parD"]<-ParD[which(TYPES==TYPE)]
  comparison[,"Covar"]<-CovarMuNu[which(TYPES==TYPE)]
  comparison[,"SignOppo"]<-round(100-comparison[,"SignCorr"]-comparison[,"NonSign"],1)
  return(comparison)
}


######################################
# (g) the summarize_output function  #
######################################
run_sims<-function(STAN, SEPARATE) {
  # Note that the STAT_WITHIN_SEPERATE method mentioned below analyzes each time series seperately and then takes the average estimate of b across all subjects, which can be seen as a two-step version of the STAT_WITHIN model.
  methodnames<-c("STAT_WITHIN","STAT_CROSS","DYN_LDVM","DYN_SEM_LAVAAN","DYN_SEM_PLUS_LAVAAN","DYN_SEM_STAN","STAT_WITHIN_SEPARATE", "sdX", "sdY", "sdXY")  # we also save sdX, sdY and SdXY fro rescaling purposes
  for (a in 1:length(Subject_trials)) {
    for (b in 1:length(Timesteps_trials)) {
    Timesteps<-Timesteps_trials[b] # length of time series
    subjects<-Subject_trials[a] # number of subjects
    Replicates<-ceiling(SampleSize/subjects)  # we need less replicates for datasets with many subsets, as they provide more information/have larger sample size. As the number of subject varies between 1-1000, this means the number of replicates varies between 50000 and 50, as sampleSize=50000. 
    # generate the data
    datalong<-generate_data(TYPE=TYPE, subjects=subjects, Timesteps=Timesteps, Replicates=Replicates, HETEROGENEITY=HETEROGENEITY, MEASUREMENT_ERROR=MEASUREMENT_ERROR, parB=ParB[which(TYPES==TYPE)], parD=ParD[which(TYPES==TYPE)], parA=ParA[which(TYPES==TYPE)], parC=ParC[which(TYPES==TYPE)],  parF=ParF[which(TYPES==TYPE)], parG=ParG[which(TYPES==TYPE)], errorEpsilon=ErrorEpsilon[which(TYPES==TYPE)], errorLambda=ErrorLambda[which(TYPES==TYPE)], errorKappa=ErrorKappa[which(TYPES==TYPE)], covarMuNu=CovarMuNu[which(TYPES==TYPE)], varianceNu=VarianceNu[which(TYPES==TYPE)], varianceMu=VarianceMu[which(TYPES==TYPE)], POP=POP)
    # construct the models for Lavaan and Stan
    lavaanModel<-lavaan_model(TYPE, Timesteps, ERROR="NONE")
    lavaanModelPlus<-lavaan_model(TYPE, Timesteps, ERROR=MEASUREMENT_ERROR)
    stanModel<-stan_model(TYPE, Timesteps, subjects)
    # we create dataframes to store the estimates, lower and upper 95% confidence intervals and statistical significance 
    estimatesB<-CIlowerB<-CIupperB<-significanceB<-as.data.frame(matrix(data = NA, nrow = Replicates, ncol = length(methodnames) , byrow = FALSE, dimnames = list(NULL,methodnames)))
    # analyze each replicate dataset seperately with the different regression methods
    for(r in 1: Replicates) {
      output<-as.data.frame(datalong[,,r]) #take one replicate dataset
      lavaanData<-lavaan_data(data=output, Timesteps=Timesteps, subjects=subjects, TYPE=TYPE)       # transform output dataset into lavaan format
      stanData<-stan_data(data=output, subject=subjects)       # transform output dataset into Stan format
      SubjectSlopes<-rep(NA, subjects)  # vector to store the estimates of b for each subject for the STAT_WITHIN_SEPERATE method
      # run the models on the simulated datasets
      if(SEPARATE==TRUE) { 
        for(i in 1:subjects) { SubjectSlopes[i]<-coef(lm(Y~X, data=subset(output,output$SubjectID==i)))[2] }
      }
      if(subjects>1) {
        DYN_SEM_LAVAANmodel<-lavaan(lavaanModel, data=lavaanData, missing = 'ML', int.ov.free = F,  int.lv.free = F,auto.fix.first = F, auto.fix.single = F,auto.cov.lv.x = F, auto.cov.y = F, auto.var = F)
        if(MEASUREMENT_ERROR!="NONE") { DYN_SEM_PLUS_LAVAANmodel<-lavaan(lavaanModelPlus, data=lavaanData, missing = 'ML', int.ov.free = F,  int.lv.free = F,auto.fix.first = F, auto.fix.single = F,auto.cov.lv.x = F, auto.cov.y = F, auto.var = F) }          
        if(STAN==TRUE) { DYN_SEM_STANmodel<-stan(model_code = stanModel,  data = stanData, iter = StanIter, chains = StanChains,  refresh = 0) }
        STAT_CROSSmodel<-lmer(Y~X+(1|SubjectID), data=output)
        STAT_WITHINmodel<-lmer(Y~deviationX+(1|SubjectID), data=output)
        DYN_LDVMmodel<-lmer(Y~deviationX+deviationYlagged+(1|SubjectID), data=output)
      } else {
        DYN_SEM_STANmodel<-stan(model_code = stanModel,  data = stanData, iter = StanIter, chains = StanChains,  refresh = 0, control = list(adapt_delta = 0.9999))
        STAT_CROSSmodel<-STAT_WITHINmodel<-lm(Y~X, data=output)
        DYN_LDVMmodel<-lm(Y~X+Ylagged, data=output)
      }
      estimatesB[r, "sdX"]<-output$sdX[1]
      estimatesB[r, "sdY"]<-output$sdY[1]
      estimatesB[r, "sdXY"]<-output$sdXY[1]
      # estimates for STAT_WITHIN_SEPARATE
      if(SEPARATE==TRUE ) {
        estimatesB[r,"STAT_WITHIN_SEPARATE"]<-median(SubjectSlopes)
        ifelse(subjects==1,CIlowerB[r,"STAT_WITHIN_SEPARATE"]<-NA, CIlowerB[r,"STAT_WITHIN_SEPARATE"]<-median(SubjectSlopes)-1.96*sd(SubjectSlopes)/length(SubjectSlopes))
        ifelse(subjects==1,CIupperB[r,"STAT_WITHIN_SEPARATE"]<-NA, CIupperB[r,"STAT_WITHIN_SEPARATE"]<-median(SubjectSlopes)+1.96*sd(SubjectSlopes)/length(SubjectSlopes))
        significanceB[r,"STAT_WITHIN_SEPARATE"]<-sign(ParB[which(TYPES==TYPE)])*(sign(ParB[which(TYPES==TYPE)])+sign(CIlowerB[r,"STAT_WITHIN_SEPARATE"])+sign(CIupperB[r,"STAT_WITHIN_SEPARATE"]))
      } else {
        estimatesB[r,"STAT_WITHIN_SEPARATE"]<-CIlowerB[r,"STAT_WITHIN_SEPARATE"]<-CIupperB[r,"STAT_WITHIN_SEPARATE"]<- significanceB[r,"STAT_WITHIN_SEPARATE"]<-NA
      }
      # estimates for DYN_SEM_LAVAAN
      ifelse(subjects==1, estimatesB[r,"DYN_SEM_LAVAAN"]<-NA,   estimatesB[r,"DYN_SEM_LAVAAN"]<-parameterestimates(DYN_SEM_LAVAANmodel)$est[which(parameterestimates(DYN_SEM_LAVAANmodel)$label=="b")[1]])
      ifelse(subjects==1, CIlowerB[r,"DYN_SEM_LAVAAN"]<-NA,   CIlowerB[r,"DYN_SEM_LAVAAN"]<-parameterestimates(DYN_SEM_LAVAANmodel)$ci.lower[which(parameterestimates(DYN_SEM_LAVAANmodel)$label=="b")[1]])
      ifelse(subjects==1, CIupperB[r,"DYN_SEM_LAVAAN"]<-NA,   CIupperB[r,"DYN_SEM_LAVAAN"]<-parameterestimates(DYN_SEM_LAVAANmodel)$ci.upper[which(parameterestimates(DYN_SEM_LAVAANmodel)$label=="b")[1]])
      significanceB[r,"DYN_SEM_LAVAAN"]<-sign(ParB[which(TYPES==TYPE)])*(sign(ParB[which(TYPES==TYPE)])+sign(CIlowerB[r,"DYN_SEM_LAVAAN"])+sign(CIupperB[r,"DYN_SEM_LAVAAN"]))
      # estimates for DYN_SEM_PLUS_LAVAAN
      if(MEASUREMENT_ERROR!="NONE") { 
        ifelse(subjects==1, estimatesB[r,"DYN_SEM_PLUS_LAVAAN"]<-NA,   estimatesB[r,"DYN_SEM_PLUS_LAVAAN"]<-parameterestimates(DYN_SEM_PLUS_LAVAANmodel)$est[which(parameterestimates(DYN_SEM_PLUS_LAVAANmodel)$label=="b")[1]])
        ifelse(subjects==1, CIlowerB[r,"DYN_SEM_PLUS_LAVAAN"]<-NA,   CIlowerB[r,"DYN_SEM_PLUS_LAVAAN"]<-parameterestimates(DYN_SEM_PLUS_LAVAANmodel)$ci.lower[which(parameterestimates(DYN_SEM_PLUS_LAVAANmodel)$label=="b")[1]])
        ifelse(subjects==1, CIupperB[r,"DYN_SEM_PLUS_LAVAAN"]<-NA,   CIupperB[r,"DYN_SEM_PLUS_LAVAAN"]<-parameterestimates(DYN_SEM_PLUS_LAVAANmodel)$ci.upper[which(parameterestimates(DYN_SEM_PLUS_LAVAANmodel)$label=="b")[1]])
        significanceB[r,"DYN_SEM_PLUS_LAVAAN"]<-sign(ParB[which(TYPES==TYPE)])*(sign(ParB[which(TYPES==TYPE)])+sign(CIlowerB[r,"DYN_SEM_LAVAAN"])+sign(CIupperB[r,"DYN_SEM_LAVAAN"]))
      }      
       # estimates for DYN_SEM_STAN
      if(STAN==TRUE | subjects==1) { 
        estimatesB[r,"DYN_SEM_STAN"]<-summary(DYN_SEM_STANmodel, pars = "beta[2]")$summary[6] 
        CIlowerB[r,"DYN_SEM_STAN"]<-summary(DYN_SEM_STANmodel, pars = "beta[2]")$summary[4]
        CIupperB[r,"DYN_SEM_STAN"]<-summary(DYN_SEM_STANmodel, pars = "beta[2]")$summary[8]
        significanceB[r,"DYN_SEM_STAN"]<-sign(ParB[which(TYPES==TYPE)])*(sign(ParB[which(TYPES==TYPE)])+sign(CIlowerB[r,"DYN_SEM_STAN"])+sign(CIupperB[r,"DYN_SEM_STAN"]))
      } else {
        estimatesB[r,"DYN_SEM_STAN"]<-CIlowerB[r,"DYN_SEM_STAN"]<-CIupperB[r,"DYN_SEM_STAN"]<- significanceB[r,"DYN_SEM_STAN"]<-NA
      }
      # estimates for STAT_WITHIN
      estimatesB[r,"STAT_WITHIN"]<-summary(STAT_WITHINmodel)$coefficients[2,1] 
      CIlowerB[r,"STAT_WITHIN"]<-estimatesB[r,"STAT_WITHIN"]-1.96*summary(STAT_WITHINmodel)$coefficients[2,2]
      CIupperB[r,"STAT_WITHIN"]<-estimatesB[r,"STAT_WITHIN"]+1.96*summary(STAT_WITHINmodel)$coefficients[2,2]
      significanceB[r,"STAT_WITHIN"]<-sign(ParB[which(TYPES==TYPE)])*(sign(ParB[which(TYPES==TYPE)])+sign(CIlowerB[r,"STAT_WITHIN"])+sign(CIupperB[r,"STAT_WITHIN"]))
      # estimates for STAT_CROSS
      estimatesB[r,"STAT_CROSS"]<-summary(STAT_CROSSmodel)$coefficients[2,1] 
      CIlowerB[r,"STAT_CROSS"]<-estimatesB[r,"STAT_CROSS"]-1.96*summary(STAT_CROSSmodel)$coefficients[2,2]
      CIupperB[r,"STAT_CROSS"]<-estimatesB[r,"STAT_CROSS"]+1.96*summary(STAT_CROSSmodel)$coefficients[2,2]
      significanceB[r,"STAT_CROSS"]<-sign(ParB[which(TYPES==TYPE)])*(sign(ParB[which(TYPES==TYPE)])+sign(CIlowerB[r,"STAT_CROSS"])+sign(CIupperB[r,"STAT_CROSS"]))
      # estimates for DYN_LDVM
      estimatesB[r,"DYN_LDVM"]<-summary(DYN_LDVMmodel)$coefficients[2,1]
      CIlowerB[r,"DYN_LDVM"]<-estimatesB[r,"DYN_LDVM"]-1.96*summary(DYN_LDVMmodel)$coefficients[2,2]
      CIupperB[r,"DYN_LDVM"]<-estimatesB[r,"DYN_LDVM"]+1.96*summary(DYN_LDVMmodel)$coefficients[2,2]
      significanceB[r,"DYN_LDVM"]<-sign(ParB[which(TYPES==TYPE)])*(sign(ParB[which(TYPES==TYPE)])+sign(CIlowerB[r,"DYN_LDVM"])+sign(CIupperB[r,"DYN_LDVM"]))
      # print an update to screen to monitor progress
      print(paste0("t",Timesteps,"_","s",subjects,"_", round(r/Replicates*100,0), "%", "_",TYPE))
    }
    #summarize the model output from different methods applied to datasets
    comparison<-summarize_output(estimatesB=estimatesB,significanceB=significanceB, methodnames=methodnames, TYPES=TYPES, TYPE=TYPE,  subjects=subjects, Timesteps=Timesteps, Replicates=Replicates, HETEROGENEITY=HETEROGENEITY, MEASUREMENT_ERROR=MEASUREMENT_ERROR, parB=parB, parD=parD, covarMuNu=covarMuNu, POP=POP)  
    ifelse( (b==1 && a==1), overview<-comparison, overview<-rbind(overview, comparison)) # append the output for each parameter combination of subects and time steps
  }
 }
return(overview)
}
