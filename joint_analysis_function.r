#Joint analysis function defined here.
#Input: Precomputed Z-scores from each data set
#Output:
#Joint:posterior probability of a gene being DE computed by joint analysis
#Single:posterior probability of a gene being DE computed by single data set analysis
#piout:estimated prior probability from ztable
#status:The DE status index of each prior probability
joint_analysis_zscore<-function(ztable,df_fit=8){
  #ztable:first column always named "GeneID"
  no_genes=dim(ztable)[1]#total number of genes measured on the microarray platform
  N_dis=dim(ztable)[2]-1#number of diseases used for joint analysis,substract 1 means that one column contains the name of genes
  print(paste(no_genes,"genes and",N_dis,"diseases will be analyzed in joint analysis"))
  zmat <- matrix(0,nrow=no_genes,ncol=N_dis)# contain the z-value pre-computed from multiple diseases
  for (i in 1:dim(zmat)[2]){
    zmat[,i]=ztable[,i+1]
  }#ztable is converted to a matrix object
  colnames(zmat)<-colnames(ztable)[-1]
  rownames(zmat)<-ztable$GeneID
  ### local fdr calculation for each disease dataset to estimate f(Z|D=1)/f(Z|D=0) density values###
  #Ultimate goal of joint anlysis: calculate posterior probability: Pr(Dn=1|z1,z2,z3....zn) for each gene
  #i.e. out of all possible ways of generating observed z1,z2,z3....zn, what is the proportion that is generated when Dn=1
  #So we need to estimate Pr(D=0|Z=z), Pr(D=1|Z=z),f(z|D=0),f(z|D=1) based on Bayes Theory
  #D=0, this gene is not disease causing gene
  #D=1, this gene is disease causing gene
  #Generally speacking, significant higher/lower expression of a gene in disease group will have a higher chance of being a disease causing gene.
  
  #1.Initialize the matrix object
  lfdr=matrix(0,nrow=no_genes,ncol=N_dis)  #Local fdr vector, posterior probability Pr(D=0|Z), for each gene in each disease
  p0=rep(0,N_dis)         #Null prior probability, Pr(D=0)   
  f0=matrix(0,no_genes,N_dis)    #Null conditional density f(Z|D=0) if the gene is NOT associated with the disease
  p1f1=matrix(0,no_genes,N_dis)    #Non-null conditional density * p1 f(Z|D=1)
  #2. calculate local false discovery rate of each gene in each dataset, i.e. posterior probability Pr(D=0|Z=z)
  for(i in 1:N_dis) {#calculate the density value f(Z|D) within each disease dataset
    if (i==1){
      Mfdr = locfdr(zmat[,i],df=df_fit)     
    }else if (i==2){
      Mfdr = locfdr(zmat[,i],df=df_fit,nulltype = 2)#central matching
    }
    lfdr[,i]=Mfdr$fdr;#This is posterior probability Pr(D=0|Z=z) for each gene
    if (i==1){#some tricks, in HD case, MLE does not work
      p0[i]=Mfdr$fp0[3,3]#Pr(D=0), prior probability estimate with MLE method
    }else{
      p0[i]=Mfdr$fp0[5,3]#Pr(D=0), parameter estimated from central matching method
    }
    f0[,i]=predict(interpSpline(Mfdr$mat[,1],Mfdr$mat[,6]),zmat[,i])$y #f(Z|D=0),need to compute it for every gene Zi score
    p1f1[,i]=predict(interpSpline( Mfdr$mat[,1],Mfdr$mat[,11]),zmat[,i])$y#f(Z|D=1)*Pr(D=1), non-null distribution density value*(1-p0) for each gene if it is associated with the disease
  }
  #Estimate f(Z|D=1) and f(Z|D=0), conditional probability density value within each dataset.
  #These density values computed from Z-scores are fixed in joint analysis
  density_0 <- matrix(0, nrow = no_genes, ncol = N_dis)
  #f(Z)=Pr(D=0)*f(Z|D=0)+Pr(D=1)*f(Z|D=1)
  for(i in 1:N_dis){
    density_0[,i] <- lfdr[,i]/p0[i]*(p0[i]*f0[,i]+p1f1[,i])#for a given z: f(Z|D=0)=Pr(D=0|Z)*f(Z)/Pr(D=0)
  }
  density_1 <- matrix(0, nrow = no_genes, ncol = N_dis)
  for(i in 1:N_dis){
    density_1[,i] <- (1-lfdr[,i])/(1-p0[i])*(p0[i]*f0[,i]+p1f1[,i])#f(Z|D=1)=(1-Pr(D=0|Z))*f(Z)/(1-Pr(D=0)) 
  }
  #3. Finally calculate posterior probablity Pr(D=1 or 0|Z1,Z2...Zn), joint analysis step###
  ### initialization ###
  #After we obtain the distribution of f(Z|D=0 or 1) in each disease dataset, we move on to joint analysis
  total_no_status=2 #means we have two conditions: DE or non-DE genes
  J <- total_no_status^N_dis#the total number of prior probability
  pivec <- rep(1/J,J)#initilization of prior probablity, initial guess for prior probability
  #likelihood variable will contain the joint density value of f(z1,z2|D1,D2)=f(z1|D1)*f(z2|D2) for each gene (each row)
  #dim=[no_genes,J]
  #each column of likelihood corresponds to each row in status
  likelihood <- matrix(0,nrow=no_genes, ncol=J)
  #generate the all combinations of prior probability
  #status dim=c(J,N_disease)
  status = data.frame(c(rep(0,J/2),rep(1,J/2)))
  for(coln in 2:N_dis){
    status <- cbind(status, rep(c(rep(0,J/(2^coln)),rep(1,J/(2^coln))),2^(coln-1)))
  }
  colnames(status) <- c(1:N_dis)
  pioutput <- array(c(0, rep(1/J,J)))
  for(i in 1:nrow(likelihood)) {#for each gene, calculate joint density value of f(Z1,Z2|D1,D2)
    output <- matrix(0, nrow=J, ncol=N_dis)#output stores all f(z|d) density, based on status variable
    for(td in 1:N_dis){#for each disease
      output[status[,td]==0,td] <- density_0[i,td]#density value of f(z|D_td=0),status[,td] picks the index where D=0
      output[status[,td]==1,td] <- density_1[i,td]#density value of f(Z|D_td=1)
    }
    #Independent assumption of joint probability density: f(z1,z2|D1,D2)=f(z1|D1)*f(z2|D2)
    #These values are fixed in joint analysis
    #compute the joint density value for each combination f(Z1|D1)*f(Z2|D2)*f(Z3|D3)...
    likelihood[i,] <- apply(output,1, prod)#The joint density value of each gene is ordered according to the row of status variable
    #Move on to next gene  
  }
  ### EM algorithm to get maximum likelihood estimate of prior probability###
  stopyn <- FALSE
  itnumber <- 1
  limit <- 100
  likelihood_increase=c()
  while(!stopyn & itnumber <=  limit){
    result <- myiteration(pivec,N=no_genes,N_dis=N_dis,likelihood = likelihood,status_matrix = status)#run one EM Update of prior probability in this line
    pivec <- result$pivec#update the new pivec
    if(result$stopyn == "y") stopyn <- TRUE
    flush.console()
    print(itnumber)
    #print(result$pivec)
    #print(result$new_total_likelihood)
    likelihood_increase=c(likelihood_increase,result$new_total_likelihood)
    pioutput <- rbind(pioutput, c(itnumber, result$pivec))#record the new estimate of prior probability
    itnumber <- itnumber + 1
  }
  mar <- result$marginalprobs#get the posterior probability for each gene i.e Pr(D1=1|Z1,Z2...Zn),Pr(D2=1|Z1,Z2...Zn)
  rownames(mar)=ztable$GeneID
  colnames(mar)=colnames(ztable)[2:3]
  single<-1-lfdr
  rownames(single)=ztable$GeneID
  colnames(single)=colnames(ztable)[2:3]
  #plot(likelihood_increase,xlab = "iteration",ylab="log likelihood")
  return(list("joint"=mar,"single"=single,"piout"=pioutput,"status_matrix"=status))
}
#one iteration of EM algorithm defined here
myiteration <- function(pivec, tol=.0001,N=no_genes,N_dis,likelihood,status_matrix){#perform one single iteration to update the pivec through EM algorithm
  #This function computes one update of EM algorithm of prior probability
  #pivec is the prior probability we want to estimate based on all observed z-scores
  #First compute posterior probability for EACH gene based on the observed Z-score in each dataset and pivec of each iteration
  #i.e. Pr(D1=1,D2=0|Z1,Z2,pivec) or Pr(D1=0,D2=1|Z1,Z2,pivec) for gene i,this is the E-step
  #Then form of M-step is derived as well
  #E-step
  #what is likelihood object:f(Z1=z1|D1)*f(Z2=z2|D2), conditional joint density pre-computed given observed Z1 and Z2
  #what is likelihoodnew: i.e. f(Z1,Z2|D1=1,D2=0)*pi(D1=1,D2=0), each row is a gene, and contains all 2^d situations
  #here pivec is the estimated value in last iteration
  likelihoodnew <- likelihood*matrix(rep(pivec,N),byrow=TRUE,nrow=N)#based on estimated prior probability of each iteration, calculate f(Z1,Z2|D1,D2)*pivec  
  #Posterior probability based on pivec: i.e. Pr(D1=1,D2=0|Z1,Z2,pivec) for each column
  conditionalprobsnew <- likelihoodnew/apply(likelihoodnew,1,sum)#compute the posterior probability for each gene (row)
  #marginalprobsnew contains Pr(Dn=1|Z1,Z2...Zn) for each gene in each disease dataset, this is the inference results
  marginalprobsnew <- matrix(0, nrow=N, ncol=N_dis)
  #calculate the marginal posterior probability: e.g. Pr(D1=1|Z1,Z2)=Pr(D1=1,D2=0|Z1,Z2)+Pr(D1=1,D2=1|Z1,Z2) for each gene
  for(j in 1:N_dis){
    marginalprobsnew[,j] <- apply(conditionalprobsnew[,status_matrix[,j]==1],1,sum)#marginalprobsnew is the output
  }
  #M-step
  #The maximizer of the log expectation function is found to be the average of the posterior probability for each gene after we observe Z. i.e. Average(Pr(Di=1|Z,old parameter))
  #This maximizer was found by Langrange multiplier
  pivecnew <- (1/N)*apply(conditionalprobsnew,2,sum)#update the prior probability, average over all genes, "2" means by average by column (gene)
  stopyn = "n"
  #output the new likelihood value of all observed sample values after prior probaility is updated
  #to validate if the updated prior probability increase the total likelihood
  new_sample_likelihood=likelihood*matrix(rep(pivec,N),byrow=TRUE,nrow=N)
  new_sample_loglikelihood=log(apply(new_sample_likelihood,1,sum))
  new_total_likelihood=sum(new_sample_loglikelihood)
  if(max(abs(pivec - pivecnew))<tol) stopyn = "y"
  return(list(new_total_likelihood=new_total_likelihood,likelihoodnew=likelihoodnew,conditionalprobs= conditionalprobsnew, marginalprobs=marginalprobsnew, pivec = pivecnew, stopyn = stopyn))
}
