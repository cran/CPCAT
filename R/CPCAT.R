#generate hypotheses represented as 0-1 contrast matrices
#' Create main hypotheses and intersection hypotheses.
#'
#' @param n The number of treatments exclusive of the control group. That is, if we have, e.g., a control group and 3 treatments choose n=3.
#'
#' @return A contrast matrix reflecting all main and intersection hypotheses induced by the closure principle (CP). Please note that the leading column of treatment 0 (i.e., the control group) is not displayed because it always contains only ones.
#' @export
#'
#' @examples hypotheses(3)
#' @examples hypotheses(n=2)
#' @importFrom utils combn
hypotheses=function(n){
  help=list()
  for(j in 1:n){
    help[[j]]=combn(1:n,j)
  }
  help2=list()
  for(j in 1:n){
    help2[[j]]=matrix(0,ncol=n,nrow=(2^(n-1)))
    colnames(help2[[j]])=paste("treatment ",1:n)
    rowindex=1
    for(i in 1:length(help)){
      for(k in 1:ncol(help[[i]])){
        if(any(help[[i]][,k]==j)){
          help2[[j]][rowindex,help[[i]][,k]]=1
          rowindex=rowindex+1
        }
      }
    }
  }
  names(help2)<-paste("H0: mu_0=mu_",1:n,sep="")
  help2
}

#' Perform a computational approach test (CAT).
#'
#' @param dat A list. The first list element contains the Poisson data of the control group, the second list element contains the Poisson data of the first treatment group etc..
#' @param contrast A matrix consisting of one row and ncol=number of treatments. The values must be either 0 or 1. 1 (0) includes (excludes) the corresponding treatment from the CAT procedure. For example consider 1 control group and 4 treatment groups. For testing H0: mu_0=mu_1=mu_4 choose contrast=matrix(c(1,0,0,1),nrow = 1).
#' @param M The number of parametric bootstrap simulations. Defaults to M=10000.
#'
#' @return A p-value for testing one intersection hypothesis.
#' @export
#'
#' @examples data(testdata)
#' @examples #Extract the group names from testdata.
#' @examples conc<-unique(testdata[,2])	#concentration levels of the test data
#' @examples #Extract the number of treatment groups.
#' @examples N=length(conc)-1
#' @examples #Devide the testdata into a list where the first element contains
#' @examples #the Poisson data of the control group and the following contain
#' @examples #the Poisson data of the treatment groups, respectively.
#' @examples dat<-list()
#' @examples #Create the final list.
#' @examples for(j in 1:length(conc))\{
#' @examples index<-which(testdata[,2]==conc[j])
#' @examples dat[[j]]<-testdata[index,]
#' @examples \}
#' @examples #Show the data list.
#' @examples dat
#' @examples #Generate the intersection hypotheses of H01: mu_0=mu_1 using the
#' @examples #hypotheses function.
#' @examples C=hypotheses(N)[[1]]
#' @examples #Show the first intersection hypothesis.
#' @examples C[1,]
#' @examples #Test the first intersection hypothesis.
#' @examples poisson.sub.test(dat=dat,contrast=C[1,],M=10000)
#' @importFrom stats rpois
poisson.sub.test<-function(dat,contrast,M=10000){
  #poisson.sub.test performs the computational approach test (CAT)
  #contrast=contrast matrix for specific H_j: mu_0=mu_j and corresponding intersection hypotheses
  #M=number of bootstrap resamples
  datsheets<-c(1,which(contrast==1)+1)
  dat2<-list()
  for(l in 1:length(datsheets)){
    dat2[[l]]<-dat[[datsheets[l]]]
  }

  #step 3
  musML<-ns<-xs<-rep(0,length(dat2))
  for(i in 1:length(musML)){
    ns[i]<-nrow(dat2[[i]])
    xs[i]<-sum(dat2[[i]][,1])
    musML[i]<-xs[i]/ns[i]
  }

  #step 4; musML[1]=control; musML[-1]=all treatments
  etaML<-sum( (sqrt(musML[-1])-sqrt(musML[1]))^2 )

  #step 5
  n<-sum(ns)
  x<-sum(xs)
  #eta0ML<-x/n
  mu0RML<-x/n

  #step 6 and 7
  artificialdata<-pseudomus<-list()
  for(j in 1:length(dat2)){
    #artificialdata[[j]]<-rpois(M,ns[j]*eta0ML)
    artificialdata[[j]]<-rpois(M,ns[j]*mu0RML)
    pseudomus[[j]]<-artificialdata[[j]]/ns[j]
  }

  #step 8
  pseudoetasML<-rep(0,M)
  for(l in 1:M){
    pseudomushelp<-vector()
    for(i in 1:length(pseudomus)){
      pseudomushelp[i]<-pseudomus[[i]][l]
    }
    pseudoetasML[l]<-sum( (sqrt(pseudomushelp[-1])-sqrt(pseudomushelp[1]))^2 )
  }

  #step 9 and 10
  pvalue<-length(which(pseudoetasML>etaML))/M
  pvalue
}

#' Perform the closure principle computational approach test (CPCAT) for one main hypothesis H0i: mu_0=mu_i.
#'
#' @param Data The data matrix.
#' @param contrastmatrix The contrasts according to the closure principle induced intersection hypotheses.
#' @param M The number of parametric bootstrap simulations. Defaults to M=10000.
#'
#' @return The set of p-values according to the intersection hypotheses and the maximum p-value.
#' @export
#'
#' @examples #Consider a data set of one control group and 3 treatment groups
#' @examples data(testdata2)
#' @examples #Test the main hypothesis H0: mu_0=mu_1 using M=10000
#' @examples #simulation runs
#' @examples C=hypotheses(3)[[1]] #Generate the set of intersection hypotheses
#' @examples #according to H0
#' @examples poisson.test(testdata2,contrastmatrix=C,M=10000)
poisson.test<-function(Data,contrastmatrix,M=10000){
  #poisson.test applies the CAT iteratively to test all intersection hypotheses corresponding to H_j: mu_0=mu_j
  #contrastmatrix=constrast matrix for specific H_j and corresponding intersection hypotheses
  #M=number of bootstrap resamples
  dat<-list()
  conc<-unique(Data[,2])	#concentration levels
  for(j in 1:length(conc)){
    index<-which(Data[,2]==conc[j])
    dat[[j]]<-Data[index,]
  }
  pvalues<-rep(0,nrow(contrastmatrix))
  for(j in 1:nrow(contrastmatrix)){
    pvalues[j]<-poisson.sub.test(dat=dat,contrast=contrastmatrix[j,],M=M)
  }
  list(contrastmatrix_pvalues=cbind(contrastmatrix,pvalues),maxpvalue=max(pvalues))
}

#' Performs the closure principle computational approach test (CPCAT).
#'
#' @param z The data set to be used. One column of z must contain the numeric Poisson data and one must contain the factor variable. The first level of the factor variable is assumed to be the control group. Factor levels (i.e., groups) should be in ascending order (e.g. increasing concentration of a test substance). If the data frame contains more than one numeric column and/or more than one factor variable the CPCAT is applied to the first numeric column and the corresponding Poisson data are grouped according the first factor variable.
#' @param M The number of parametric bootstrap simulations. Defaults to M=10000.
#'
#' @return A p-value for each main hypotheses H_0i: mu_0=mu_i "control vs. treatment i".
#' @export
#'
#' @examples data(testdata)
#' @examples CPCAT(testdata)
#' @examples CPCAT(z=testdata)
#' @examples CPCAT(z=testdata, M=1000)
CPCAT=function(z,M=10000){
  j=1
  while(is.factor((z[,j]))==FALSE & j<ncol(z)){
    j=j+1
  }
  if(is.factor(z[,j])==FALSE){
    errormessage1="Group variable missing. The group variable has to be of the data type 'factor'."
    warning(errormessage1)
  }
  k=1
  while(is.numeric((z[,k]))==FALSE & k<ncol(z)){
    k=k+1
  }
  if(is.numeric(z[,k])==FALSE){
    errormessage2="Numeric Poisson distributed data missing."
    warning(errormessage2)
  }

  if(is.factor(z[,j])==TRUE && is.numeric(z[,k])==TRUE){
    z=z[,c(k,j)]
    #if the data frame contains more than one grouping variable and more than one numeric variable
    #the first grouping variable and the first numeric variable is selected for the CPCAT procedure
    #if(is.factor(z[,j]) & is.numeric(z[,k])){
    #  if(length(unique(z[,j]))>2){#there are at least 3 groups
    #    z=data.frame(Poissondata=z[,k],Group=z[,j])
    #  }else{#there are only 2 groups --> the CPCAT is programmed for >2 groups
    #        #so add a 3rd dummy group which yields guaranteed significances and later delete the dummy results
    #    z=data.frame(Poissondata=c(z[,k],rep(100*max(z[,k]),3)),Group=factor(c(z[,j],rep(999,3))))
    #  }
    #  z[,1]=as.numeric(z[,1])
    #  z[,2]=as.factor(z[,2])
    #}

    #CPCAT applies the complete CPCAT procedure, i.e., it tests all hypotheses H_j: mu_0=mu_j (j=1,...,number of treatments) using the poisson.test function
    #the CPCAT function avoids double/triple etc. testing of intersection hypotheses which belong to different main hypotheses H_j: mu_0=mu_j
    #z is a data frame with z[,1] containing poisson data and z[,2] containing corresponding treatment names (incl. the control group)
    #the first group in z[,2] is the control group and group names should be in ascending treatment order (e.g., ascending substance concentration)
    #M=number of parametric bootstrap samples
    groupsnames=levels(factor(z[,2]))
    numbersofreplicates=rep(0,length(groupsnames))
    ntype=matrix(numbersofreplicates,nrow=1)
    for(l in 1:length(numbersofreplicates)){
      numbersofreplicates[l]=length(which(z[,2]==groupsnames[l]))
    }
    if(length(groupsnames)>2){
      allhypotheses<-hypotheses(length(numbersofreplicates)-1)
      allhypothesescompact<-numeric()
      for(l in 1:length(allhypotheses)){
        allhypothesescompact<-rbind(allhypothesescompact,allhypotheses[[l]])
      }
    allhypothesescompact<-unique(allhypothesescompact)


    #flag all hypotheses which have already been tested by assigning a p-value, else p-value=-9999
      flagpvalues<-matrix(-9999,nrow=nrow(allhypothesescompact),ncol=ncol(allhypothesescompact))
      results=list()
      pvalsCPCAT=vector()
      alreadyflaggedindex=vector()
      for(j in 1:(length(groupsnames)-1)){
        #identify contrast a p-value!=-9999 has been assigned to - these contrast must not be tested again
        contrasts<-hypotheses((ncol(ntype)-1))[[j]]
        matchingrows<-numeric()
        for(i in 1:nrow(contrasts)){
          matchingrows<-c(matchingrows,which(apply(allhypothesescompact,1,identical,contrasts[i,])))
        }
        alreadyflaggedindex<-which(flagpvalues[matchingrows,j]!=-9999)

        #shorten contrasts to be tested by elimantion of already tested contrasts
        if(length(alreadyflaggedindex)>0){
          contrasts<-contrasts[-alreadyflaggedindex,]
        }
        #in the last step the contrastmatrix reduces to a vector
        #make it a matrix consisting of nrow=1
        if(is.matrix(contrasts)==FALSE){
          contrasts<-matrix(contrasts,nrow=1)
        }
        notflaggedindex<-which(flagpvalues[matchingrows,j]==-9999)
        #flagpvalues which are still -9999
        #after CPCAT corresponding pvalues will be !=-9999
        tobeflagged<-matchingrows[notflaggedindex]

        results[[j]]<-poisson.test(z,contrastmatrix=contrasts,M=M)[[1]]
        #write obtained p-values into column j of flagpvalues and find max p-value
        pvalshelp<-results[[j]][,ncol(results[[j]])]
        flagpvalues[tobeflagged,j]<-pvalshelp
        #put together new p-values of reduced contrast matrix and relevant pvalues in flagpvalues[,j]
        if(j>1){	#in step j=1 all flagpvalues equal -9999
          pvalshelp2<-c(pvalshelp,flagpvalues[matchingrows[-notflaggedindex],j])
        }else{
          pvalshelp2<-pvalshelp
        }
        pvalsCPCAT[j]<-max(pvalshelp2)
        #copy p-values obtained so far to the next column of flagpvalues
        if(j<(length(groupsnames)-1)){
          flagpvalues[,j+1]<-flagpvalues[,j]
        }
        #assign flagpvalues corresponding to the hypotheses tested in step j
      }
      #if(any(groupsnames==999)){
      #  results=data.frame(Hypothesis=paste("H0: mu_",groupsnames[1],"=mu_",groupsnames[2],sep=""), pvalue=pvalsCPCAT[1])
      #}else{
        results=data.frame(Hypothesis=paste("H0: mu_",groupsnames[1],"=mu_",groupsnames[2:length(groupsnames)],sep=""), pvalue=pvalsCPCAT)
      #}
  }
    if(length(groupsnames)==2){
    allhypothesescompact=matrix(1,nrow=1, ncol=1)
      colnames(allhypothesescompact)="treatment 1"
      dat=list()
      IDsgroup1=which(z[,2]==groupsnames[1])
      IDsgroup2=which(z[,2]==groupsnames[2])
      dat[[1]]=z[IDsgroup1,]
      dat[[2]]=z[IDsgroup2,]
      pvalsCPCAT=poisson.sub.test(dat,contrast=1,M=M)[[1]]
      results=data.frame(Hypothesis=paste("H0: mu_",groupsnames[1],"=mu_",groupsnames[2],sep=""), pvalue=pvalsCPCAT)
    }
    structure(results)
  }
}
