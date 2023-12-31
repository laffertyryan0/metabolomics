setwd(".")

#read all fixed things from a json file init.json

params = rjson::fromJSON(file="./data/init.json")

cov_files = list.files(params$covariate_files)
K = length(cov_files)


cov_list = list()
for(i in 1:K){
  cov_list[[i]] = read.csv2(paste(params$covariate_files,"/",cov_files[i],sep=""),
                            sep=",",header=FALSE)
  cov_list[[i]] = matrix(data = as.numeric(unname(as.matrix(cov_list[[i]]))),
                         nrow = dim(unname(as.matrix(cov_list[[i]])))[1])
}

N = dim(cov_list[[1]])[1] #labs

P = dim(cov_list[[1]])[2] #covariates

K = length(cov_list) #metabolites

covariates = array(0,dim=c(N,P,K))
for(k in 1:K){
  covariates[,,k] = cov_list[[k]]
}

beta = params$true_beta

b = params$true_b

pr = array(rep(0,N*K),dim = c(N,K))
missing = array(rep(0,N*K),dim = c(N,K))
for(lab_num in 1:N){
  for(metabolite in 1:K){
    pr[lab_num,metabolite] = pnorm(-1 + covariates[lab_num,,metabolite]%*%beta + b[lab_num])
    missing[lab_num,metabolite] = rbinom(1,1,pr[lab_num,metabolite])
  }
}

#For each metabolite, randomly choose one lab and force it to be non-missing
#The model is still approximately correct
#this makes the mets.in.some.lab.idx no longer necessary, but we will keep it 
#for now
for(metabolite in 1:K){
  random_lab = sample(1:N,1)
  missing[random_lab,metabolite] = 0
}

mets.in.some.lab.idx = (colSums(1-missing) > 0)  #boolean index showing TRUE where metabolite exists
                                                    #in at least one lab

data = data.frame(y = c(t(missing)), 
                  lab_ind = seq(1,N)%x%rep(1,K))  #this is not the metabolite data, but the data for fitting missingness model
              #    X1 = c(t(covariates[,1,]))#,
              #    X2 = c(t(covariates[,2,])),
               #   X3 = c(t(covariates[,3,])),
                #  X4 = c(t(covariates[,4,]))
                 # )

varsumstr = ""
for(i in 1:P){
  data = cbind(data,c(t(covariates[,i,])))
  varname = paste("X",i,sep = "")
  names(data)[i+2] = varname
  varsumstr = paste(varsumstr,varname,"+")
}

# y~X1 + X2 +.. +X_P + (1|lab_ind)

mod<-lme4::glmer(as.formula(paste("y ~",varsumstr," (1|lab_ind)")), 
           data = data, family=binomial("probit"))

samples.per.lab = params$samples_per_lab
mean = params$mean
corr = do.call("cbind",params$corr)
cov = diag(K)%*%corr%*%diag(K)
max.samp = max(samples.per.lab)

met.data = list()
X = list()
lam = list()
for(i in 1:N){
  met.data[[i]] = mvtnorm::rmvnorm(samples.per.lab[i],
                          mean=mean,
                          sigma = cov)
  
  #ignore data for metabolites that are missing in every lab
  met.data[[i]] = met.data[[i]][,mets.in.some.lab.idx]
  
  X[[i]] = cor(met.data[[i]],method = "spearman")
  lam[[i]] = 1/samples.per.lab[i]
  
}


K.dropped = sum(mets.in.some.lab.idx)

p = array(0,dim=c(N,K,K))
delta = array(0,dim=c(N,K,K))
bihat = c(lme4::ranef(mod)$lab_ind[[1]])
for(i in 1:N){
  betas = coef(mod)[1]$lab_ind[i,]
  
  one.minus.pr.est.vec.i = 1-pnorm(as.double(
    bihat[i]  + t(covariates[i,,])%*%as.double(betas[-1])))
  #p is the probability of being non-missing
  p[i,,] = outer(one.minus.pr.est.vec.i,one.minus.pr.est.vec.i)
  diag(p[i,,]) = sqrt(diag(p[i,,]))
  delta[i,,] = outer(1-missing[i,],1-missing[i,])
}


delta.dropped = list()
p.dropped = list()

for(i in 1:N){
  delta.dropped[[i]] = delta[i,mets.in.some.lab.idx,mets.in.some.lab.idx]
  p.dropped[[i]] = p[i,mets.in.some.lab.idx,mets.in.some.lab.idx]
}


X.tilde = matrix(rep(0,K.dropped*K.dropped),nrow=K.dropped)
H = matrix(rep(0,K.dropped*K.dropped),nrow=K.dropped)

for(j1 in 1:K.dropped){
  for(j2 in 1:K.dropped){
    sum.x.num = 0
    sum.x.denom = 0
    sum.2.h = 0
    for(i in 1:N){
      #Only contribute to the sum if probability of non-missing is > 0
      if(p.dropped[[i]][j1,j2]>0){
        sum.x.num = sum.x.num + 
          lam[[i]]*delta.dropped[[i]][j1,j2]*X[[i]][j1,j2]/p.dropped[[i]][j1,j2]
        sum.x.denom = sum.x.denom + 
          lam[[i]]*delta.dropped[[i]][j1,j2]/p.dropped[[i]][j1,j2]
        sum.2.h = sum.2.h + 
          lam[[i]]*delta.dropped[[i]][j1,j2]/p.dropped[[i]][j1,j2]
      }
      
    }
    if(sum.x.denom!=0){  
      X.tilde[j1,j2] = sum.x.num/sum.x.denom 
    }
    else{  #in case a lab has no metabolites, it shouldn't contribute to the sum
      X.tilde[j1,j2] = 0
    }
    
    H[j1,j2] = sqrt(sum.2.h)
  }
}

#print(X.tilde)
#print(H)

stacked.met.data = do.call(rbind,met.data)

tryCatch({
write.table(x=mets.in.some.lab.idx+0,file="./mets_in_some_lab_index.csv",sep=",",row.names = F, col.names = F)
write.table(x=stacked.met.data,file="./stacked_met_data.csv",sep=",",row.names = F,col.names = F)
write.table(x=X.tilde,file="./x_matrix.csv",sep = ",",row.names = F,col.names = F)
write.table(x=H,file="./h_matrix.csv",sep = ",",row.names = F,col.names = F)

write.table(x=t(c(pr)),file="./pr_vecs.csv",sep = ",",row.names = F,col.names = F,append=T)
write.table(x=t(c(missing)),file="./missing_vecs.csv",sep = ",",row.names = F,col.names = F,append=T)
write.table(x=min(eigen(X.tilde)$values),file="./min_eigs.csv",sep = ",",row.names = F,col.names = F,append=T)
},
error = function(e){
  print(e$message)
})


print(min(eigen(X.tilde)$values))
