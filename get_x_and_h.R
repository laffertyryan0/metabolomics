setwd("/home/ryan/projects/metabolomics")

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

data = data.frame(y = c(t(missing)), 
                  lab_ind = seq(1,N)%x%rep(1,K))
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
  X[[i]] = cor(met.data[[i]],method = "spearman")
  lam[[i]] = 1/samples.per.lab[i]
}


delta = list()
p = list()
bihat = c(lme4::ranef(mod)$lab_ind[[1]])
for(i in 1:N){
  betas = coef(mod)[1]$lab_ind[i,]
  delta[[i]] = matrix(rep(0,K*K),nrow=K)
  p[[i]] = matrix(rep(0,K*K),nrow=K)
  for(j1 in 1:K){
    for(j2 in 1:K){
      pr.est1 = pnorm(as.double(
        bihat[i]  + covariates[i,,j1]%*%as.double(betas[-1])))
      pr.est2 = pnorm(as.double(
        bihat[i]  + covariates[i,,j2]%*%as.double(betas[-1])))
                       
      if(j1 == j2){
        p[[i]][j1,j2] = 1 - pr.est1
        delta[[i]][j1,j2] = as.integer((1 - missing[i,j1]))
      }
      else{
        delta[[i]][j1,j2] = as.integer((1 - missing[i,j1])&(1-missing[i,j2]))
        p[[i]][j1,j2] = (1 - pr.est1)*(1 - pr.est2)
      }
    }
  }
}

X.tilde = matrix(rep(0,K*K),nrow=K)
H = matrix(rep(0,K*K),nrow=K)
for(j1 in 1:K){
  for(j2 in 1:K){
    sum.x.num = 0
    sum.x.denom = 0
    sum.2.h = 0
    for(i in 1:N){
      sum.x.num = sum.x.num + 
        lam[[i]]*delta[[i]][j1,j2]*X[[i]][j1,j2]/p[[i]][j1,j2]
      sum.x.denom = sum.x.denom + 
        lam[[i]]*delta[[i]][j1,j2]/p[[i]][j1,j2]
      sum.2.h = sum.2.h + 
        lam[[i]]*delta[[i]][j1,j2]/p[[i]][j1,j2]
    }
    X.tilde[j1,j2] = sum.x.num/sum.x.denom
    H[j1,j2] = sqrt(sum.2.h)
  }
}

#print(X.tilde)
#print(H)

write.table(x=X.tilde,file="x_matrix.csv",sep = ",",row.names = F,col.names = F)
write.table(x=H,file="h_matrix.csv",sep = ",",row.names = F,col.names = F)

write.table(x=t(c(pr)),file="pr_vecs.csv",sep = ",",row.names = F,col.names = F,append=T)
write.table(x=t(c(missing)),file="missing_vecs.csv",sep = ",",row.names = F,col.names = F,append=T)
write.table(x=min(eigen(X.tilde)$values),file="min_eigs.csv",sep = ",",row.names = F,col.names = F,append=T)

print(min(eigen(X.tilde)$values))
