###################################################
#         SCRIPT TO RUN the ANALYSIS 
###################################################

MC=16 # Numbers of cores to parralellize the computations

#--------------------------------------------
# Estimate of FD for classic simple trees (Table 1)
#--------------------------------------------

# Parameters 
nsp=c(32,64)
repet=1000
TreeType=c("yule",'coal')
savedsp=c(8,16)
sig=1
n.trt=c(1,2,4)
theta=0
k=0

#Simple BM 
alpha=0
beta=0
params <- expand.grid(model="BM1",n=nsp,k=k,transfo='no', rep=1:repet, m=savedsp,TreeType=TreeType,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=NA)

#Simple OU 
alpha=c(1.4,7)
params <- rbind(params,expand.grid(model='OU1',n=nsp,k=k,transfo='no', rep=1:repet, m=savedsp,TreeType=TreeType,theta=theta, alpha=alpha, sigma=sig, beta=0, n.trt=n.trt,kappa=NA))

#Simple EB
alpha=0
beta=c(-5,-1)
params <- rbind(params,expand.grid(model='EB',n=nsp,k=k,transfo='no', rep=1:repet, m=savedsp,TreeType=TreeType,theta=theta, alpha=0, sigma=sig, beta=beta, n.trt=n.trt,kappa=NA))

#MK model
k=4
params <- rbind(params,expand.grid(model="DISC_MC",n=nsp,k=k,transfo='no', rep=1:repet, m=savedsp,TreeType=TreeType,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=NA))

head(params)
nn=dim(params)[1]

#Runs
FDall=t(mcmapply(FUN=SimpleSimulatedTreeAnalysis,1:nn,MoreArgs = list(params),mc.cores=MC))
Res=data.frame(FDall)
indx <- colnames(Res)[!colnames(Res)%in%c('Tree','model','transfo')]
Res[indx] <- lapply(Res[indx], function(x) as.numeric(as.character(x)))

#Saving
save(Res,file='Res_Table1_FINAL_2.Rdata')


#----------------------------------------------------
# Estimate of FD for more imbalanced trees (Figure 1)
#----------------------------------------------------

# 'Archetypal' imbalanced and balanced Trees
#-------------------------------------------
MC=30
# Parameters 
nsp=c(32,64)
repet=1000
TreeType=c("fullImbal",'fullbal')
savedsp=c(8,16)
sig=1
n.trt=c(1,2,4)
theta=0
k=0
kappa=c(seq(from=0,to=1,by=.1))
transfo='kappa'

#Simple BM 
alpha=0
beta=0
params <- expand.grid(model="BM1",n=nsp,k=k,transfo=transfo, rep=1:repet, m=savedsp,TreeType=TreeType,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=kappa)
nn=dim(params)[1]

#MK model
savedsp=c(8)
nsp=c(64)
k=4
sig=c(.1,.5,1)
repet=1000
TreeType=c("fullImbal",'fullbal')
sig=1
n.trt=c(1,2,4)
params <- rbind(params,expand.grid(model="DISC_MC",n=nsp,k=k,transfo=transfo, rep=1:repet, m=savedsp,TreeType=TreeType,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=kappa))
nn=dim(params)[1]

#Runs
MC=16
FDall=t(mcmapply(FUN=SimpleSimulatedTreeAnalysis,1:nn,MoreArgs = list(params),mc.cores=MC))
Res=data.frame(FDall)
indx <- colnames(Res)[!colnames(Res)%in%c('Tree','model','transfo')]
Res[indx] <- lapply(Res[indx], function(x) as.numeric(as.character(x)))

save(Res,file='Results_ArchetypalTrees_Kappa_Discrete_Continuous_1_2_4traits_32_64sp_8_16savedBIS2.Rdata')


#MK model
savedsp=c(8)
nsp=c(64)
k=4
sig=c(.1,.5,1)
params <- rbind(expand.grid(model="DISC_MC",n=nsp,k=k,transfo=transfo, rep=1:repet, m=savedsp,TreeType=TreeType,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=kappa))
nn=dim(params)[1]
MC=16
FDall=t(mcmapply(FUN=SimpleSimulatedTreeAnalysis,1:nn,MoreArgs = list(params),mc.cores=MC))
Res=data.frame(FDall)
indx <- colnames(Res)[!colnames(Res)%in%c('Tree','model','transfo')]
Res[indx] <- lapply(Res[indx], function(x) as.numeric(as.character(x)))
save(Res,file='Results_ArchetypalTrees_Kappa_Discrete_1_2_4traits_32_64sp_8_16saved.Rdata')



# Intermediate imbalanced and balanced Trees
#-------------------------------------------

# We first produced the trees 
# using J. Beaulieu scripts that can be found here http://www.jeremybeaulieu.org/r.html
# They come from this publication: Beaulieu, J.M., and B.C. O'Meara (2015). Extinction can be estimated from moderately sized molecular phylogenies. Evolution. 
MC=30

# 64 species

nsp=64
repet=100
sig=1
n.trt=c(1,2,4)
alpha=0
beta=0
k=0
TreeFileAdress='/home/mazel/SCAP/BeaulieuSelec/'
savedsp=c(8,16) 
transfo='kappa'
kappa=c(seq(from=0,to=1,by=.1))
listTree=list.files(TreeFileAdress)
theta=0

#Simple BM 
params <- expand.grid(model="BM1",n=nsp,k=k,transfo=transfo, rep=1:repet, m=savedsp,TreeType=listTree,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=kappa)

#MK model
savedsp=c(8)
k=4
n.trt=c(2)
sig=c(.1,1)
params <- rbind(params,expand.grid(model="DISC_MC",n=nsp,k=k,transfo=transfo, rep=1:repet, m=savedsp,TreeType=listTree,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=kappa))

nn=dim(params)[1] # 792 000 simulations
FDall=t(mcmapply(FUN=TreeToLoadAnalysis,1:nn,MoreArgs = list(params,TreeFileAdress),mc.cores=MC))
save(FDall,file='Raw_IntermediateBalance_64sp_1_2_4_traits_m8_16.Rdata')
Res=data.frame(FDall)
indx <- colnames(Res)[!colnames(Res)%in%c('Tree','model','transfo')]
Res[indx] <- lapply(Res[indx], function(x) as.numeric(as.character(x)))
save(Res,file='IntermediateBalance_64sp_1_2_4_traits_m8_16.Rdata')

# 32 species
nsp=32
repet=100
sig=1
n.trt=c(1,2,4)
alpha=0
beta=0
k=0
TreeFileAdress='/home/mazel/SCAP/BeaulieuSelec32/'
savedsp=c(8,16) 
transfo='kappa'
kappa=c(seq(from=0,to=1,by=.1))
listTree=list.files(TreeFileAdress)
theta=0

#Simple BM 
params <- expand.grid(model="BM1",n=nsp,k=k,transfo=transfo, rep=1:repet, m=savedsp,TreeType=listTree,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=kappa)

#MK model
savedsp=c(8)
k=4
n.trt=c(2)
sig=c(.1,1)
params <- rbind(params,expand.grid(model="DISC_MC",n=nsp,k=k,transfo=transfo, rep=1:repet, m=savedsp,TreeType=listTree,theta=theta, alpha=0, sigma=sig, beta=0, n.trt=n.trt,kappa=kappa))

nn=dim(params)[1] # 792 000 simulations
FDall=t(mcmapply(FUN=TreeToLoadAnalysis,1:nn,MoreArgs = list(params,TreeFileAdress),mc.cores=MC))
save(FDall,file='Raw_IntermediateBalance_32sp_1_2_4_traits_m8_16.Rdata')
Res=data.frame(FDall)
indx <- colnames(Res)[!colnames(Res)%in%c('Tree','model','transfo')]
Res[indx] <- lapply(Res[indx], function(x) as.numeric(as.character(x)))
save(Res,file='IntermediateBalance_32sp_1_2_4_traits_m8_16.Rdata')


