
rm(list=ls())
library(phytools)
library(mvMORPH)
library(apTreeshape)
library(ade4)
library(geometry)
library(geiger)

# One of the functions ('evolveRateTre') below has been writen 
# by j. Beaulieu and can be found here: http://www.jeremybeaulieu.org/r.html

CreateSimpleTrees=function(n,type)
{
  if (type=="yule")
  {Tree=sim.bdtree(n=n)}
  
  if (type=="fullbal"){Tree=compute.brlen(stree(n=n,type="balanced"))
  #Tree=rescale(Tree,model='delta',delta=7.2)
  }
  
  if (type=="fullImbal")
  {Tree=compute.brlen(stree(n=n,type="left"))
  #Tree=rescale(Tree,model='delta',delta=.358)
  }
  
  if (type=="star")
  {
    Tree1=sim.bdtree(n=n)
    Tree=stree(n=n, type = "star", tip.label = Tree1$tip.label)
    Tree$edge.length=rep(1,n)
  }
  
  if (type=="coal")
  {Tree <- rcoal(n)}
  
  return(Tree)
}

SimpleSimulatedTreeAnalysis=function(e,params)
{
  print(e)
  # Simulate the Tree
  TreeType=as.character(params[e,'TreeType'])
  Tree=CreateSimpleTrees(n=params[e,'n'],type=TreeType)
  
  #Run the analysis
  Res=try(FDrdPdsets(m=params[e,"m"],Tree=Tree,nameTree=TreeType,model=as.character(params[e,"model"]),n.trt=params[e,"n.trt"],delta=params[e,"delta"],k=params[e,"k"],transfo=as.character(params[e,'transfo']),kappa=params[e,'kappa'],sigma=params[e,"sigma"],beta=params[e,"beta"],alpha=params[e,"alpha"],theta=params[e,"theta"]))
  if (class(Res)=="try-error") {Res=rep(NA,22)}
  
  return(Res)
}

TreeToLoadAnalysis=function(e,params,TreeFileAdress)
{
  print(e)
  # Load the Tree
  TreeType=as.character(params[e,'TreeType'])
  Tree=read.tree(paste(TreeFileAdress,TreeType,sep=''))

  #Run the analysis
  re=try(FDrdPdsets(m=params[e,"m"],Tree=Tree,nameTree=TreeType,model=as.character(params[e,"model"]),n.trt=params[e,"n.trt"],delta=params[e,"delta"],k=params[e,"k"],transfo=params[e,'transfo'],kappa=params[e,'kappa'],sigma=params[e,"sigma"],beta=params[e,"beta"],alpha=params[e,"alpha"],theta=params[e,"theta"]))
  if (class(re)=="try-error") {re=rep(NA,22)}
  
  return(re)
}


FDrdPdsets=function(m,Tree,nameTree='',model='BM1',n.trt=1,k,sigma=1,theta=0,delta,alpha,beta,transfo,kappa)
{

n=length(Tree$tip.label)

#################
# Handle Trees  #
#################
  
  #standardize Trees
  H=max(nodeHeights(Tree))
  Tree$edge.length=Tree$edge.length/H
  
  TreeIni=Tree #keep the initial phylo tree
  
  #Get Measures of TreeShape
  Ic=colless(as.treeshape(Tree))
  Gamma=gammaStat(Tree)
  Beta=maxlik.betasplit(as.treeshape(Tree))$max_lik
  
  #Tree transformations
  if (transfo=="kappa"){Tree=rescale(Tree,model='kappa',kappa=kappa)}
  #standardize Trees

  
#################
# Create Traits #
#################

  #Create BM/OU/EB Traits
  if (model%in%c('BM1','OU1','EB'))
  {
    opt <- list(ntraits=n.trt, sigma=diag(sigma,n.trt), alpha=diag(alpha,n.trt), beta=rep(beta,n.trt), theta=rep(theta,n.trt))
    Traits <- try(mvSIM(Tree, n=1, model=model, param=opt))
    if (class(Traits)=='try-error') {Traits=matrix(NA,ncol=n.trt,nrow=n)}
    if (n.trt>1) {colnames(Traits)=c(1:n.trt)} 
    rownames(Traits)=Tree$tip.label
  }
  
  #Create Discrete Traits
  if (model=="DISC_MC")
  {
    rates=matrix(sigma,ncol=k,nrow=k)
    diag(rates)=0
    t1=rTraitDisc(phy = Tree,model = rates,root.value = sample(x=1:k,size=1))
    Traits=data.frame(t1)
    if (n.trt>1)
    {for (j in 2:n.trt){Traits[[paste('t',j,sep="")]]=rTraitDisc(phy = Tree,model = rates,root.value =  sample(x=1:k,size=1))}}
  }
  
#################
#  Phy. signal  #
#################  
  
  # 1. Mantel
  if (model%in%c("DISC_MC","DISC_Thres"))
  {TraitDis=as.matrix(daisy(Traits,metric='gower'))} else {TraitDis=as.matrix(dist(Traits,method='euclidian'))}
  PhyDist=cophenetic(TreeIni)
  MantelSpearman=try(mantel(TraitDis,PhyDist[colnames(TraitDis),colnames(TraitDis)],method="spearman",perm=0)$statistic)
  if (class(MantelSpearman)=='try-error') {MantelSpearman=NA}
  
  # 2. Bloomberg K
  if (model%in%c("DISC_MC")) {BloomK=Lambda=NA} 
  if (!model%in%c("DISC_MC")) 
  {
    if (n.trt==1) 
    {BloomK<-phylosig(TreeIni,Traits)
    Lambda=phylosig(TreeIni,Traits,method="lambda")$lambda}
    else {
      BloomK=Lambda=c()
      for (f in 1:n.trt)
      {
      BloomK<-c(BloomK,phylosig(TreeIni,Traits[,f]))
      Lambda<-c(Lambda,phylosig(TreeIni,Traits[,f],method="lambda")$lambda)
      }
      BloomK=mean(BloomK)
      Lambda=mean(Lambda)
      }
  } 
    
############################
#  Create SETS of species  #
############################  
  
  setPD=GreedyMMD(TreeIni,m)   #Define ste of PD maximising species
  setRD=sample(Tree$tip.label,m) #Define set of random sp
  
  sitesp=matrix(0,n,ncol=2)
  rownames(sitesp)=Tree$tip.label
  colnames(sitesp)=c("Rd","PD")
  sitesp[setRD,"Rd"]=1
  sitesp[setPD,"PD"]=1
  sitesp=sitesp[rownames(Traits),]
  sitesp=data.frame(sitesp)
  
  
################
#  Compute FD  #
################
sp=rownames(Traits)

if (!model%in%c("DISC_MC","DISC_Thres"))
{Traits=apply(Traits,2,scale)
rownames(Traits)=sp} #resclaing of continuous traits


if (!model%in%c("DISC_MC","DISC_Thres"))  # Continuous Traits
  {
    # Convex Hull metric
    if ((n.trt>1)&(n.trt<m)&(sum(is.na(Traits[setPD,]))==0))
    {   
      rawFRic_PD=convhulln(Traits[setPD,],'FA')$vol
      rawFRic_RD=convhulln(Traits[setRD,],'FA')$vol
    } else if (n.trt==1){
      rawFRic_PD=max(Traits[setPD,])-min(Traits[setPD,])
      rawFRic_RD=max(Traits[setRD,])-min(Traits[setRD,])
    } else {rawFRic_PD=rawFRic_RD=NA}
    
    # Rao metric
    distTraits=dist(Traits)
    distTraitsSd=dist(Traits)/n.trt #following Botka (2005)
    
    div=divc(sitesp,distTraits)
    Rao_PD=div['PD',]
    Rao_RD=div['Rd',]
    
    div=divc(sitesp,distTraitsSd)
    RaoSd_PD=div['PD',]
    RaoSd_RD=div['Rd',]
  }
  
if (model%in%c("DISC_MC")){  # Discrete Traits
  if (n.trt>1)
    {   
    rawFRic_PD=nrow(unique(Traits[setPD,]))
    rawFRic_RD=nrow(unique(Traits[setRD,]))}
    else {
    rawFRic_PD=length(unique(Traits[setPD,]))
    rawFRic_RD=length(unique(Traits[setRD,]))
    }
  RaoSd_PD=RaoSd_RD=Rao_PD=Rao_RD=NA
  }
  
  FD=c(rawFRic_PD,rawFRic_RD,RaoSd_PD,RaoSd_RD,Rao_PD,Rao_RD,Ic,Gamma,Beta,BloomK,MantelSpearman,Lambda)
  names(FD)=c('FRich_PD','FRich_RD','RaoSd_PD','RaoSd_RD','Rao_PD','Rao_RD','Ic','Gamma','Beta','BloomK','MantelSpearman','Lambda')
  
  para=c(k=k,n=n,m=m,model=model,n.trt=n.trt,sigma=sigma,theta=theta,alpha=alpha,beta=beta,transfo=transfo,kappa=kappa,Tree=nameTree)
  FD=c(para,FD)
  return(FD)
}

#Greedy ALGO
GreedyMMD=function(Tree,k,tol = 1e-08)
{
  D=cophenetic(Tree)
  #1. Initialize 2 first species
  #set1=arrayInd(which.max(D),dim(D))[1:2]
  set1=arrayInd(which(abs(c(D)-max(D))<tol),dim(D)) #list of pairs
  set1=set1[sample(x=1:dim(set1)[1],1),] #randomly choose one
  set=colnames(D)[set1]
  c=2
  #2. Loop
  while (c<k)
  {
    #print(c)
    MinDist=apply(D[set,],2,min)
    #set=c(set,names(MinDist)[(MinDist==max(MinDist))][1])
    set2=names(MinDist)[abs(MinDist-max(MinDist))<tol] #find all other sp that would maximise PD
    set=c(set,sample(set2,1)) #randomly select 1
    c=c+1
  }
  return(set) 
}


#epsvec <- c(0, .25, .5, .75);
#lamvec <- c(.06437, 0.07834, 0.1018213, 0.1512680);
#sdvec <- c(0.001, seq(0.01, .06, by=0.01));
#REPS <- 2000;

evolveRateTree <- function (b, stdev, time.stop = 0, mintax = 3, maxtax = 5000) 
{
  
  if (time.stop == 0) 
    stop("Must have stopping criterion\n");
  
  breaker <- 0;
  while (1) {
    edge <- rbind(c(1, 2), c(1, 3));
    
    birth.time <- rep(0, maxtax);
    end.time <- rep(-1, maxtax);
    lambda <- rep(b, maxtax);
    edge.status <- rep(FALSE, maxtax);
    
    current <- 2;
    next.node <- 4;
    
    currvec <- current;
    
    while (sum(edge.status[1:nrow(edge)]) < nrow(edge)){
      
      current <- edge[,2][edge.status[1:nrow(edge)] == FALSE][1]; #new current
      t <- birth.time[1:nrow(edge)][edge[,2] == current];
      
      repeat{
        
        currvec <- c(currvec, current);
        
        dt <- rexp(1, lambda[1:nrow(edge)][edge[,2] == current]);
        t <- t + dt;
        print(t)
        if (t >= time.stop) {
          end.time[1:nrow(edge)][edge[,2] == current] <- time.stop;
          edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
          #cat('Tmax exceeded:', t, 'curr:', current, '\n', sep='\t');
          #print((1:nrow(edge))[edge[,2] == current]);
          #print(edge.status[edge[,2] == current]); 
          #print(edge.status[(1:nrow(edge))[edge[,2] == current]]);
          #print(edge.status[1:nrow(edge)][edge[,2] == current]) 
          break;
        }
        ##
        
        edge.status[1:nrow(edge)][edge[,2] == current] <- TRUE;
        
        edge <- rbind(edge, c(current, next.node), c(current, next.node + 1));
        
        birth.time[1:nrow(edge)][edge[,1] == current] <- t; 
        
        end.time[1:nrow(edge)][edge[,2] == current] <- t;
        
        ## evolving lambda:
        
        lambda[1:nrow(edge)][edge[,1] == current] <- rlnorm(1, log(lambda[1:nrow(edge)][edge[,2] == current]), sqrt(dt)*stdev);
        
        
        
        #print(end.time[1:nrow(edge)][edge[,2 ] == current]);
        
        #cat('curr:', current, 'next:', next.node, 'time:', t, '\n', sep='\t');
        #z <- data.frame(v1=edge[,1], v2=edge[,2], es=edge.status[1:nrow(edge)], bt=birth.time[1:nrow(edge)], et=end.time[1:nrow(edge)]);
        #print(z);
        
        
        current <- next.node;
        next.node <- next.node + 2;
        
        if (nrow(edge) >= maxtax)
          break;
        
        #breaker <- breaker + 1;
        #if (breaker > 5000)
        #	stop('breaker exceeded\n');
        
      }#repeat
      
      if (nrow(edge) >= maxtax){
        print('maxtax exceeded');
        break;	
      }
      
      
    }#while (sum(edge.status))
    
    birth.time <- birth.time[1:nrow(edge)];
    end.time <- end.time[1:nrow(edge)];
    lambda <- lambda[1:nrow(edge)];
    
    edge.length <- end.time - birth.time;
    
    if (nrow(edge) >= 2*mintax & nrow(edge) < maxtax)
      break;
    
  }# while (1)
  
  
  n <- -1
  for (i in 1:max(edge)) {
    if (any(edge[, 1] == i)) {
      edge[which(edge[, 1] == i), 1] <- n;
      edge[which(edge[, 2] == i), 2] <- n;
      n <- n - 1;
    }
  }
  edge[edge > 0] <- 1:sum(edge > 0);
  tip.label <- 1:sum(edge > 0);
  mode(edge) <- "character";
  mode(tip.label) <- "character";
  obj <- list(edge = edge, edge.length = edge.length, tip.label = tip.label, birth.time= birth.time, end.time = end.time, lambda=lambda, currvec = currvec);
  class(obj) <- "phylo";
  obj <- old2new.phylo(obj);
  obj;
  
}



