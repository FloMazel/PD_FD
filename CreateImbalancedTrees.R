library(ape)
library(geiger)
library(apTreeshape)

# LOAD R function by J. Beaulieu 
path="/home/" # where to save the trees

b=0.06587368
eps=0
stdev=0.06
NN=64 # Nb of wanted species
seqs=rev(seq(from=-2,to=2,by=0.4))
intNb=length(seqs)-1
for (i in 1:intNb)
{
  su=0
  while (su<10)
  {
    try(mm<-evolveRateTree.eps(b=b,eps=eps,stdev=stdev,maxtax=5000),silent=TRUE)
    try(Tree<-drop.extinct(mm),silent=TRUE)
    ntips=length(Tree$tip.label)
    print(ntips)
    if (ntips>NN) #we want tree with 64 species only
    {
      todrop=length(Tree$tip.label)-NN
      Tree=drop.tip(Tree,sample(Tree$tip.label,todrop))
      try(Beta<-maxlik.betasplit(as.treeshape(Tree))$max_lik)
      
      if ((Beta<seqs[i])&(Beta>seqs[i+1])) {
        su=su+1
        write.tree(Tree,file=paste(path,"BeaulieuTREE_No",su,"GammaInf",seqs[i],"_Sup",seqs[i+1],".trees",sep=''),append=TRUE)
      }
    }
 
  }
  print(i)
}

