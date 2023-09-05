library(ConsensusClusterPlus)

select_optimal_k <- function(results,max_k){
  # result : ConsensusClustering output
  # max_k : integer
  Kvec = 2:max_k
  x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC = rep(NA,length(Kvec)) 
  names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
  for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }#end for i
  # The optimal K
  optK = Kvec[which.min(PAC)]
  return(optK)
}

random_consensus <- function(exp_m,n_target,max_k,consensus_path){
  # exp_m : matrix; row is feature, col is sample
  # n_target : integer; number of random sampling features
  # max_k : integer; max number of cluster number
  # consensus_path : output path of consensus clustering result
  random_genes = sample(rownames(exp_m),n_target)
  sub_m = exp_m[random_genes,]
  results = ConsensusClusterPlus(sub_m,maxK=max_k,reps=50,pItem=0.8,pFeature=1,title=consensus_path,clusterAlg="pam",distance="euclidean",seed=1262118388.71279,plot='png')
  optimal_k = select_optimal_k(results,max_k)
  label = results[[optimal_k]]$consensusClass
  return(label)
}

probiotic_consensus <- function(exp_m,n_target,consensus_path,output_path,max_k =10){
  # exp_m : matrix; row is feature, col is sample
  # n_target : integer; number of consensus clustering
  # max_k : integer; max number of cluster number
  # consensus_path : output path of consensus clustering result
  dir.create(consensus_path, showWarnings = FALSE)

  results = ConsensusClusterPlus(exp_m,maxK=max_k,reps=50,pItem=0.8,pFeature=1,
                                 title=consensus_path,clusterAlg="pam",distance="euclidean",seed=1262118388.71279,plot='png')
  label = data.frame(cluster = results[[n_target]]$consensusClass)
  write.table(label,output_path,sep='\t',quote = F)
}