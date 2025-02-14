setwd('/Users/hzc973/Library/Mobile Documents/com~apple~CloudDocs/Undervisning/Ornithology 2025/Project') #desktop
setwd("/Users/jespersonne/Library/Mobile Documents/com~apple~CloudDocs/Undervisning/Ornithology 2025/Project") #laptop

require(bipartite)
require(scales)
require(mipfp)

####loading some homemade functions from another script
source("functions.R")

####loading the interaction matrix
net=read.csv("Network/data/Cajanuma/interactionmatrix.csv")[,-1]
n_hums=ncol(net)
n_plants=nrow(net)

####loading ecological attributes
hum_morph=read.csv("Network/data/Cajanuma/birdmorphology.csv")[,-1]
plant_morph=read.csv("Network/data/Cajanuma/plantmorphology.csv")[,-1]

hum_abund=read.csv("Network/data/Cajanuma/birdabundance.csv")[,-1]
plant_abund=read.csv("Network/data/Cajanuma/plantabundance.csv")[,-1]

plant_phenol=read.csv("Network/data/Cajanuma/plantphenology.csv")[,-1]
plant_phenol=rowSums(plant_phenol)

##### For convenience, let's order the species according to their trait sizes. NB all data files have to follow the same order (i.e. according to their morphologies)
pl_ord=order(plant_morph) #order of plant traits
hum_ord=order(hum_morph) #ordier of bird traits

net=net[pl_ord,] #reordering rows
net=net[,hum_ord] #reordering columns

#reordering morphologies
hum_morph=hum_morph[hum_ord]
plant_morph=plant_morph[pl_ord]

#reordering abundances
hum_abund=hum_abund[hum_ord]
plant_abund=plant_abund[pl_ord]



#####'spashing' vectors into matrix format to make them easier to work with
hum_morph_mat=splash_network(net,hum_morph)
plant_morph_mat=splash_network(net,plant_morph,by_col = T)

hum_abund_mat=splash_network(net,hum_abund)
plant_abund_mat=splash_network(net,plant_abund,by_col = T)

hum_abund_mat=splash_network(var=hum_abund)
plant_abund_mat=splash_network(var=plant_abund,by_col = T)

### Now, write here your models for morphological matching
absdif=function(h=hum_morph_mat,p=plant_morph_mat,tounge=1.8){
  m=1/(abs(h*tounge-p))
  return(as.matrix(m))
}

barrier_long_flowers=function(h=hum_morph_mat,p=plant_morph_mat,tounge=1.8){
  m=h*tounge-p
  
  b=m;b[]=1
  b[which(m<0)]=0
  
  return(as.matrix(b))
}

barrier_short_flowers=function(h=hum_morph_mat,p=plant_morph_mat,prop=0.2){
  m=p/h
  
  b=m;b[]=1
  b[which(m<prop)]=0
  
  return(as.matrix(b))
}

compt_load=function(exp=1){1/splash_network(IN=net,var,by_col = T)^exp}
  
###Write a model that represents your null hypothesis. Here is a suggestion based on the species abundances
abundance_model=function(h=hum_abund_mat,p=plant_abund_mat){
  m=h*p
  return(m)
}


######simulate interaction networks based on morphological matching

#state here number of interactions simulated per hummingbird species
n_hum_int=rep(100,n_hums) # assuming hummingbirds have the same total number of interactions. For now let us not assume anything regarding the plants' total number of interactions


#generating matrices of morphological matching
matching_matrix=absdif(tounge=1.8)
barrier_matrix=barrier_long_flowers(tounge=2)

#simulate a network
sim_net=simulate_ZI_matrix(p = NULL, 
                           h = n_hum_int, 
                           n = NULL, 
                           W_bin=barrier_matrix,
                           W_freq=matching_matrix,
                           comb_method_bin = "product",
                           comb_method_freq =  "product",
                           normalize = TRUE,
                           tol = 1e-20) 


#plot it
plotweb(sim_net,method="normal",empty = F)


#### alternative models based on species abundances

#generating matrices of without morphological matching - all interactions are equally likely morphologically speaking
unit=net;unit[]=1
abundance_matrix=abundance_model()


#simulate a network
sim_net=simulate_ZI_matrix(p = NULL, 
                           h = NULL, 
                           n = 500, 
                           W_bin=unit,
                           W_freq=abundance_matrix,
                           comb_method_bin = "product",
                           comb_method_freq =  "product",
                           normalize = TRUE,
                           tol = 1e-8) 

#plot it
plotweb(sim_net,method="normal",empty = F)

#### for clarity, lets order the network based on the species abundances
sim_net=sim_net[order(plant_abund),]
sim_net=sim_net[,order(hum_abund)]

#plot it again
plotweb(sim_net,method="normal",empty = F)




#### combining models based on abundances and morphologies
#state here number of interactions simulated per hummingbird species
sim_net=simulate_ZI_matrix(p = NULL, 
                           h = NULL, 
                           n = 500, 
                           W_bin=barrier_matrix,
                           W_freq=list(abundance_matrix,matching_matrix),
                           comb_method_bin = "product",
                           comb_method_freq =  "product",
                           normalize = TRUE,
                           tol = 1e-8) 

#plot it
plotweb(sim_net,method="normal",empty = F)

#compare it to the empirical matrix
plotweb(net,method="normal",empty = F)






