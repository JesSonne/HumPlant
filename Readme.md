HumPlant R package
================
2025-02-14

- [Visualising community structure aming interacting hummingbirds and
  plants](#visualising-community-structure-aming-interacting-hummingbirds-and-plants)
  - [Installing the package](#installing-the-package)
  - [Taking a look at the data](#taking-a-look-at-the-data)
  - [Defining models for how morphologies influence the species’
    interactions](#defining-models-for-how-morphologies-influence-the-species-interactions)
  - [Simulate an interaction network based on the two morphological
    models](#simulate-an-interaction-network-based-on-the-two-morphological-models)
  - [Simulate an interaction network based on the species’
    abundances.](#simulate-an-interaction-network-based-on-the-species-abundances)
  - [Combining models based on abundances and
    morphologies](#combining-models-based-on-abundances-and-morphologies)

## Visualising community structure aming interacting hummingbirds and plants

The following script will guide you through examples of how to model
mechanisms, such as morphological matching, and visualise their effects
on species interactions.

### Installing the package

You will use the functions from the ‘HumPlant’ R package, which is
hosted on this GitHub page. We use another package called ‘devtools’ to
install it, which we first have to install and load in your R session.
Remove the hashtag to enable the installation.

``` r
#install.packages('devtools') #you only need to run this once
require(devtools)
```

I will refine the content throughout the course with your feedback and
suggestions. Therefore, you must reinstall the package using the code
below after each update. Remove the hashtag to enable the installation.

``` r
#devtools::install_github("JesSonne/HumPlant")
```

Now load the package in your current R session.

``` r
library(HumPlant)
```

### Taking a look at the data

The object ‘Cajanuma’ contains the data I collected at the
high-elevation site in Southern Ecuador. It contains the interaction
network along with the ecological attributes we will use for the
modelling.

Below, we define a series of objects with the data we are going to use.
Notice that the birds and plants are sorted according to the length of
their bill/flower

``` r
#loading the interaction matrix
net=Cajanuma$Network

#loading ecological attributes
hum_morph=Cajanuma$Hummingbird_morphologies
plant_morph=Cajanuma$Plant_morphologies

hum_abund=Cajanuma$Hummingbird_abundances
plant_abund=Cajanuma$Plant_abundances

plant_phenol=Cajanuma$Plant_phenologies

#number of hummingbirds and plants in the comunity
n_hums=ncol(net)
n_plants=nrow(net)
```

Use the ‘plotweb’ function to visualise the network. Notice that the
species are sorted according to the length of their bill/flowers from
the shortest bills/flowers on the left to the longest on the right.

``` r
plotweb(net,method="normal",empty = F)
```

![](Readme_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

The following code converts or ‘tile’ the vectors on the ecological
attributes into a matrix format. This is just to make them easier to
work with.

``` r
hum_morph_mat=tile_vector(net,hum_morph)
plant_morph_mat=tile_vector(net,plant_morph,by_col = T)

hum_abund_mat=tile_vector(net,hum_abund)
plant_abund_mat=tile_vector(net,plant_abund,by_col = T)

hum_abund_mat=tile_vector(var=hum_abund)
plant_abund_mat=tile_vector(var=plant_abund,by_col = T)
```

### Defining models for how morphologies influence the species’ interactions

The code below defines two initial models for how the morphologies of
species might influence the species’ interactions. These are just two
examples to get you started. In the project, you should come up with
additional models that perhaps better capture the structure of the
community. I will help you write the code as long as you can
specifically formulate your hypotheses.

The first is based on the ‘forbidden links’ concept in which species are
prevented from visiting certain partners in the community. In this case,
hummingbirds are prevented from visiting flowers that are longer than
their bill + tounge length. The tounge length is defined as a fraction
of the total bill length.

The second model is derived from optimal foraging theory, where
hummingbirds should prefer visiting flowers that are morphologically
similar to their bills. Here, we simply subtract the two measurements
from each other and take the absolute value. 0 will reflect a perfect
match, and higher values reflect increasing mismatch. For convenience,
we take inverse, such that higher values represent greater morphological
matching.

``` r
barrier_long_flowers=function(h=hum_morph_mat,p=plant_morph_mat,tounge=1.8){
  m=h*tounge-p
  
  b=m;b[]=1
  b[which(m<0)]=0
  
  return(as.matrix(b))
}

absdif=function(h=hum_morph_mat,p=plant_morph_mat,tounge=1.8){
  m=1/(abs(h*tounge-p))
  return(as.matrix(m))
}
```

### Simulate an interaction network based on the two morphological models

Start by stating the number of interactions you simulated per
hummingbird species. To start with, let’s assume they have the same
total number of interactions.

``` r
n_hum_int=rep(100,n_hums) 
```

Now, generate matrices of morphological matching. Notice that the two
models use different tongue lengths. Why do you think that is, and does
it make sense?

``` r
matching_matrix=absdif(tounge=1)
barrier_matrix=barrier_long_flowers(tounge=1.8)
```

The function below simulates an artificial network based on your input.
The purpose is to directly visualise how different hypotheses/models
affect interactions in the network. The function takes several input
files, where some can be left empty (denoted ‘NULL’).

The simulation itself is a two-stage process.

- First, we simulate whether a pair of species have any or no
  Interactions. As such, this is the binary part of the simulation. In
  this example, we use the model of the morphological barrier.

- Secondly, we simulate how often a pair of species interact, should
  they have any interactions with each other. As such, this is the
  quantitative part of the simulation. In this example, we use the model
  of the morphological matching.

When the simulation concludes, it generates two indices of network
structure. One index measures the degree of network nestedness, while
the other measures the degree of specialisation. Use these two indices
to discuss the implications of your hypotheses for niche differentiation
in the community.

NB. In the code, we make no specific assumptions about how many
interactions the plants have, but we could. Be aware that when
specifying the number of interactions for both plants and hummingbirds,
they must sum up to the same number (i.e. the total number of
interactions in the network). If ‘p’ and ‘h’ are both NULL, you must
specify the total number of interactions in the network ‘n’.

``` r
sim_net_morph=simulate_ZI_matrix(
                           #A vector stating the total number of interactions for each plant species
                           p = NULL,                     
                           #A vector stating the total number of interactions for each hummingibird species
                           h = n_hum_int,                
                           #Alternatively state the total number of interactions in the network 
                           n = NULL,                     
                           #A matrix of probabilities for two individual species to have any interactions
                           W_bin=barrier_matrix,         
                           #A matrix of weights proportional the species interaction frequencies
                           W_freq=matching_matrix,       
                           ) 
```

    ##                   Nestedness Complementary specialization 
    ##                   20.0164723                    0.3150969

Now, let’s take a look at the simulated network. Does the structure
coincide with your expectations and what could be improved?

``` r
plotweb(sim_net_morph,method="normal",empty = F)
```

![](Readme_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### Simulate an interaction network based on the species’ abundances.

Common species are more easily encountered than rare species. If species
interact by random encounters, the interaction probability between two
species (i and j) equals their joint ecounter probability = P(i) x P
(j). We expect the encounter probabilities of species to be strongly
correlated with their abundances. Therefore, an intuitive null model
would assume the number of interactions between two species are
proportional to the product of their abundances. We will use it for the
quantitative part of the next simulation.

``` r
abundance_model=function(h=hum_abund_mat,p=plant_abund_mat){
  m=h*p
  return(m)
}
```

Under the assumption of random interactions, there are no processes
preventing two species from having any interactions. Therefore, for the
binary part of the simulation, we use a matrix stating that all
interactions are equally likely.

``` r
unit=net;unit[]=1
```

Now, we are ready to simulate a network without the influence of species
morphologies where species interact as they randomly encounter each
other.

NB. In the code, we now make no specific assumptions about how many
interactions the plants and hummingbirds have. We only state the total
number of interactions we want simulated ‘n’.

``` r
sim_net_abund=simulate_ZI_matrix(
                           #A vector stating the total number of interactions for each plant species
                           p = NULL,                     
                           #A vector stating the total number of interactions for each hummingibird species
                           h = NULL,                
                           #Alternatively state the total number of interactions in the network 
                           n = 500,                     
                           #A matrix of probabilities for two individual species to have any interactions 
                           W_bin=unit,         
                           #A matrix of weights proportional the species interaction frequencies
                           W_freq=abundance_model(),       
                           ) 
```

    ##                   Nestedness Complementary specialization 
    ##                     56.39239                      0.00000

Let’s plot it! Do you recognize this structure? For clarity, we orde the
species in the network according to their abundance. The rarest left,
and the most common species is on the right.

``` r
sim_net_abund=sim_net_abund[order(plant_abund),]
sim_net_abund=sim_net_abund[,order(hum_abund)]

plotweb(sim_net_abund,method="normal",empty = F)
```

![](Readme_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

### Combining models based on abundances and morphologies

We can combine the models based on morphologies and abundances and
visualize how they combined influence the structure of the network.

We just need to tell the function how to aggregate the weight matrices.
For most applications, it makes the most sense to use the product

``` r
sim_net_morph_abund=simulate_ZI_matrix(
                           #A vector stating the total number of interactions for each plant species
                           p = NULL,                     
                           #A vector stating the total number of interactions for each hummingibird species
                           h = NULL,                
                           #Alternatively state the total number of interactions in the network 
                           n = 500,                     
                           #A matrix of probabilities for two individual species to have any interactions 
                           W_bin=barrier_matrix,         
                           #A matrix of weights proportional the species interaction frequencies
                           W_freq=list(abundance_model(),matching_matrix),
                           #how to combine the weight matrices for the binary part of the simulation
                           comb_method_bin = "product",
                           #how to combine the weight matrices for the quantitative part of the simulation
                           comb_method_freq =  "product",
                           #should the matrices be normalised before they are combined - always use 'TRUE'
                           normalize = TRUE,
                           ) 
```

    ##                   Nestedness Complementary specialization 
    ##                   32.1723331                    0.5278625

``` r
plotweb(sim_net_morph_abund,method="normal",empty = F)
```

![](Readme_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->
