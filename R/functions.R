library(scales)
library(bipartite)
library(mipfp)


zscore <- function(x){(x-mean(x))/sd(x)}
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

tile_vector=function(IN=net,var,by_col=F){
  fill=IN;fill[]=NA
  if(by_col){fill=t(fill)}
  
  for(i in 1:nrow(fill)){
    fill[i,]=var  
  }
  if(by_col){fill=t(fill)}
  return(as.matrix(fill))
}

generate_interaction_frequcies=function(var,n_int_tot){
  sam=sample(1:length(var),size = n_int_tot,replace = T,prob = var)
  tab=table(sam)
  
  out=rep(0,length(var))
  nam=as.numeric(names(tab))
  out[nam]=tab
  return(out)
}


networklevel_sim=function(IN=net,NI=tot_int,pred_mat_freq,pred_mat_bin,null_model=NULL){
  #flatten matrics
  unit=as.vector(unlist(pred_mat_freq));unit[]=1
  freq=as.vector(unlist(pred_mat_freq))
  bin=as.vector(unlist(pred_mat_bin))
  
  ###combining morphological matching with the null model
  if(!is.null(null_model)){
    freq=freq/sum(freq)
    null_model=null_model/sum(null_model)
    freq=freq*null_model
  }
  
  freq
  
  ###adding zeros from the binary component
  zeros=rep(1,length(freq))
  sam=rbinom(length(zeros), 1,bin)
  freq[which(sam==0)]=0
  

  
  #simulating interactions from the frequency component
  ints=rep(0,length(freq))
  sam=sample(1:length(freq),size = NI,prob = freq,replace = T)
  tab=table(sam)
  ints[as.numeric(names(tab))]=tab
  
  result=net
  result[]=ints
  return(result)
}


hummingbird_foraging_sim=function(IN=net,NIps=rep(100,ncol(net)),pred_mat_freq,pred_mat_bin){
  
  ###adding zeros from the binary component
  zero=as.matrix(IN);zero[]=1
  for(i in 1:nrow(IN)){
    sam=rbinom(8, 1,  as.vector(unlist(pred_mat_bin[i,])))
    which(sam==0)
    zero[i,which(sam==0)]=0
  }
  pred_mat_freq=pred_mat_freq*zero
  
  
  #simulating interactions from the frequency component
  zero=as.matrix(IN);zero[]=0
  for(i in 1:ncol(IN)){
    sam=sample(1:nrow(IN),size = NIps[i],prob = as.vector(pred_mat_freq[,i]),replace = T)
    tab=table(sam)
    zero[as.numeric(names(tab)),i]=tab
  }
  
  
  return(zero)
}


simulate_matrix <- function(h = NULL, p = NULL, pred_mat_freq, pred_mat_bin) {
  tol = 1e-20; j=p;k=h
  
  freq=as.vector(unlist(pred_mat_freq))
  bin=as.vector(unlist(pred_mat_bin))
  
  sam=rbinom(length(bin), 1,bin)
  freq[which(sam==0)]=0
  
  W=pred_mat_freq;W[]=freq
  
  
  # Check that at least one marginal is provided:
  if (is.null(j) && is.null(k)) {
    stop("At least one of j (row sums) or k (column sums) must be provided.")
  }
  
  # Check dimensions of W versus provided marginals:
  if (!is.null(j)) {
    if (length(j) != nrow(W))
      stop("Length of j (row sums) must equal the number of rows in W.")
  }
  if (!is.null(k)) {
    if (length(k) != ncol(W))
      stop("Length of k (column sums) must equal the number of columns in W.")
  }
  
  # Determine the total sum n from the provided marginal(s)
  if (!is.null(j) && !is.null(k)) {
    if (sum(j) != sum(k))
      stop("When both j and k are provided, their sums must be equal.")
    n <- sum(j)
  } else if (!is.null(j)) {
    n <- sum(j)
  } else {  # only k is provided
    n <- sum(k)
  }
  
  # Create a seed matrix from W scaled to total n.
  seed <- n * W
  
  # Load the mipfp package (install if necessary)
  if (!requireNamespace("mipfp", quietly = TRUE)) {
    stop("Package 'mipfp' is required. Please install it using install.packages('mipfp').")
  }
  library(mipfp)
  
  # Set up the target lists based on which marginals are provided.
  target.list <- list()
  target.data <- list()
  if (!is.null(j)) {
    target.list <- c(target.list, list(1))
    target.data <- c(target.data, list(j))
  }
  if (!is.null(k)) {
    target.list <- c(target.list, list(2))
    target.data <- c(target.data, list(k))
  }
  
  # Use IPF to adjust the seed so that the provided marginals are matched.
  result <- Ipfp(seed, target.list, target.data, tol = tol,iter=10000)
  continuous_table <- result$x.hat
  
  # Start by flooring the continuous solution.
  integer_table <- floor(continuous_table)
  
  # Recompute the fractional parts.
  fractional_part <- continuous_table - integer_table
  
  # Define deficits for each provided margin.
  if (!is.null(j)) {
    row_deficit <- j - rowSums(integer_table)
  }
  if (!is.null(k)) {
    col_deficit <- k - colSums(integer_table)
  }
  
  # Greedy allocation: while there is any deficit in the constrained margins.
  # We'll update the fractional_part after each allocation.
  # When both margins are constrained, both deficits must be positive for eligibility.
  # When only one margin is constrained, only that margin is used.
  allocate_extra <- function() {
    # Build an eligibility matrix based on provided constraints.
    eligible <- matrix(TRUE, nrow = nrow(integer_table), ncol = ncol(integer_table))
    if (!is.null(j)) {
      for (i in seq_len(nrow(integer_table))) {
        if (row_deficit[i] <= 0) eligible[i, ] <- FALSE
      }
    }
    if (!is.null(k)) {
      for (jcol in seq_len(ncol(integer_table))) {
        if (col_deficit[jcol] <= 0) eligible[, jcol] <- FALSE
      }
    }
    # Set non-eligible cells to -Inf so they won't be chosen.
    temp <- fractional_part
    temp[!eligible] <- -Inf
    return(which(temp == max(temp), arr.ind = TRUE))
  }
  
  # Decide which deficit to check in the while loop:
  if (!is.null(j) && !is.null(k)) {
    total_deficit <- sum(row_deficit)  # same as sum(col_deficit)
  } else if (!is.null(j)) {
    total_deficit <- sum(row_deficit)
  } else {  # only k is provided
    total_deficit <- sum(col_deficit)
  }
  
  while (total_deficit > 0) {
    idx <- allocate_extra()
    # In case of ties, pick the first.
    i_sel <- idx[1, 1]
    j_sel <- idx[1, 2]
    
    # Add one count to the selected cell.
    integer_table[i_sel, j_sel] <- integer_table[i_sel, j_sel] + 1
    
    # Update the fractional part for that cell.
    fractional_part[i_sel, j_sel] <- continuous_table[i_sel, j_sel] - integer_table[i_sel, j_sel]
    
    # Update deficits.
    if (!is.null(j)) {
      row_deficit[i_sel] <- row_deficit[i_sel] - 1
    }
    if (!is.null(k)) {
      col_deficit[j_sel] <- col_deficit[j_sel] - 1
    }
    
    # Recompute the overall deficit.
    if (!is.null(j) && !is.null(k)) {
      total_deficit <- sum(row_deficit)
    } else if (!is.null(j)) {
      total_deficit <- sum(row_deficit)
    } else {
      total_deficit <- sum(col_deficit)
    }
  }
  
  return(integer_table)
}


