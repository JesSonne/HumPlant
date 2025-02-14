# Helper function: combine one or more weight matrices/data.frames.
combine_weights <- function(W_input, comb_method = c("average", "product"), normalize = TRUE) {
  comb_method <- match.arg(comb_method)
  
  # If W_input is a data.frame, convert it to a matrix.
  if (is.data.frame(W_input)) {
    W_input <- as.matrix(W_input)
  }
  
  # If W_input is a list, convert any data.frames to matrices.
  if (is.list(W_input)) {
    W_input <- lapply(W_input, function(mat) {
      if (is.data.frame(mat)) {
        mat <- as.matrix(mat)
      }
      if (normalize) {
        mat / sum(mat)
      } else {
        mat
      }
    })
    
    # Check that all matrices have the same dimensions.
    dims <- sapply(W_input, dim)
    if (!all(apply(dims, 1, function(x) length(unique(x)) == 1))) {
      stop("All weight matrices must have the same dimensions.")
    }
    
    if (comb_method == "average") {
      combined <- Reduce("+", W_input) / length(W_input)
    } else if (comb_method == "product") {
      combined <- Reduce("*", W_input)
      # Optionally, take the geometric mean:
      # combined <- combined^(1/length(W_input))
    }
  } else {
    # If a single weight matrix/data.frame is provided:
    if (is.data.frame(W_input)) {
      W_input <- as.matrix(W_input)
    }
    combined <- if (normalize) W_input / sum(W_input) else W_input
  }
  return(combined)
}

# Function to simulate an integer matrix based on a combined weight matrix.
simulate_matrix <- function(p = NULL, h = NULL, n = NULL, 
                            W, comb_method = c("average", "product"),
                            normalize = TRUE,
                            tol = 1e-8) {
  comb_method <- match.arg(comb_method)
  
  # Accept W as a data.frame, matrix, or list of them.
  if (is.data.frame(W)) {
    W <- as.matrix(W)
  }
  
  if (is.list(W)) {
    W <- lapply(W, function(mat) {
      if (is.data.frame(mat)) mat <- as.matrix(mat)
      if (normalize) mat / sum(mat) else mat
    })
    dims <- sapply(W, dim)
    if (!all(apply(dims, 1, function(x) length(unique(x)) == 1))) {
      stop("All weight matrices must have the same dimensions.")
    }
    if (comb_method == "average") {
      combined_W <- Reduce("+", W) / length(W)
    } else if (comb_method == "product") {
      combined_W <- Reduce("*", W)
    }
  } else {
    combined_W <- if (normalize) W / sum(W) else W
  }
  
  # Check dimensions against provided marginals.
  if (!is.null(p) && length(p) != nrow(combined_W))
    stop("Length of p (row sums) must equal the number of rows in the weight matrix.")
  if (!is.null(h) && length(h) != ncol(combined_W))
    stop("Length of h (column sums) must equal the number of columns in the weight matrix.")
  
  # Determine overall total.
  if (!is.null(p) && !is.null(h)) {
    if (sum(p) != sum(h))
      stop("When both p and h are provided, their sums must be equal.")
    n_total <- sum(p)
  } else if (!is.null(p)) {
    n_total <- sum(p)
  } else if (!is.null(h)) {
    n_total <- sum(h)
  } else {
    if (is.null(n))
      stop("Provide at least one of p, h, or the overall total n.")
    n_total <- n
  }
  
  # Create seed matrix by scaling combined weights.
  seed <- n_total * combined_W
  
  # If any marginal constraints are provided, adjust the seed using IPF.
  if (!is.null(p) || !is.null(h)) {
    if (!requireNamespace("mipfp", quietly = TRUE)) {
      stop("Package 'mipfp' is required. Please install it using install.packages('mipfp').")
    }
    library(mipfp)
    
    target.list <- list()
    target.data <- list()
    if (!is.null(p)) {
      target.list <- c(target.list, list(1))
      target.data <- c(target.data, list(p))
    }
    if (!is.null(h)) {
      target.list <- c(target.list, list(2))
      target.data <- c(target.data, list(h))
    }
    
    result <- Ipfp(seed, target.list, target.data, tol = tol,iter = 100000)
    continuous_table <- result$x.hat
  } else {
    continuous_table <- seed
  }
  
  # Rounding: floor and then adjust by fractional parts.
  integer_table <- floor(continuous_table)
  fractional_part <- continuous_table - integer_table
  
  if (!is.null(p) || !is.null(h)) {
    if (!is.null(p)) row_deficit <- p - rowSums(integer_table)
    if (!is.null(h)) col_deficit <- h - colSums(integer_table)
    
    total_deficit <- if (!is.null(p) && !is.null(h)) sum(row_deficit)
    else if (!is.null(p)) sum(row_deficit)
    else sum(col_deficit)
    
    while (total_deficit > 0) {
      eligible <- matrix(TRUE, nrow = nrow(integer_table), ncol = ncol(integer_table))
      if (!is.null(p)) {
        for (i in seq_len(nrow(integer_table))) {
          if (row_deficit[i] <= 0) eligible[i, ] <- FALSE
        }
      }
      if (!is.null(h)) {
        for (jcol in seq_len(ncol(integer_table))) {
          if (col_deficit[jcol] <= 0) eligible[, jcol] <- FALSE
        }
      }
      temp <- fractional_part
      temp[!eligible] <- -Inf
      idx <- which(temp == max(temp), arr.ind = TRUE)
      i_sel <- idx[1, 1]
      j_sel <- idx[1, 2]
      
      integer_table[i_sel, j_sel] <- integer_table[i_sel, j_sel] + 1
      fractional_part[i_sel, j_sel] <- continuous_table[i_sel, j_sel] - integer_table[i_sel, j_sel]
      
      if (!is.null(p)) row_deficit[i_sel] <- row_deficit[i_sel] - 1
      if (!is.null(h)) col_deficit[j_sel] <- col_deficit[j_sel] - 1
      
      total_deficit <- if (!is.null(p) && !is.null(h)) sum(row_deficit)
      else if (!is.null(p)) sum(row_deficit)
      else sum(col_deficit)
    }
  } else {
    current_total <- sum(integer_table)
    total_deficit <- n_total - current_total
    while (total_deficit > 0) {
      temp <- fractional_part
      idx <- which(temp == max(temp), arr.ind = TRUE)
      i_sel <- idx[1, 1]
      j_sel <- idx[1, 2]
      
      integer_table[i_sel, j_sel] <- integer_table[i_sel, j_sel] + 1
      fractional_part[i_sel, j_sel] <- continuous_table[i_sel, j_sel] - integer_table[i_sel, j_sel]
      total_deficit <- n_total - sum(integer_table)
    }
  }
  
  return(integer_table)
}

# Now define the zero-inflated simulation function.
simulate_ZI_matrix <- function(p = NULL, h = NULL, n = NULL, 
                               W_bin, W_freq,
                               comb_method_bin = c("average", "product"),
                               comb_method_freq = c("average", "product"),
                               normalize = TRUE,
                               tol = 1e-8) {
  comb_method_bin <- match.arg(comb_method_bin)
  comb_method_freq <- match.arg(comb_method_freq)
  
  # Allow W_bin and W_freq to be data.frames, matrices, or lists thereof.
  combined_W_bin <- combine_weights(W_bin, comb_method_bin, normalize)
  combined_W_freq <- combine_weights(W_freq, comb_method_freq, normalize)
  
  # Check dimensions.
  nr <- nrow(combined_W_bin)
  nc <- ncol(combined_W_bin)
  if (nrow(combined_W_freq) != nr || ncol(combined_W_freq) != nc) {
    stop("The binary and frequency weight matrices must have the same dimensions.")
  }
  
  # For the binary component, convert combined weights into probabilities.
  # If the weight matrix is uniform (e.g., all ones), normalization yields a uniform probability matrix.
  prob_bin <- combined_W_bin / max(combined_W_bin)
  
  # Simulate the binary indicator matrix.
  M_bin <- matrix(rbinom(nr * nc, size = 1, prob = prob_bin), nrow = nr, ncol = nc)
  
  # Simulate the frequency component.
  M_freq <- simulate_matrix(p = p, h = h, n = n, W = combined_W_freq, 
                            comb_method = comb_method_freq, normalize = normalize, tol = tol)
  
  # Combine: force cells with a binary indicator of 0 to be zero.
  M_final <- M_freq * M_bin
  
  index=networklevel(M_final,c("H2","weighted NODF"))
  names(index)=c("Nestedness","Complementary specialization")
  print(index)
  
  return(M_final)
}
