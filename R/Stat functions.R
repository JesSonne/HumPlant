#' @export
fitinteractionModel <- function(response, countPredictors = list(), ziPredictors = list()) {
  # Check that the response matrix has row and column names
  if (!is.matrix(response) || is.null(rownames(response)) || is.null(colnames(response))) {
    stop("Response matrix must be a matrix with both row and column names.")
  }
  
  # Helper function to check predictor matrices
  check_predictors <- function(predictors, predictorType) {
    if (!is.list(predictors)) {
      stop(paste(predictorType, "predictors must be provided as a list of matrices."))
    }
    for (i in seq_along(predictors)) {
      pred_mat <- predictors[[i]]
      if (!is.matrix(pred_mat) || is.null(rownames(pred_mat)) || is.null(colnames(pred_mat))) {
        stop(paste("Each", predictorType, "predictor must be a matrix with row and column names."))
      }
      if (!all(dim(pred_mat) == dim(response))) {
        stop(paste("Each", predictorType, "predictor matrix must have the same dimensions as the response matrix."))
      }
    }
  }
  
  check_predictors(countPredictors, "count")
  check_predictors(ziPredictors, "zero-inflation")
  
  # Load required libraries
  library(glmmTMB)
  library(reshape2)
  
  # Melt the response matrix to long format with species identifiers
  df <- melt(response, varnames = c("plant", "hummingbird"), value.name = "count")
  
  # Process count predictor matrices
  if (length(countPredictors) > 0) {
    # Use provided names or assign default names
    countNames <- names(countPredictors)
    if (is.null(countNames) || any(countNames == "")) {
      countNames <- paste0("countPred", seq_along(countPredictors))
    }
    for (i in seq_along(countPredictors)) {
      temp <- melt(countPredictors[[i]], 
                   varnames = c("plant", "hummingbird"), 
                   value.name = countNames[i])
      df <- merge(df, temp, by = c("plant", "hummingbird"))
    }
    count_formula_vars <- paste(countNames, collapse = " + ")
  } else {
    count_formula_vars <- "1"
  }
  
  # Process zero-inflation predictor matrices
  if (length(ziPredictors) > 0) {
    ziNames <- names(ziPredictors)
    if (is.null(ziNames) || any(ziNames == "")) {
      ziNames <- paste0("ziPred", seq_along(ziPredictors))
    }
    for (i in seq_along(ziPredictors)) {
      temp <- melt(ziPredictors[[i]], 
                   varnames = c("plant", "hummingbird"), 
                   value.name = ziNames[i])
      df <- merge(df, temp, by = c("plant", "hummingbird"))
    }
    zi_formula_vars <- paste(ziNames, collapse = " + ")
  } else {
    zi_formula_vars <- "1"
  }
  
  # Construct the model formulas
  # Count (frequency) component: fixed effects from countPredictors and random effects for plants and hummingbirds.
  count_formula <- as.formula(paste("count ~", count_formula_vars,
                                    "+ (1 | plant) + (1 | hummingbird)"))
  
  # Zero-inflation (binary) component formula: fixed effects from ziPredictors
  zi_formula <- as.formula(paste("~", zi_formula_vars))
  
  # Fit the zero-inflated negative binomial model
  model <- glmmTMB(count_formula,
                   ziformula = zi_formula,
                   family = nbinom1,
                   data = df)
  
  # Return the fitted model and the merged data frame used for modeling
  return(list(model = model, data = df))
}



#####
#####
#####
#####


#' @export
EvaluatePredictorCombinations <- function(response, countCandidates = list(), ziCandidates = list(),maxminiter=1000) {
  # --- Check response matrix ---
  if (!is.matrix(response) || is.null(rownames(response)) || is.null(colnames(response))) {
    stop("Response matrix must be a matrix with both row and column names.")
  }
  
  # --- Helper: Convert 2D candidate to 3D array if necessary ---
  convert_to_array <- function(candidate) {
    if (is.matrix(candidate)) {
      candidate <- array(candidate, dim = c(nrow(candidate), ncol(candidate), 1))
    }
    if (!(is.array(candidate) && length(dim(candidate)) == 3)) {
      stop("Each candidate predictor must be either a 2D matrix or a 3D array with dimensions (h, p, i).")
    }
    candidate
  }
  
  # Convert candidates if needed and check dimensions.
  process_candidates <- function(candidates, predictorType) {
    if (!is.list(candidates)) {
      stop(paste(predictorType, "predictors must be provided as a list."))
    }
    for (i in seq_along(candidates)) {
      candidates[[i]] <- convert_to_array(candidates[[i]])
      if (!all(dim(candidates[[i]])[1:2] == dim(response))) {
        stop(paste("Each", predictorType, "predictor must have the same row and column dimensions as the response matrix."))
      }
    }
    candidates
  }
  
  hum_names=colnames(net);hummingbird_names=rownames(net)
  
  countCandidates <- process_candidates(countCandidates, "count")
  ziCandidates <- process_candidates(ziCandidates, "zero-inflation")
  
  # --- Load required libraries ---
  library(glmmTMB)
  library(reshape2)
  
  # --- Prepare the base long-format response data ---
  df_base <- melt(response, varnames = c("plant", "hummingbird"), value.name = "count")
  
  # --- Helper: Generate all subsets (including empty set) from a list of candidate arrays ---
  # For each candidate in the subset, generate all possible choices of its variant (i dimension).
  # Each returned combination is a list with:
  #   predictors: a list of matrices (one slice per candidate),
  #   variants: a named list with the chosen variant for each candidate,
  #   candidateIndices: the indices of the candidates selected.
  get_array_subsets <- function(candList) {
    n <- length(candList)
    subsets <- list()
    # Empty set (intercept-only)
    subsets[[1]] <- list(predictors = list(), variants = list(), candidateIndices = integer(0))
    count_index <- 2
    if (n > 0) {
      for (k in 1:n) {
        idx_combinations <- combn(seq_len(n), k, simplify = FALSE)
        for (idx_set in idx_combinations) {
          # For the chosen subset, list the available variant numbers for each candidate.
          variant_options <- lapply(idx_set, function(i) {
            seq_len(dim(candList[[i]])[3])
          })
          # Cartesian product of variant choices.
          grid <- expand.grid(variant_options, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
          for (row in seq_len(nrow(grid))) {
            chosen_variants <- as.numeric(grid[row,])
            pred_list <- list()
            var_list <- list()
            # Use candidate names if provided; otherwise, assign defaults.
            cand_names <- names(candList)[idx_set]
            if (is.null(cand_names) || any(cand_names == "")) {
              cand_names <- paste0("Pred", idx_set)
            }
            for (j in seq_along(idx_set)) {
              cand_index <- idx_set[j]
              variant_index <- chosen_variants[j]
              pred_mat <- candList[[cand_index]][, , variant_index]
              pred_list[[ cand_names[j] ]] <- pred_mat
              var_list[[ cand_names[j] ]] <- variant_index
            }
            subsets[[count_index]] <- list(predictors = pred_list,
                                           variants = var_list,
                                           candidateIndices = idx_set)
            count_index <- count_index + 1
          }
        }
      }
    }
    return(subsets)
  }
  
  # Generate all subsets for count and zero-inflation predictors.
  count_subsets <- get_array_subsets(countCandidates)
  zi_subsets <- get_array_subsets(ziCandidates)
  
  # --- Get names for all candidate predictors ---
  count_all_names <- names(countCandidates)
  if (is.null(count_all_names) || any(count_all_names == "")) {
    count_all_names <- paste0("Pred", seq_along(countCandidates))
  }
  zi_all_names <- names(ziCandidates)
  if (is.null(zi_all_names) || any(zi_all_names == "")) {
    zi_all_names <- paste0("Pred", seq_along(ziCandidates))
  }
  
  # --- Loop over all combinations of count and zero-inflation predictor subsets ---
  results <- list()
  for (i in seq_along(count_subsets)) {
    for (j in seq_along(zi_subsets)) {
      # Start with the base response data.
      df <- df_base
      
      # --- Process count predictors ---
      count_combo <- count_subsets[[i]]
      if (length(count_combo$predictors) > 0) {
        for (pred_name in names(count_combo$predictors)) {
          
          add=count_combo$predictors[[pred_name]]
          rownames(add)=hummingbird_names;colnames(add)=hum_names
          
          temp <- melt(add,
                       varnames = c("plant", "hummingbird"),
                       value.name = pred_name)
          df <- merge(df, temp, by = c("plant", "hummingbird"))
        }
        count_formula_vars <- paste(names(count_combo$predictors), collapse = " + ")
      } else {
        count_formula_vars <- "1"
      }
      
      # --- Process zero-inflation predictors ---
      zi_combo <- zi_subsets[[j]]
      if (length(zi_combo$predictors) > 0) {
        for (pred_name in names(zi_combo$predictors)) {
          add=zi_combo$predictors[[pred_name]]
          rownames(add)=hummingbird_names;colnames(add)=hum_names
          
          temp <- melt(add,
                       varnames = c("plant", "hummingbird"),
                       value.name = pred_name)
          df <- merge(df, temp, by = c("plant", "hummingbird"))
        }
        zi_formula_vars <- paste(names(zi_combo$predictors), collapse = " + ")
      } else {
        zi_formula_vars <- "1"
      }
      
      # --- Construct model formulas ---
      count_formula <- as.formula(paste("count ~", count_formula_vars,
                                        "+ (1 | plant) + (1 | hummingbird)"))
      zi_formula <- as.formula(paste("~", zi_formula_vars))
      
      # --- Fit the model and extract AIC ---
      model_fit <- tryCatch({
        mod <- glmmTMB(count_formula,
                       ziformula = zi_formula,
                       family = nbinom1,
                       data = df,
                       control = glmmTMBControl(optCtrl = list(iter.max = maxminiter)))
        AIC(mod)
      }, error = function(e) NA)
      
      # --- Prepare output: one column per candidate predictor ---
      # For count predictors:
      count_variant_out <- setNames(rep("", length(count_all_names)), paste0("count_", count_all_names))
      if (length(count_combo$variants) > 0) {
        for (nm in names(count_combo$variants)) {
          count_variant_out[paste0("count_", nm)] <- count_combo$variants[[nm]]
        }
      }
      # For zero-inflation predictors:
      zi_variant_out <- setNames(rep("", length(zi_all_names)), paste0("zi_", zi_all_names))
      if (length(zi_combo$variants) > 0) {
        for (nm in names(zi_combo$variants)) {
          zi_variant_out[paste0("zi_", nm)] <- zi_combo$variants[[nm]]
        }
      }
      
      # Save the result for this combination.
      results[[length(results) + 1]] <- c(list(AIC = model_fit),
                                          count_variant_out,
                                          zi_variant_out)
    }
  }
  
  # Convert results list to a data frame and remove entries with NA AIC.
  results_df <- do.call(rbind, lapply(results, as.data.frame, stringsAsFactors = FALSE))
  results_df <- results_df[!is.na(results_df$AIC), ]
  
  # Sort results by AIC (lowest on top) and compute delta AIC.
  results_df <- results_df[order(results_df$AIC), ]
  
  
  best_AIC <- min(results_df$AIC)
  results_df$delta_AIC <- results_df$AIC - best_AIC
  
  
  # Reset row names and return.
  rownames(results_df) <- NULL
  
  AICs=results_df[,c(1,ncol(results_df))]
  results_df=results_df[,-c(1,ncol(results_df))]
  results_df[] <- lapply(results_df, as.numeric)
  results_df=cbind(results_df,AICs)
  
  return(results_df)
} 

#' @export
repeat_function <- function(func, vary = list(), ...) {
  # func: the function to be called repeatedly
  # times: number of iterations
  # vary: a named list of vectors; each vector should have 'times' elements
  # ...: additional arguments to be passed to func (constant across iterations)
  
  # Pre-allocate a list to store each output matrix
  
  times=length(vary[[1]])
  results <- vector("list", times)
  
  for (i in seq_len(times)) {
    # Start with constant arguments
    args <- list(...)
    # Update arguments with the i-th element from each varying parameter
    for (name in names(vary)) {
      args[[name]] <- vary[[name]][i]
    }
    # Call the function with these arguments
    results[[i]] <- do.call(func, args)
  }
  
  # Combine the matrices along the third dimension into a 3D array.
  # Here we use the abind package. Make sure it's installed.
  if (!requireNamespace("abind", quietly = TRUE)) {
    stop("Package 'abind' is required but not installed. Install it using install.packages('abind').")
  }
  arr3d <- abind::abind(results, along = 3)
  
  return(arr3d)
}

#' @export
aggregate_by_category <- function(df, cat_col, cont_col,func=min) {
  # Check if the specified columns exist in the data frame
  if (!all(c(cat_col, cont_col) %in% names(df))) {
    stop("One or both specified columns do not exist in the data frame.")
  }
  
  # Use aggregate to compute the minimum of the continuous variable for each category
  result <- aggregate(df[[cont_col]], 
                      by = list(Category = df[[cat_col]]), 
                      FUN = func)
  
  # Rename the resulting column for clarity
  names(result)[1] <- paste(cat_col)
  names(result)[2] <- paste(cont_col)
  return(result)
}


#' @export
GenerateDataFrames <- function(response, countCandidates = list(), ziCandidates = list()) {
  # --- Check response matrix ---
  if (!is.matrix(response) || is.null(rownames(response)) || is.null(colnames(response))) {
    stop("Response matrix must be a matrix with both row and column names.")
  }
  
  # --- Helper: Convert 2D candidate to 3D array if necessary ---
  convert_to_array <- function(candidate) {
    if (is.matrix(candidate)) {
      candidate <- array(candidate, dim = c(nrow(candidate), ncol(candidate), 1))
    }
    if (!(is.array(candidate) && length(dim(candidate)) == 3)) {
      stop("Each candidate predictor must be either a 2D matrix or a 3D array with dimensions (h, p, i).")
    }
    candidate
  }
  
  # Process candidates: check class, convert if needed, and check dimensions.
  process_candidates <- function(candidates, predictorType) {
    if (!is.list(candidates)) {
      stop(paste(predictorType, "predictors must be provided as a list."))
    }
    for (i in seq_along(candidates)) {
      candidates[[i]] <- convert_to_array(candidates[[i]])
      if (!all(dim(candidates[[i]])[1:2] == dim(response))) {
        stop(paste("Each", predictorType, "predictor must have the same row and column dimensions as the response matrix."))
      }
    }
    candidates
  }
  
  countCandidates <- process_candidates(countCandidates, "count")
  ziCandidates <- process_candidates(ziCandidates, "zero-inflation")
  
  # --- Use the response matrix dimnames for merging ---
  plant_names <- rownames(response)
  hummingbird_names <- colnames(response)
  
  # --- Load required library ---
  library(reshape2)
  
  # --- Prepare the base long-format response data ---
  df_base <- melt(response, varnames = c("plant", "hummingbird"), value.name = "n_interactions")
  
  # --- Modified Helper: Generate only full combinations for a list of candidate arrays ---
  # For each candidate predictor, we generate all possible choices of its variant (i dimension),
  # but we only return combinations that include ALL predictors.
  get_array_subsets <- function(candList) {
    n <- length(candList)
    subsets <- list()
    if (n == 0) {
      # No predictors supplied: return an "empty" set (e.g., intercept only)
      subsets[[1]] <- list(predictors = list(), variants = list(), candidateIndices = integer(0))
    } else {
      idx_set <- seq_len(n)  # Use all candidate predictors
      # For each candidate, get available variant options
      variant_options <- lapply(idx_set, function(i) {
        seq_len(dim(candList[[i]])[3])
      })
      # Cartesian product of variant choices.
      grid <- expand.grid(variant_options, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      count_index <- 1
      for (row in seq_len(nrow(grid))) {
        chosen_variants <- as.numeric(grid[row,])
        pred_list <- list()
        var_list <- list()
        # Use candidate names if provided; otherwise assign default names.
        cand_names <- names(candList)
        if (is.null(cand_names) || any(cand_names == "")) {
          cand_names <- paste0("Pred", idx_set)
        }
        for (j in seq_along(idx_set)) {
          cand_index <- idx_set[j]
          variant_index <- chosen_variants[j]
          pred_mat <- candList[[cand_index]][, , variant_index]
          pred_list[[ cand_names[j] ]] <- pred_mat
          var_list[[ cand_names[j] ]] <- variant_index
        }
        subsets[[count_index]] <- list(predictors = pred_list,
                                       variants = var_list,
                                       candidateIndices = idx_set)
        count_index <- count_index + 1
      }
    }
    return(subsets)
  }
  
  # Generate full combinations for count and zero-inflation predictors.
  count_subsets <- get_array_subsets(countCandidates)
  zi_subsets <- get_array_subsets(ziCandidates)
  
  # --- Get names for all candidate predictors (for later metadata) ---
  count_all_names <- names(countCandidates)
  if (is.null(count_all_names) || any(count_all_names == "")) {
    count_all_names <- paste0("Pred", seq_along(countCandidates))
  }
  zi_all_names <- names(ziCandidates)
  if (is.null(zi_all_names) || any(zi_all_names == "")) {
    zi_all_names <- paste0("Pred", seq_along(ziCandidates))
  }
  
  # --- Loop over all combinations of count and zero-inflation predictor subsets ---
  # Here, for each combination we merge the predictors with the base response data.
  result_dfs <- list()
  for (i in seq_along(count_subsets)) {
    for (j in seq_along(zi_subsets)) {
      # Start with the base response data.
      df <- df_base
      
      # --- Process count predictors ---
      count_combo <- count_subsets[[i]]
      if (length(count_combo$predictors) > 0) {
        for (pred_name in names(count_combo$predictors)) {
          add <- count_combo$predictors[[pred_name]]
          # Set row and column names using response names.
          rownames(add) <- plant_names
          colnames(add) <- hummingbird_names
          temp <- melt(add, varnames = c("plant", "hummingbird"), value.name = pred_name)
          df <- merge(df, temp, by = c("plant", "hummingbird"))
        }
      }
      
      # --- Process zero-inflation predictors ---
      zi_combo <- zi_subsets[[j]]
      if (length(zi_combo$predictors) > 0) {
        for (pred_name in names(zi_combo$predictors)) {
          add <- zi_combo$predictors[[pred_name]]
          rownames(add) <- plant_names
          colnames(add) <- hummingbird_names
          temp <- melt(add, varnames = c("plant", "hummingbird"), value.name = pred_name)
          df <- merge(df, temp, by = c("plant", "hummingbird"))
        }
      }
      
      # --- Optionally, add metadata columns to record the combination ---
      if (length(count_combo$variants) > 0) {
        for (nm in names(count_combo$variants)) {
          df[[paste0("count_variant_", nm)]] <- count_combo$variants[[nm]]
        }
      }
      if (length(zi_combo$variants) > 0) {
        for (nm in names(zi_combo$variants)) {
          df[[paste0("zi_variant_", nm)]] <- zi_combo$variants[[nm]]
        }
      }
      
      # Save the data frame for this combination.
      result_dfs[[length(result_dfs) + 1]] <- df
    }
  }
  
  return(result_dfs)
}


# Generate all subsets for count and zero-inflation predictors.
count_subsets <- get_array_subsets(countCandidates)
zi_subsets <- get_array_subsets(ziCandidates)

# --- Get names for all candidate predictors (for later metadata) ---
count_all_names <- names(countCandidates)
if (is.null(count_all_names) || any(count_all_names == "")) {
  count_all_names <- paste0("Pred", seq_along(countCandidates))
}
zi_all_names <- names(ziCandidates)
if (is.null(zi_all_names) || any(zi_all_names == "")) {
  zi_all_names <- paste0("Pred", seq_along(ziCandidates))
}

# --- Loop over all combinations of count and zero-inflation predictor subsets ---
# Instead of fitting a model, we simply build the merged data frame for each combination.
result_dfs <- list()
for (i in seq_along(count_subsets)) {
  for (j in seq_along(zi_subsets)) {
    # Start with the base response data.
    df <- df_base
    
    # --- Process count predictors ---
    count_combo <- count_subsets[[i]]
    if (length(count_combo$predictors) > 0) {
      for (pred_name in names(count_combo$predictors)) {
        add <- count_combo$predictors[[pred_name]]
        # Set row and column names using response names.
        rownames(add) <- plant_names
        colnames(add) <- hummingbird_names
        temp <- melt(add, varnames = c("plant", "hummingbird"), value.name = pred_name)
        df <- merge(df, temp, by = c("plant", "hummingbird"))
      }
    }
    
    # --- Process zero-inflation predictors ---
    zi_combo <- zi_subsets[[j]]
    if (length(zi_combo$predictors) > 0) {
      for (pred_name in names(zi_combo$predictors)) {
        add <- zi_combo$predictors[[pred_name]]
        rownames(add) <- plant_names
        colnames(add) <- hummingbird_names
        temp <- melt(add, varnames = c("plant", "hummingbird"), value.name = pred_name)
        df <- merge(df, temp, by = c("plant", "hummingbird"))
      }
    }
    
    # --- Optionally, add metadata columns to record the combination ---
    # Here we attach constant columns to every row in the data frame.
    if (length(count_combo$variants) > 0) {
      for (nm in names(count_combo$variants)) {
        df[[paste0("count_variant_", nm)]] <- count_combo$variants[[nm]]
      }
    }
    if (length(zi_combo$variants) > 0) {
      for (nm in names(zi_combo$variants)) {
        df[[paste0("zi_variant_", nm)]] <- zi_combo$variants[[nm]]
      }
    }
    
    # Save the data frame for this combination.
    result_dfs[[length(result_dfs) + 1]] <- df
  }
}

return(result_dfs)
}

#' @export
run_glmmTMB_models <- function(df_list, formula, ziformula) {
  # Load the required package
  library(glmmTMB)
  
  # Preallocate a results data.frame with one row per model.
  results <- data.frame(
    ID = character(length(df_list)),
    AIC = numeric(length(df_list)),
    delta_AIC = numeric(length(df_list)),
    stringsAsFactors = FALSE
  )
  
  # Loop through each data frame in the list
  for (i in seq_along(df_list)) {
    # If the list has names, use them; otherwise create an ID based on the index.
    df_id <- if (!is.null(names(df_list)) && names(df_list)[i] != "") {
      names(df_list)[i]
    } else {
      paste0("DF", i)
    }
    
    # Fit the model with family = nbinom1
    mod <- glmmTMB(
      formula = formula,
      data = df_list[[i]],
      ziformula = ziformula,
      family = nbinom1()
    )
    
    # Store the model's AIC and its ID
    results$ID[i] <- i
    results$AIC[i] <- AIC(mod)
  }
  
  # Calculate the delta AIC (difference from the best AIC)
  best_AIC <- min(results$AIC,na.rm=T)
  results$delta_AIC <- results$AIC - best_AIC
  
  # Sort results by ascending AIC
  results <- results[order(results$AIC), ]
  results= results[complete.cases(results),]
  
  return(results)
}

#' @export
prep_predictors=function(preds){
  n_mat=length(preds)
  if(n_mat>1){
    preds=Reduce(`*`, lapply(preds,function(x) x / sum(x)))
  } else preds = preds[[1]]
  
  pos=which(preds==0)
  preds[pos]=sort(unique(c(preds)))[2]/10
  
  preds=preds/sum(preds)
  return(preds)
}
#' @export
calc_AIC=function(pred_mat,net){
  n_mat=length(preds)
  if(n_mat>1){
    preds=Reduce(`*`, lapply(preds,function(x) x / sum(x)))
  } else preds = preds[[1]]
  
  pos=which(preds==0)
  preds[pos]=sort(unique(c(preds)))[2]/10
  
  preds=preds/sum(preds)
  
  AIC_X<- -2*dmultinom(c(net), prob=c(preds), log=T) + 2*(sum(dim(net))*n_mat)
  
  return(AIC_X) 
}
#' @export
all_combinations <- function(vec_list) {
  # Get the Cartesian product as a data.frame
  result <- do.call(expand.grid, vec_list)
  
  # If the list has names, assign them as column names for the result.
  if (!is.null(names(vec_list))) {
    colnames(result) <- names(vec_list)
  }
  
  return(result)
}



