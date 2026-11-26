
# Nopp-Mayr et. al
# From ingestion to egestion:
# Benchmarking FTIR spectroscopy and microscopy for identifying plant dietary components
# in feces from feeding trials 

# R Script for Band Aggregation and Random Forest Classification of plant parts in spectra of fecal samples

##### Packages #####
library(ranger)
library(caret)
library(tidyverse)
library(data.table)
library(rsample)


##### Functions for Feature Aggregation  #####
#Helper Functions
trapz <- function(x, y) {
  nx <- length(x)
  if (nx < 2) return(NA_real_)
  sum((x[-1] - x[-nx]) * (y[-1] + y[-nx])) / 2
}

lin_baseline <- function(x, y) {
  y0 <- y[1]; y1 <- y[length(y)]
  x0 <- x[1]; x1 <- x[length(x)]
  if (!is.finite(y0) || !is.finite(y1) || x1 == x0) return(y)
  m  <- (y1 - y0) / (x1 - x0)
  b  <- y0 - m * x0
  y - (m * x + b)
}

centroid <- function(x, y_pos) {
  S <- sum(y_pos)
  if (!is.finite(S) || S <= 0) return(NA_real_)
  sum(x * y_pos) / S
}

fwhm <- function(x, y) {
  ymax <- max(y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) return(NA_real_)
  half <- ymax / 2
  above <- which(y >= half)
  if (length(above) < 2) return(NA_real_)
  x[above[length(above)]] - x[above[1]]
}

band_features_one <- function(x_wn, y_abs, bands, calc_fwhm = FALSE,
                              min_points_area = 2,
                              min_points_height = 1,
                              adapt_expand = TRUE,
                              expand_step = 2,
                              max_expands = 3) {
  o <- order(x_wn)
  x <- x_wn[o]; y <- y_abs[o]
  out <- list()
  for (b in bands) {
    lo <- min(b$lower, b$upper); hi <- max(b$lower, b$upper)
    sel <- which(x >= lo & x <= hi)
    if (adapt_expand && length(sel) < min_points_height) {
      ex <- 0
      while (length(sel) < min_points_height && ex < max_expands) {
        ex <- ex + 1
        lo <- lo - expand_step
        hi <- hi + expand_step
        sel <- which(x >= lo & x <= hi)
      }
    }
    nm <- b$name
    
    A <- NA_real_; H <- NA_real_; Cen <- NA_real_; F <- NA_real_
    if (length(sel) >= min_points_height) {
      xx <- x[sel]; yy <- y[sel]
      yy_bc  <- lin_baseline(xx, yy)
      H <- max(yy_bc, na.rm = TRUE)
      
      if (length(sel) >= min_points_area) {
        yy_pos <- pmax(yy_bc, 0)
        A <- trapz(xx, yy_pos)
        S <- sum(yy_pos)
        if (is.finite(S) && S > 0) Cen <- centroid(xx, yy_pos)
        if (calc_fwhm) F <- fwhm(xx, yy_bc)
      } else {
        A <- 0
        Cen <- NA_real_
        if (calc_fwhm) F <- NA_real_
      }
    }
    vals <- c(A, H, Cen)
    names(vals) <- c(paste0(nm, "_area"),
                     paste0(nm, "_height"),
                     paste0(nm, "_centroid"))
    if (calc_fwhm) vals <- c(vals, setNames(F, paste0(nm, "_fwhm")))
    out[[nm]] <- vals
  }
  unlist(out)
}

# Band Aggregates
make_default_bands <- function() {
  list(
    # OH
    list(name = "OH_3490_3080", lower = 3080, upper = 3490), 
    
    # Lipids
    list(name = "TriGly_1750_1738",   lower = 1738, upper = 1750),
    list(name = "Lipids_1470_1445",    lower = 1445,  upper = 1470),
    
    # Proteine 
    list(name = "AmideI_1670_1625", lower = 1625, upper = 1670),
    list(name = "AmideII_1580_1520",lower = 1520, upper = 1580),
    
    #Aromatics
    list(name = "Phenolics_1606_1600", lower = 1600, upper = 1606),  
    list(name = "Lignin_1515_1500",   lower = 1500, upper = 1515),
    list(name = "Phenolics_1495_1475",  lower = 1475, upper = 1495), 
    
    # Carbohydrates
    list(name = "Cellu_1060_1050", lower = 1150, upper = 1060), 
    list(name = "Carb_1035_1020", lower = 1020, upper = 1035),
    list(name = "Cellu_1011_1000", lower = 1000, upper = 1011), 
    list(name = "Cellu_997_980", lower = 980, upper = 997),

    # Oxalat
    list(name = "Oxalate_1320_1310", lower = 1312, upper = 1320),
    list(name = "Fructose and Oxalate_780_773", lower = 773, upper = 779),

    # Tannins
    list(name = "Tannin_1724_1716",  lower = 1716, upper = 1724)
  )
}

#Main Functions
aggregate_bands_df <- function(df,
                               meta_cols = 1:6,
                               remove_region = c(2500, 1800),
                               bands = make_default_bands(),
                               calc_fwhm = FALSE) {
  stopifnot(all(meta_cols %in% seq_len(ncol(df))))
  spec_cols <- setdiff(seq_len(ncol(df)), meta_cols)
  if (length(spec_cols) < 3) stop("Not enough columns found")
  
  # Extract wavenumbers from column names
  coln <- colnames(df)[spec_cols]
  wn_vals <- suppressWarnings(as.numeric(coln))
  if (anyNA(wn_vals)) {
    wn_vals <- suppressWarnings(as.numeric(sub(".*?(\\d+\\.?\\d*).*", "\\1", coln)))
  }
  if (anyNA(wn_vals)) stop("Colnames must include wavenumbers ('3304' or 'wn_3304').")
  
  # Remove CO2 and crystal regions
  keep <- rep(TRUE, length(wn_vals))
  if (!is.null(remove_region) && length(remove_region) == 2) {
    lo <- min(remove_region); hi <- max(remove_region)
    keep <- !(wn_vals <= hi & wn_vals >= lo)
  }
  wn_keep <- wn_vals[keep]
  if (length(wn_keep) < 3) stop("Not enough wavenumbers left after artefact regions were removed")
  
  # Spectral matrix
  mat <- as.matrix(df[, spec_cols, drop = FALSE])[, keep, drop = FALSE]
  if (!is.numeric(mat)) stop("Spectral data is not numeric.")
  
  # Calculate band features 
  feat_list <- lapply(seq_len(nrow(mat)), function(i) {
    band_features_one(wn_keep, mat[i, ], bands, calc_fwhm)
  })
  lens <- vapply(feat_list, length, integer(1))
  if (length(unique(lens)) != 1) {
    stop("Band feature vectors have different lengths. Check band windows.")
  }
  feat_mat <- do.call(rbind, feat_list)
  feat_df  <- as.data.frame(feat_mat, check.names = TRUE)
  rownames(feat_df) <- rownames(df)
  
  cbind(df[, meta_cols, drop = FALSE], feat_df)
}


##### Application of Feature Aggregation ####

data <- fread("FK_raw_FTIR_data_10072026.csv", header=TRUE, dec=",", sep=";", stringsAsFactors = FALSE)
data <- fread("IBK_raw_FTIR_data_10072026.csv", header=TRUE, dec=",", sep=";", stringsAsFactors = FALSE)

#Create classification subset with needles, berries and leaves
data <- data[data$part != "leaves/catkins",]
data <- data[data$part != "berries/leaves",]
data <- data[data$part != "berries/shoots",]
data <- data[data$part != "leaves/shoots",]
data <- data[data$part != "leaves/shoots/berries",]
data <- data[data$part != "needles/berries",]
data <- data[data$part != "needles/buds",]

#Optional: Remove Rubus sp. and Quercus 
#data <- data[data$species != "Rubus sp.",]
#data <- data[data$species != "Quercus",]

table(data$part)
table(data$species)

#Prepare spectral matrix and metadata columns 
spectra <- as.matrix(data[, 7:ncol(data)])
y <- as.character(data[[4]])
ID <- as.character(data[[1]])

##### Standard normal variate normalization #############
snv <- function(x) (x - mean(x)) / sd(x)
spectra_snv <- t(apply(spectra, 1, snv))
spectra_snv <- as.data.frame(spectra_snv)
data_snv <- cbind( ID, y, spectra_snv)
colnames(data_snv)[1] <- "ID"
colnames(data_snv)[2] <- "part"

data_snv <- data_snv

##### Aggregate Features ####
result_agg <- aggregate_bands_df(data_snv,
                                 meta_cols    = 1:2,
                                 remove_region= c(2500, 1800),
                                 bands        = make_default_bands(),
                                 calc_fwhm    = FALSE)
result <- result_agg

meta_cols    <- 1:2
target_col   <- "part"   
id_col       <- "ID"         

cand_cols <- setdiff(seq_along(result), c(meta_cols, which(names(result) == target_col)))

# Select Band-Features (_area, _height, _centroid) 
feature_mask <- grepl("(?:_height)$", names(result)) |
  grepl("^RAW_", names(result))
feature_cols <- intersect(which(feature_mask), cand_cols)

is_constant <- vapply(result[ , feature_cols, drop=FALSE], function(v) {
  v_allna <- all(is.na(v))
  v_allna || (length(unique(v[!is.na(v)])) <= 1)
}, logical(1))
feature_cols <- feature_cols[!is_constant]

X_all <- as.matrix(result[ , feature_cols, drop = FALSE])

storage.mode(X_all) <- "numeric"

part <- as.factor(result[[target_col]])
id <- as.factor(result[[id_col]])


data_agg <- cbind.data.frame( id= id, part = part , as.data.frame(X_all), stringsAsFactors = FALSE)


##### Split Training/Testdata  #########
# Identical samples IDs are not separated in different data splits
set.seed(7)

id_part <- aggregate(part ~ ID, data = result, FUN = function(x) unique(x)[1])

id_groups <- split(id_part$ID, id_part$part)
train_ids <- unlist(lapply(id_groups, function(v) {
  v[sample.int(length(v), size = max(1, round(0.8 * length(v))))]
}), use.names = FALSE)

train_idx <- result$ID %in% train_ids
test_idx  <- !train_idx

X_train <- X_all[train_idx, , drop = FALSE]
X_test  <- X_all[test_idx, , drop = FALSE]
y_train <- part[train_idx]
y_test  <- part[test_idx]


train_set <- cbind.data.frame(y_train, X_train)
test_set  <- cbind.data.frame(y_test,  X_test)

train_set <- train_set %>% rename(part = y_train)
test_set <- test_set %>% rename(part = y_test)


##### Optional: SNV normalized data of raw wavenumbers without feature aggregation ####
meta_cols  <- 1:2
target_col <- "part"
id_col     <- "ID"

# Remove CO2 and crystal regions
wn_cols <- setdiff(names(data_snv), names(data_snv)[meta_cols])
wn_cols <- setdiff(wn_cols, c(target_col, id_col))
wn_numeric <- as.numeric(wn_cols)
keep <- is.na(wn_numeric) | wn_numeric < 1800 | wn_numeric > 2500
wn_cols <- wn_cols[keep]

is_constant <- vapply(data_snv[, wn_cols, drop = FALSE], function(v) {
  length(unique(v[!is.na(v)])) <= 1
}, logical(1))
wn_cols <- wn_cols[!is_constant]

X_all_r <- as.matrix(data_snv[, wn_cols, drop = FALSE])
storage.mode(X_all_r) <- "numeric"

part_r <- as.factor(data_snv[[target_col]])

#Split Training/Testdata
set.seed(7)

id_part_r   <- aggregate(
  as.formula(paste(target_col, "~", id_col)), 
  data = data_snv, 
  FUN  = function(x) unique(x)[1]
)

id_groups_r <- split(id_part_r[[id_col]], id_part_r[[target_col]])
train_ids_r <- unlist(lapply(id_groups_r, function(v) {
  v[sample.int(length(v), size = max(1, round(0.8 * length(v))))]
}), use.names = FALSE)

train_idx_r <- data_snv[[id_col]] %in% train_ids_r
test_idx_r  <- !train_idx_r

X_train_r <- X_all_r[train_idx_r, , drop = FALSE]
X_test_r  <- X_all_r[test_idx_r,  , drop = FALSE]
y_train_r <- part_r[train_idx_r]
y_test_r  <- part_r[test_idx_r]

train_set_r <- cbind.data.frame(y_train_r = y_train_r, as.data.frame(X_train_r))
test_set_r  <- cbind.data.frame(y_test_r  = y_test_r,  as.data.frame(X_test_r))

train_set_r <- train_set_r %>% rename(part = y_train_r)
test_set_r <- test_set_r %>% rename(part = y_test_r)


#### Hyperparameter Tuning ####

folds <- vfold_cv(train_set, v = 10)#

tune_grid <- expand.grid(
  num.trees = c(100, 300,500,1000),
  mtry = c(1,2,5,10,15),
  max.depth = c(10, 50, 100),
  min.node.size = c(1,2, 5)
)

cv_results <- tune_grid %>%
  purrr::pmap_dfr(function(num.trees, mtry, max.depth, min.node.size) {
    
    fold_accs <- map_dbl(folds$splits, function(split) {
      train_fold <- analysis(split)
      val_fold <- assessment(split)
      
      model <- ranger(
        dependent.variable.name = "y_train",
        data = train_fold,
        num.trees = num.trees,
        mtry = mtry,
        max.depth = max.depth,
        min.node.size = min.node.size
      )
      
      preds <- predict(model, data = val_fold)$predictions
      mean(preds == val_fold$y_train)
    })
    
    tibble(
      num.trees = num.trees,
      mtry = mtry,
      max.depth = max.depth,
      min.node.size = min.node.size,
      mean_accuracy = mean(fold_accs),
      sd_accuracy = sd(fold_accs)
    )
  })

cv_results %>% arrange(desc(mean_accuracy))

print(cv_results)


#### RandomForest Classification of plant parts in spectra of fecal samples ####
table(train_set$part)
table(test_set$part)

train_set$part <- as.factor(train_set$part)
test_set$part <- as.factor(test_set$part)

#RF model training
set.seed(7)

rf_model <- ranger(
  dependent.variable.name = "part",
  data = train_set,
  num.trees = 1000,
  mtry = 5,
  max.depth = 100,
  min.node.size = 1,
  importance = "permutation"
)

print(rf_model)

# RF model prediction
predictions <- predict(rf_model, data = test_set)

# Confusion matrix
conf_matrix <- confusionMatrix(
  predictions$predictions, 
  test_set$y_test
)

print(conf_matrix)  


#### Evaluate variable importance ##### 
n_runs <- 50
importance_list <- map_dfr(1:n_runs, function(i) {
  
  rf_vi <- ranger(
    dependent.variable.name = "part",
    data = train_set,
    num.trees = 1000,
    mtry = 5,
    max.depth = 50,
    min.node.size = 1,
    importance = "permutation"
  )
  
  tibble(
    run        = i,
    variable   = names(rf_vi$variable.importance),
    importance = unname(rf_vi$variable.importance)
  )
})

# Average runs
importance_summary <- importance_list %>%
  group_by(variable) %>%
  summarise(
    mean_importance = mean(importance),
    sd_importance   = sd(importance),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))

print(importance_summary, n = 15)

write.csv(importance_summary,
          "RF_part_variable_importance.csv",
          row.names = FALSE)


# Forward Feature selection by cv accuracy (OOB error)
sorted_features <- importance_summary$variable

cv_acc <- map_dbl(seq_along(sorted_features), function(k) {
  feat <- sorted_features[1:k]
  sub_train <- train_set[, c("part", feat)]
  
  model <- ranger(
    dependent.variable.name = "part",
    data = sub_train, num.trees = 1000, mtry = min(k, 15)
  )
  1 - model$prediction.error  
})

tibble(n_features = seq_along(sorted_features), accuracy = cv_acc) %>%
  ggplot(aes(n_features, accuracy)) + 
  geom_line() + geom_point() +
  labs(title = "OOB Accuracy vs. Number of Features")

