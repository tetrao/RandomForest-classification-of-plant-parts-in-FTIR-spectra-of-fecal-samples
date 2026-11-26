# Revision - Nopp-Mayr et. al
# From ingestion to egestion:
# Benchmarking FTIR spectroscopy and microscopy for identifying plant dietary components
# in feces from feeding trials 

# R Script for Feature Aggregation and Random Forest Classification of plant parts in spectra of fecal samples 
# including sample-grouped cross-validation for hyperparameter optimzation and variable importance analysis
# including carryover analysis between consecutive feeding blocks (feature level and full-spectrum level)


rm(list=ls())

setwd("https://raw.githubusercontent.com/tetrao/RandomForest-classification-of-plant-parts-in-FTIR-spectra-of-fecal-samples/refs/heads/main/Revision") # Adapt path

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


band_features_one <- function(x_wn, y_abs, bands, 
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
    out[[nm]] <- vals
  }
  unlist(out)
}

# Band Aggregates
make_default_bands <- function() {
  list(
    list(name = "OH_3490_3080", lower = 3080, upper = 3490), 
    list(name = "TriGly_1750_1738",   lower = 1738, upper = 1750),
    list(name = "Lipids_1470_1445",    lower = 1445,  upper = 1470),
    list(name = "AmideII_1580_1520",lower = 1520, upper = 1580),
    list(name = "AmideI_1653_1645", lower = 1645, upper = 1653),
    list(name = "Phenolics_1606_1600", lower = 1600, upper = 1606),  
    list(name = "Lignin_1515_1500",   lower = 1500, upper = 1515),
    list(name = "Cellu_1060_1050", lower = 1050, upper = 1060), 
    list(name = "Carb_1035_1020", lower = 1020, upper = 1035),
    list(name = "Cellu_1011_1000", lower = 1000, upper = 1011), 
    list(name = "Cellu_997_980", lower = 980, upper = 997),
    list(name = "FructoseOx_780_773", lower = 773, upper = 780),
    list(name = "Oxalate_1320_1312", lower = 1312, upper = 1320),
    list(name = "Tannin_1724_1716",  lower = 1716, upper = 1724),
    list(name = "Resin_Phenolics_1495_1475",  lower = 1475, upper = 1495) 
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
  
  coln <- colnames(df)[spec_cols]
  wn_vals <- suppressWarnings(as.numeric(coln))
  if (anyNA(wn_vals)) {
    wn_vals <- suppressWarnings(as.numeric(sub(".*?(\\d+\\.?\\d*).*", "\\1", coln)))
  }
  if (anyNA(wn_vals)) stop("Colnames must include wavenumbers ('3304' or 'wn_3304').")
  
  keep <- rep(TRUE, length(wn_vals))
  if (!is.null(remove_region) && length(remove_region) == 2) {
    lo <- min(remove_region); hi <- max(remove_region)
    keep <- !(wn_vals <= hi & wn_vals >= lo)
  }
  wn_keep <- wn_vals[keep]
  if (length(wn_keep) < 3) stop("Not enough wavenumbers left after artefact regions were removed")
  
  mat <- as.matrix(df[, spec_cols, drop = FALSE])[, keep, drop = FALSE]
  if (!is.numeric(mat)) stop("Spectral data is not numeric.")
  
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

data <- fread("https://raw.githubusercontent.com/tetrao/RandomForest-classification-of-plant-parts-in-FTIR-spectra-of-fecal-samples/refs/heads/main/Revision/IBK_L_raw_FTIR_data.csv", header=TRUE, dec=",", sep=";", stringsAsFactors = FALSE)
data <- as.data.frame(data)

# Ungefilterte Kopie (inkl. Mischfutter) fuer den Carry-over-Test am Ende des Skripts
data_all <- data

#Create classification subset with needles, berries and leaves
data <- data[data$part != "leaves/catkins",]
data <- data[data$part != "berries/leaves",]
data <- data[data$part != "berries/shoots",]
data <- data[data$part != "leaves/shoots",]
data <- data[data$part != "leaves/shoots/berries",]
data <- data[data$part != "needles/berries",]

table(data$part)
table(data$species)

#Prepare spectral matrix and metadata columns 
spec_cols <- grep("^X?[0-9]+([.,][0-9]+)?$", names(data), value = TRUE)

wn_all <- as.numeric(gsub(",", ".", sub("^X", "", spec_cols)))
if (anyNA(wn_all)) stop("Wellenzahlen nicht lesbar: ", paste(head(spec_cols[is.na(wn_all)]), collapse = ", "))
spec_cols <- spec_cols[!(wn_all >= 1800 & wn_all <= 2500)]
cat("Spektralspalten nach Entfernen der Artefaktregion:", length(spec_cols), "\n")

spectra <- as.matrix(data[, spec_cols])
y <- as.character(data[["part"]])
ID <- as.character(data[["ID"]])

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



#### Set hyperparameters automatically ####
USE_HPO_BEST    <- FALSE    
TIE_TOL         <- 0.005   
MANUAL_MTRY     <- 1    # optimized 
MANUAL_MIN_NODE <- 1    # selected from HPO list (performs better than the value from the HPO  (50), leaves are rarely classified with a high min.node.size)
N_TREES_FINAL   <- 1000
MODEL_FILE      <- "RF_final_model.rds"


#### Hyperparameter Tuning ####
train_set_cv <- cbind.data.frame(ID = result$ID[train_idx], train_set)
folds <- group_vfold_cv(train_set_cv, group = ID, v = 10)

tune_grid <- expand.grid(
  mtry = c(1,2,5,10,15),
  min.node.size = c(1,2, 5, 10, 50,75)
)


set.seed(7)   
cv_results <- tune_grid %>%
  purrr::pmap_dfr(function(mtry, min.node.size) {
    
    fold_accs <- map_dbl(folds$splits, function(split) {
      train_fold <- analysis(split)   %>% select(-ID)   # ID nur fuer die Fold-Bildung, kein Feature
      val_fold   <- assessment(split) %>% select(-ID)
      
      model <- ranger(
        dependent.variable.name = "part",
        data = train_fold,
        mtry = mtry,
        min.node.size = min.node.size
      )
      
      preds <- predict(model, data = val_fold)$predictions
      mean(preds == val_fold$part)
    })
    
    tibble(mtry = mtry, min.node.size = min.node.size,
           mean_accuracy = mean(fold_accs), sd_accuracy = sd(fold_accs))
  }) %>%
  arrange(desc(mean_accuracy))
print(cv_results, n = Inf)

best_hp <- cv_results %>%
  filter(mean_accuracy >= max(mean_accuracy) - TIE_TOL) %>%
  arrange(desc(min.node.size), mtry, desc(mean_accuracy)) %>%
  slice(1)

if (USE_HPO_BEST) {
  final_mtry     <- best_hp$mtry
  final_min_node <- best_hp$min.node.size
} else {
  final_mtry     <- MANUAL_MTRY
  final_min_node <- MANUAL_MIN_NODE
}
cat("Finale Hyperparameter:", if (USE_HPO_BEST) "(aus HPO)" else "(manuell)",
    "| mtry =", final_mtry, "| min.node.size =", final_min_node,
    "| num.trees =", N_TREES_FINAL, "\n")



#### RandomForest classification of plant parts in spectra of fecal samples ####
table(train_set$part)
table(test_set$part)

train_set$part <- as.factor(train_set$part)
test_set$part <- as.factor(test_set$part)

#RF model training
set.seed(7)
rf_model <- ranger(
  dependent.variable.name = "part",
  data = train_set,
  num.trees = N_TREES_FINAL,
  mtry = final_mtry,
  min.node.size = final_min_node,
  importance = "permutation"
)
print(rf_model)

saveRDS(rf_model, MODEL_FILE)   

# RF model prediction
predictions <- predict(rf_model, data = test_set)

# Confusion matrix
conf_matrix <- confusionMatrix(
  predictions$predictions, 
  test_set$part
)

print(conf_matrix)  

saveRDS(list(model = rf_model, features = colnames(X_all), classes = levels(part),
             mtry = final_mtry, min.node.size = final_min_node, num.trees = N_TREES_FINAL,
             bands = make_default_bands(), train_ids = train_ids, hpo = cv_results),
        MODEL_FILE)


#### Evaluate variable importance ##### 

vi_num_trees     <- N_TREES_FINAL
vi_mtry          <- final_mtry
vi_min_node_size <- final_min_node

n_repeats <- 5  
n_perm    <- 10

train_set_cv_vi <- cbind.data.frame(ID = result$ID[train_idx], train_set)

vi_grouped <- map_dfr(seq_len(n_repeats), function(rep) {
  
  set.seed(1000 + rep)
  folds_rep <- group_vfold_cv(train_set_cv_vi, group = ID, v = 10)
  
  map_dfr(seq_along(folds_rep$splits), function(k) {
    
    split      <- folds_rep$splits[[k]]
    train_fold <- analysis(split)   %>% select(-ID)
    val_fold   <- assessment(split) %>% select(-ID)
    
    model <- ranger(
      dependent.variable.name = "part",
      data = train_fold,
      num.trees     = vi_num_trees,
      mtry          = vi_mtry,
      min.node.size = vi_min_node_size
    )
    
    base_acc <- mean(predict(model, val_fold)$predictions == val_fold$part)
    
    feats <- setdiff(names(val_fold), "part")
    
    map_dfr(feats, function(f) {
      
      drops <- map_dbl(seq_len(n_perm), function(p) {
        vf <- val_fold
        vf[[f]] <- sample(vf[[f]])                 # Feature permutieren
        perm_acc <- mean(predict(model, vf)$predictions == vf$part)
        base_acc - perm_acc                        # Accuracy-Abfall = Importance
      })
      
      tibble(
        repeat_id  = rep,
        fold       = k,
        variable   = f,
        importance = mean(drops)
      )
    })
  })
})


importance_summary <- vi_grouped %>%
  group_by(variable) %>%
  summarise(
    mean_importance = mean(importance),
    sd_importance   = sd(importance),
    n               = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_importance))

print(importance_summary, n = 50)

write.csv(importance_summary,
          "RF_part_variable_importance.csv",
          row.names = FALSE)



#### Carry-over auf Feature-Ebene ####

MAX_GAP_DAYS <- 1  
CF_DIR <- file.path(getwd(), "carry_over_features")
dir.create(CF_DIR, showWarnings = FALSE)

parse_date <- function(x) {
  x <- trimws(as.character(x))
  for (f in c("%Y-%m-%d", "%d.%m.%Y", "%d.%m.%y", "%d/%m/%Y", "%Y/%m/%d", "%m/%d/%Y")) {
    d <- as.Date(x, format = f)
    if (all(!is.na(d[x != ""]))) return(d)
  }
  stop("Datum nicht lesbar. Beispiele: ", paste(head(unique(x)), collapse = ", "))
}

data_all$.date <- parse_date(data_all$date)
data_all$feeding_block <- trimws(as.character(data_all$feeding_block))

blocks <- data_all %>%
  distinct(ID = as.character(ID), block = feeding_block, date = .date,
           feed_sp = trimws(sp), feed_part = trimws(part)) %>%
  group_by(block) %>%
  summarise(saison    = sub("_[0-9]+$", "", block[1]),
            start     = min(date), ende = max(date),
            n_tage    = as.integer(max(date) - min(date)) + 1L,
            n_ids     = n_distinct(ID),
            feed_sp   = paste(unique(feed_sp), collapse = " + "),
            feed_part = paste(unique(feed_part), collapse = " + "), .groups = "drop") %>%
  arrange(saison, start) %>%
  group_by(saison) %>%
  mutate(vorg_block     = lag(block),
         vorg_feed_sp   = lag(feed_sp),
         vorg_feed_part = lag(feed_part),
         luecke_tage    = as.integer(start - lag(ende))) %>%
  ungroup() %>%
  mutate(vorg_gueltig = !is.na(vorg_block) & !is.na(luecke_tage) &
           luecke_tage >= 1 & luecke_tage <= MAX_GAP_DAYS)
print(as.data.frame(blocks %>% select(block, start, ende, n_tage, n_ids, feed_sp, feed_part,
                                      vorg_block, luecke_tage, vorg_gueltig)), row.names = FALSE)
write.csv(blocks, file.path(CF_DIR, "fuetterungsbloecke.csv"), row.names = FALSE)


ID_all      <- as.character(data_all[["ID"]])
spectra_all <- as.matrix(data_all[, spec_cols])
snv_all     <- as.data.frame(t(apply(spectra_all, 1, snv)))
data_snv_all <- cbind(ID = ID_all, part = as.character(data_all[["part"]]), snv_all)
agg_all <- aggregate_bands_df(data_snv_all,
                              meta_cols     = 1:2,
                              remove_region = c(2500, 1800),
                              bands         = make_default_bands(),
                              calc_fwhm     = FALSE)
F_all <- scale(as.matrix(agg_all[, names(result)[feature_cols], drop = FALSE]))   # gleiche Features wie im Modell
F_all <- F_all[, apply(F_all, 2, function(v) all(is.finite(v))), drop = FALSE]


prof <- rowsum(F_all, ID_all)
prof <- prof / as.vector(table(ID_all)[rownames(prof)])
id_info <- data_all %>%
  distinct(ID = as.character(ID), block = feeding_block, date = .date) %>%
  mutate(tag = as.integer(date - blocks$start[match(block, blocks$block)]) + 1L)

eucl <- function(a, b) sqrt(sum((a - b)^2))
bl_eval <- blocks %>%
  filter(vorg_gueltig, paste(feed_sp, feed_part) != paste(vorg_feed_sp, vorg_feed_part))

cf_id <- map_dfr(seq_len(nrow(bl_eval)), function(k) {
  b  <- bl_eval$block[k]; vb <- bl_eval$vorg_block[k]
  prev_ids <- id_info %>% filter(block == vb) %>% filter(date == max(date)) %>% pull(ID)
  own_ids  <- id_info %>% filter(block == b) %>% pull(ID)
  if (length(prev_ids) == 0 || length(own_ids) < 2) return(NULL)
  map_dfr(own_ids, function(i) {
    a_i <- mean(vapply(setdiff(own_ids, i), function(j) eucl(prof[i, ], prof[j, ]), numeric(1)))
    b_i <- mean(vapply(prev_ids,            function(j) eucl(prof[i, ], prof[j, ]), numeric(1)))
    tibble(ID = i, block = b, tag = id_info$tag[id_info$ID == i][1],
           futter = paste(bl_eval$feed_sp[k], bl_eval$feed_part[k]),
           vorg_futter = paste(bl_eval$vorg_feed_sp[k], bl_eval$vorg_feed_part[k]),
           a_eigen = a_i, b_vorg = b_i, s = (b_i - a_i) / max(a_i, b_i))
  })
}) %>%
  mutate(tag_gruppe = factor(ifelse(tag == 1, "Tag 1", ifelse(tag == 2, "Tag 2", "Tag 3+")),
                             levels = c("Tag 1", "Tag 2", "Tag 3+")),
         tag1 = factor(ifelse(tag == 1, "Tag 1", "Tag 2+"), levels = c("Tag 2+", "Tag 1")))
write.csv(cf_id, file.path(CF_DIR, "feature_carry_over_je_ID.csv"), row.names = FALSE)

cf_tag <- cf_id %>% group_by(tag_gruppe) %>%
  summarise(n_ids = n(), n_bloecke = n_distinct(block),
            s_mittel = mean(s), s_median = median(s), s_sd = sd(s),
            a_eigen_mittel = mean(a_eigen), b_vorg_mittel = mean(b_vorg), .groups = "drop")
print(as.data.frame(cf_tag), digits = 3, row.names = FALSE)
write.csv(cf_tag, file.path(CF_DIR, "feature_carry_over_je_tag.csv"), row.names = FALSE)


cf_tests <- list()
paar_f <- cf_id %>% group_by(block, tag1) %>% summarise(s = mean(s), .groups = "drop") %>%
  pivot_wider(names_from = tag1, values_from = s) %>% filter(!is.na(`Tag 1`), !is.na(`Tag 2+`))
if (nrow(paar_f) >= 3) {
  w <- suppressWarnings(wilcox.test(paar_f$`Tag 1`, paar_f$`Tag 2+`, paired = TRUE, exact = FALSE))
  cf_tests$paar <- tibble(test = paste0("s Tag 1 vs. Tag 2+ gepaart je Block (Wilcoxon; n = ", nrow(paar_f), ")"),
                          differenz = mean(paar_f$`Tag 1` - paar_f$`Tag 2+`),
                          anteil_bloecke_tag1_niedriger = mean(paar_f$`Tag 1` < paar_f$`Tag 2+`),
                          chisq = NA_real_, p = w$p.value)
}
if (requireNamespace("lme4", quietly = TRUE) && n_distinct(cf_id$block) >= 3) {
  m <- tryCatch(lme4::lmer(s ~ tag1 + (1 | block), data = cf_id, REML = FALSE),
                error = function(e) { message("LMM: ", conditionMessage(e)); NULL })
  if (!is.null(m)) {
    a <- anova(update(m, . ~ . - tag1), m)
    cf_tests$lmm <- tibble(test = "s ~ Tag1 + (1|Block) LMM (Likelihood-Ratio-Test)",
                           differenz = lme4::fixef(m)[["tag1Tag 1"]],
                           anteil_bloecke_tag1_niedriger = NA_real_,
                           chisq = a$Chisq[2], p = a$`Pr(>Chisq)`[2])
  }
}
cf_tests_df <- bind_rows(cf_tests)
print(as.data.frame(cf_tests_df), digits = 3, row.names = FALSE)
write.csv(cf_tests_df, file.path(CF_DIR, "feature_carry_over_tests.csv"), row.names = FALSE)

cf_plot <- ggplot(cf_id, aes(x = tag_gruppe, y = s)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, fill = "grey90", width = 0.5) +
  geom_line(aes(group = block), colour = "grey60", alpha = 0.6) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(title = "Carry-over auf Feature-Ebene",
       subtitle = paste0("+1 = wie eigenes Futter, 0 = gleich aehnlich, -1 = wie Vorgaengerfutter\n",
                         n_distinct(cf_id$block), " Bloecke; Linien verbinden Proben desselben Blocks"),
       x = NULL, y = "Silhouette s = (b - a) / max(a, b)") +
  theme_minimal()
print(cf_plot)
ggsave(file.path(CF_DIR, "feature_carry_over.png"), cf_plot, width = 7, height = 5, dpi = 300, bg = "white")

confint(m, parm = "tag1Tag 1", method = "profile")  



#### Optional: SNV normalized data of raw wavenumbers without feature aggregation ####
meta_cols  <- 1:2
target_col <- "part"
id_col     <- "ID"


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



#### Carry-over auf Vollspektrum-Ebene ####

CS_DIR <- file.path(getwd(), "carry_over_spectrum")
dir.create(CS_DIR, showWarnings = FALSE)

spec_cols_all <- grep("^X?[0-9]+([.,][0-9]+)?$", names(data_all), value = TRUE)
wn_all_cs <- as.numeric(gsub(",", ".", sub("^X", "", spec_cols_all)))
if (anyNA(wn_all_cs)) stop("Wellenzahlen (data_all) nicht lesbar: ",
                           paste(head(spec_cols_all[is.na(wn_all_cs)]), collapse = ", "))
spec_cols_all <- spec_cols_all[!(wn_all_cs >= 1800 & wn_all_cs <= 2500)]
cat("Vollspektrum-Spalten nach Entfernen der Artefaktregion:", length(spec_cols_all), "\n")

spectra_full <- as.matrix(data_all[, spec_cols_all])
storage.mode(spectra_full) <- "numeric"
snv_full <- t(apply(spectra_full, 1, snv))

S_all <- scale(snv_full)
S_all <- S_all[, apply(S_all, 2, function(v) all(is.finite(v))), drop = FALSE]
cat("Vollspektrum-Spalten nach Entfernen nicht-endlicher Variablen:", ncol(S_all), "\n")

ID_all_cs <- as.character(data_all[["ID"]])

prof_s <- rowsum(S_all, ID_all_cs)
prof_s <- prof_s / as.vector(table(ID_all_cs)[rownames(prof_s)])


cs_id <- map_dfr(seq_len(nrow(bl_eval)), function(k) {
  b  <- bl_eval$block[k]; vb <- bl_eval$vorg_block[k]
  prev_ids <- id_info %>% filter(block == vb) %>% filter(date == max(date)) %>% pull(ID)
  own_ids  <- id_info %>% filter(block == b) %>% pull(ID)

  prev_ids <- intersect(prev_ids, rownames(prof_s))
  own_ids  <- intersect(own_ids,  rownames(prof_s))
  if (length(prev_ids) == 0 || length(own_ids) < 2) return(NULL)
  map_dfr(own_ids, function(i) {
    a_i <- mean(vapply(setdiff(own_ids, i), function(j) eucl(prof_s[i, ], prof_s[j, ]), numeric(1)))
    b_i <- mean(vapply(prev_ids,            function(j) eucl(prof_s[i, ], prof_s[j, ]), numeric(1)))
    tibble(ID = i, block = b, tag = id_info$tag[id_info$ID == i][1],
           futter = paste(bl_eval$feed_sp[k], bl_eval$feed_part[k]),
           vorg_futter = paste(bl_eval$vorg_feed_sp[k], bl_eval$vorg_feed_part[k]),
           a_eigen = a_i, b_vorg = b_i, s = (b_i - a_i) / max(a_i, b_i))
  })
}) %>%
  mutate(tag_gruppe = factor(ifelse(tag == 1, "Tag 1", ifelse(tag == 2, "Tag 2", "Tag 3+")),
                             levels = c("Tag 1", "Tag 2", "Tag 3+")),
         tag1 = factor(ifelse(tag == 1, "Tag 1", "Tag 2+"), levels = c("Tag 2+", "Tag 1")))

write.csv(cs_id, file.path(CS_DIR, "spektrum_carry_over_je_ID.csv"), row.names = FALSE)

cs_tag <- cs_id %>% group_by(tag_gruppe) %>%
  summarise(n_ids = n(), n_bloecke = n_distinct(block),
            s_mittel = mean(s), s_median = median(s), s_sd = sd(s),
            a_eigen_mittel = mean(a_eigen), b_vorg_mittel = mean(b_vorg), .groups = "drop")
print(as.data.frame(cs_tag), digits = 3, row.names = FALSE)
write.csv(cs_tag, file.path(CS_DIR, "spektrum_carry_over_je_tag.csv"), row.names = FALSE)


cs_tests <- list()
paar_s <- cs_id %>% group_by(block, tag1) %>% summarise(s = mean(s), .groups = "drop") %>%
  pivot_wider(names_from = tag1, values_from = s) %>% filter(!is.na(`Tag 1`), !is.na(`Tag 2+`))
if (nrow(paar_s) >= 3) {
  w <- suppressWarnings(wilcox.test(paar_s$`Tag 1`, paar_s$`Tag 2+`, paired = TRUE, exact = FALSE))
  cs_tests$paar <- tibble(test = paste0("s Tag 1 vs. Tag 2+ gepaart je Block (Wilcoxon; n = ", nrow(paar_s), ")"),
                          differenz = mean(paar_s$`Tag 1` - paar_s$`Tag 2+`),
                          anteil_bloecke_tag1_niedriger = mean(paar_s$`Tag 1` < paar_s$`Tag 2+`),
                          chisq = NA_real_, p = w$p.value)
}
if (requireNamespace("lme4", quietly = TRUE) && n_distinct(cs_id$block) >= 3) {
  m <- tryCatch(lme4::lmer(s ~ tag1 + (1 | block), data = cs_id, REML = FALSE),
                error = function(e) { message("LMM (Vollspektrum): ", conditionMessage(e)); NULL })
  if (!is.null(m)) {
    a <- anova(update(m, . ~ . - tag1), m)
    cs_tests$lmm <- tibble(test = "s ~ Tag1 + (1|Block) LMM (Likelihood-Ratio-Test)",
                           differenz = lme4::fixef(m)[["tag1Tag 1"]],
                           anteil_bloecke_tag1_niedriger = NA_real_,
                           chisq = a$Chisq[2], p = a$`Pr(>Chisq)`[2])
  }
}
cs_tests_df <- bind_rows(cs_tests)
print(as.data.frame(cs_tests_df), digits = 3, row.names = FALSE)
write.csv(cs_tests_df, file.path(CS_DIR, "spektrum_carry_over_tests.csv"), row.names = FALSE)

cs_plot <- ggplot(cs_id, aes(x = tag_gruppe, y = s)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, fill = "grey90", width = 0.5) +
  geom_line(aes(group = block), colour = "grey60", alpha = 0.6) +
  geom_point(alpha = 0.7, size = 1.8) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(title = "Carry-over auf Vollspektrum-Ebene",
       subtitle = paste0("+1 = wie eigenes Futter, 0 = gleich aehnlich, -1 = wie Vorgaengerfutter\n",
                         n_distinct(cs_id$block), " Bloecke; SNV, Artefaktregion entfernt; ",
                         ncol(S_all), " Wellenzahlen"),
       x = NULL, y = "Silhouette s = (b - a) / max(a, b)") +
  theme_minimal()
print(cs_plot)
ggsave(file.path(CS_DIR, "spektrum_carry_over.png"), cs_plot, width = 7, height = 5, dpi = 300, bg = "white")

confint(m, parm = "tag1Tag 1", method = "profile")  
