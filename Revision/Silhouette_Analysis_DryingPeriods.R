#______________________________________________________________________________#
# Systematische Unterschiede zwischen zwei Fütterungsperioden
# Minimalversion: Permutationstest auf Feature- UND Vollspektrum-Ebene (ohne Plots)
#______________________________________________________________________________#

rm(list = ls())
library(tidyverse)
library(cluster)

##### Einstellungen ############################################################
REMOVE_REGION <- c(1800, 2500)
ID_COL        <- "ID"
PART_COL      <- "part"
SPECIES_COL   <- "species"
FEED_SP_COL   <- "sp"
DATE_COL      <- "date"
BLOCK_COL     <- "feeding_block"
MAX_GAP_DAYS  <- 3
RUN_LABEL     <- "Trocknungstest"
SEED          <- 7

PERIOD_MAP <- list(
  P1 = c("autumn_22", "summer_23"),
  P2 = c("autumn_23", "summer_24", "autumn_24")
)
N_PERM_TEST        <- 999
MIN_IDS_PER_PERIOD <- 2

####### Feature Aggregation ####################################################
trapz <- function(x, y) { nx <- length(x); if (nx < 2) return(NA_real_)
sum((x[-1] - x[-nx]) * (y[-1] + y[-nx])) / 2 }

lin_baseline <- function(x, y) {
  y0 <- y[1]; y1 <- y[length(y)]; x0 <- x[1]; x1 <- x[length(x)]
  if (!is.finite(y0) || !is.finite(y1) || x1 == x0) return(y)
  m <- (y1 - y0) / (x1 - x0); y - (m * x + (y0 - m * x0))
}
centroid <- function(x, y_pos) { S <- sum(y_pos)
if (!is.finite(S) || S <= 0) return(NA_real_); sum(x * y_pos) / S }

band_features_one <- function(x_wn, y_abs, bands,
                              min_points_area = 2, min_points_height = 1,
                              expand_step = 2, max_expands = 3) {
  o <- order(x_wn); x <- x_wn[o]; y <- y_abs[o]; out <- list()
  for (b in bands) {
    lo <- min(b$lower, b$upper); hi <- max(b$lower, b$upper)
    sel <- which(x >= lo & x <= hi); ex <- 0
    while (length(sel) < min_points_height && ex < max_expands) {
      ex <- ex + 1; lo <- lo - expand_step; hi <- hi + expand_step
      sel <- which(x >= lo & x <= hi)
    }
    A <- NA_real_; H <- NA_real_; Cen <- NA_real_
    if (length(sel) >= min_points_height) {
      xx <- x[sel]; yy_bc <- lin_baseline(xx, y[sel]); H <- max(yy_bc, na.rm = TRUE)
      if (length(sel) >= min_points_area) {
        yy_pos <- pmax(yy_bc, 0); A <- trapz(xx, yy_pos); S <- sum(yy_pos)
        if (is.finite(S) && S > 0) Cen <- centroid(xx, yy_pos)
      } else A <- 0
    }
    vals <- c(A, H, Cen)
    names(vals) <- paste0(b$name, c("_area", "_height", "_centroid"))
    out[[b$name]] <- vals
  }
  unlist(unname(out))
}

make_default_bands <- function() list(
  list(name = "OH_3490_3080",            lower = 3080, upper = 3490),
  list(name = "TriGly_1750_1738",        lower = 1738, upper = 1750),
  list(name = "Lipids_1470_1445",        lower = 1445, upper = 1470),
  list(name = "AmideII_1580_1520",       lower = 1520, upper = 1580),
  list(name = "AmideI_1653_1645",        lower = 1645, upper = 1653),
  list(name = "Phenolics_1606_1600",     lower = 1600, upper = 1606),
  list(name = "Lignin_1515_1500",        lower = 1500, upper = 1515),
  list(name = "Cellu_1060_1050",         lower = 1050, upper = 1060),
  list(name = "Carb_1035_1020",          lower = 1020, upper = 1035),
  list(name = "Cellu_1011_1000",         lower = 1000, upper = 1011),
  list(name = "Cellu_997_980",           lower = 980,  upper = 997),
  list(name = "FructoseOx_780_773",      lower = 773,  upper = 780),
  list(name = "Oxalate_1320_1310",       lower = 1310, upper = 1320),
  list(name = "Tannin_1724_1716",        lower = 1716, upper = 1724),
  list(name = "Resin_Phenolics_1495_1475", lower = 1475, upper = 1495)
)

aggregate_bands_df <- function(df, meta_cols, bands = make_default_bands()) {
  spec_cols <- setdiff(seq_len(ncol(df)), meta_cols)
  wn_vals <- suppressWarnings(as.numeric(sub(".*?(\\d+\\.?\\d*).*", "\\1",
                                             colnames(df)[spec_cols])))
  mat <- as.matrix(df[, spec_cols, drop = FALSE])
  feat_mat <- do.call(rbind, lapply(seq_len(nrow(mat)),
                                    function(i) band_features_one(wn_vals, mat[i, ], bands)))
  cbind(df[, meta_cols, drop = FALSE], as.data.frame(feat_mat, check.names = TRUE))
}

##### Daten laden ##############################################################
setwd("C:/Users/schae/Desktop/FWF Grouse/Feeding Experiment FK IBK/Revision/data")
data <- read.csv("IBK_L_raw_FTIR_data.csv", header = TRUE, dec = ",",
                 sep = ";", stringsAsFactors = FALSE)

keep_parts   <- c("berries", "needles", "leaves")
drop_species <- character(0)
data$part    <- trimws(data$part)
data$species <- trimws(data$species)

##### Fütterungsblöcke ########################################################
parse_date <- function(x) {
  x <- trimws(as.character(x))
  for (f in c("%Y-%m-%d","%d.%m.%Y","%d.%m.%y","%d/%m/%Y","%Y/%m/%d","%m/%d/%Y")) {
    d <- as.Date(x, format = f); if (all(!is.na(d[x != ""]))) return(d)
  }
  stop("Datum nicht lesbar.")
}
data[[BLOCK_COL]]   <- trimws(as.character(data[[BLOCK_COL]]))
data[[FEED_SP_COL]] <- trimws(as.character(data[[FEED_SP_COL]]))
data$.date <- parse_date(data[[DATE_COL]])

blocks <- data %>%
  distinct(ID = .data[[ID_COL]], block = .data[[BLOCK_COL]], date = .date,
           feed_sp = .data[[FEED_SP_COL]], feed_part = .data[[PART_COL]]) %>%
  group_by(block) %>%
  summarise(saison = sub("_[0-9]+$", "", block[1]),
            feed_sp = paste(unique(feed_sp), collapse = " + "),
            feed_part = paste(unique(feed_part), collapse = " + "), .groups = "drop") %>%
  mutate(futter = ifelse(endsWith(feed_sp, paste0("_", feed_part)) |
                           endsWith(feed_sp, paste0(" ", feed_part)),
                         feed_sp, paste(feed_sp, feed_part)))

norm_sp <- function(x) sub("\\.$", "", tolower(gsub("\\s+", " ", trimws(x))))
is_drop <- norm_sp(data$species) %in% norm_sp(drop_species)
data <- data[data$part %in% keep_parts & !is_drop, ]
stopifnot(nrow(data) > 0)
cat("Verbleibende Spektren:", nrow(data), "\n")

##### Spektren: Artefaktregion entfernen, dann SNV ############################
spec_cols <- grep("^X?[0-9]+([.,][0-9]+)?$", names(data), value = TRUE)
spectra <- as.matrix(data[, spec_cols]); storage.mode(spectra) <- "numeric"
wn <- suppressWarnings(as.numeric(gsub(",", ".", sub("^[^0-9]*", "", colnames(spectra)))))
keep_wn <- !(wn >= min(REMOVE_REGION) & wn <= max(REMOVE_REGION))
spectra <- spectra[, keep_wn, drop = FALSE]; wn <- wn[keep_wn]

snv <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
spectra_snv <- t(apply(spectra, 1, snv))
colnames(spectra_snv) <- paste0("wn_", wn)

ID   <- as.character(data[[ID_COL]])
part <- factor(data[[PART_COL]])

##### Band-Feature-Matrix #####################################################
df_snv <- data.frame(ID = ID, part = part, spectra_snv, check.names = FALSE)
agg    <- aggregate_bands_df(df_snv, meta_cols = 1:2)
feat   <- agg[, grepl("_height$", names(agg)), drop = FALSE]
bad    <- vapply(feat, function(v) anyNA(v) || length(unique(v)) <= 1, logical(1))
X_all  <- as.matrix(feat[, !bad, drop = FALSE]); storage.mode(X_all) <- "numeric"
cat("Band-Features:", ncol(X_all), "\n")

##### Perioden zuordnen #######################################################
saison_map <- unlist(lapply(names(PERIOD_MAP),
                            function(p) setNames(rep(p, length(PERIOD_MAP[[p]])), PERIOD_MAP[[p]])))
blocks$periode <- unname(saison_map[blocks$saison])

id_info <- tibble(ID = ID, block = data[[BLOCK_COL]], part = as.character(part)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  left_join(select(blocks, block, futter, periode, saison), by = "block") %>%
  filter(!is.na(periode))

# Futter in beiden Perioden mit genügend Proben
uebersicht <- id_info %>% count(futter, periode) %>%
  pivot_wider(names_from = periode, values_from = n, values_fill = 0, names_prefix = "ids_")
futter_beide <- uebersicht %>%
  filter(ids_P1 >= MIN_IDS_PER_PERIOD, ids_P2 >= MIN_IDS_PER_PERIOD) %>% pull(futter)
cat("Futter in beiden Perioden:", length(futter_beide), "-",
    paste(futter_beide, collapse = ", "), "\n")
if (length(futter_beide) == 0) stop("Kein Futter in beiden Perioden.")

####### Generische Analyse-Funktion (ohne Plots) ##############################
sil_futter_g <- function(prof_mat, info, futter_vec) {
  map_dfr(futter_vec, function(fu) {
    sub <- info %>% filter(futter == fu)
    s <- silhouette(as.integer(factor(sub$periode)), dist(prof_mat[sub$ID, , drop = FALSE]))
    tibble(futter = fu, ID = sub$ID, periode = sub$periode, block = sub$block,
           s = s[, "sil_width"])
  })
}
period_means_g <- function(prof_mat, info, futter_vec) {
  map(futter_vec, function(fu) {
    sub <- info %>% filter(futter == fu)
    list(m1 = colMeans(prof_mat[sub$ID[sub$periode == "P1"], , drop = FALSE]),
         m2 = colMeans(prof_mat[sub$ID[sub$periode == "P2"], , drop = FALSE]))
  })
}
stat_mean_g  <- function(sil) mean(sil$s)
stat_dist_g  <- function(pm)  mean(vapply(pm, function(x) sqrt(sum((x$m2 - x$m1)^2)), numeric(1)))
stat_bands_g <- function(pm)  Reduce(`+`, lapply(pm, function(x) x$m2 - x$m1)) / length(pm)

run_period_analysis <- function(X_feat, ID_vec, info, futter_vec, label) {
  cat("\n### Periodentest:", label, "(", ncol(X_feat), "Features )\n")
  
  X_z <- scale(X_feat)
  X_z <- X_z[, apply(X_z, 2, function(v) all(is.finite(v))), drop = FALSE]
  cnt <- as.vector(table(ID_vec)[rownames(rowsum(X_z, ID_vec))])
  prof_mat <- (rowsum(X_z, ID_vec) / cnt)[info$ID, , drop = FALSE]
  
  sil_obs   <- sil_futter_g(prof_mat, info, futter_vec)
  obs       <- stat_mean_g(sil_obs)
  pm_obs    <- period_means_g(prof_mat, info, futter_vec)
  obs_dist  <- stat_dist_g(pm_obs)
  obs_bands <- stat_bands_g(pm_obs)
  
  set.seed(SEED)
  blk_per    <- info %>% filter(futter %in% futter_vec) %>% distinct(futter, block, periode)
  info_beide <- info %>% filter(futter %in% futter_vec)
  perm <- map(seq_len(N_PERM_TEST), function(i) {
    bp <- blk_per %>% group_by(futter) %>% mutate(periode = sample(periode)) %>% ungroup()
    info_p <- info_beide %>% select(-periode) %>% left_join(select(bp, block, periode), by = "block")
    pm <- period_means_g(prof_mat, info_p, futter_vec)
    list(sil = stat_mean_g(sil_futter_g(prof_mat, info_p, futter_vec)),
         dist = stat_dist_g(pm), bands = stat_bands_g(pm))
  })
  perm_sil   <- vapply(perm, `[[`, numeric(1), "sil")
  perm_dist  <- vapply(perm, `[[`, numeric(1), "dist")
  perm_bands <- do.call(rbind, lapply(perm, `[[`, "bands"))
  pval <- function(o, null) (1 + sum(null >= o)) / (length(null) + 1)
  
  test <- tibble(
    statistik = c("Silhouette nach Periode (innerhalb Futter)",
                  "Abstand der Periodenmittelwerte (innerhalb Futter)"),
    beobachtet = c(obs, obs_dist),
    null_mittel = c(mean(perm_sil), mean(perm_dist)),
    null_q95 = c(quantile(perm_sil, 0.95), quantile(perm_dist, 0.95)),
    p_permutation = c(pval(obs, perm_sil), pval(obs_dist, perm_dist)),
    n_futter = length(futter_vec), n_ids = nrow(sil_obs),
    n_bloecke = n_distinct(sil_obs$block), ebene = label)
  print(as.data.frame(test), digits = 3, row.names = FALSE)
  
  band_diff <- map_dfr(futter_vec, function(fu) {
    sub <- info %>% filter(futter == fu)
    m1 <- colMeans(prof_mat[sub$ID[sub$periode == "P1"], , drop = FALSE])
    m2 <- colMeans(prof_mat[sub$ID[sub$periode == "P2"], , drop = FALSE])
    tibble(futter = fu, feature = colnames(prof_mat), differenz_P2_minus_P1 = m2 - m1)
  })
  band_p <- vapply(colnames(prof_mat), function(f)
    (1 + sum(abs(perm_bands[, f]) >= abs(obs_bands[f]))) / (N_PERM_TEST + 1), numeric(1))
  band_sum <- band_diff %>% group_by(feature) %>%
    summarise(mittlere_differenz = mean(differenz_P2_minus_P1),
              anteil_P2_hoeher = mean(differenz_P2_minus_P1 > 0),
              n_futter = n(), .groups = "drop") %>%
    mutate(p_permutation = band_p[feature],
           p_bh = p.adjust(p_permutation, method = "BH")) %>%
    arrange(p_permutation, desc(abs(mittlere_differenz)))
  cat("\nFeatures: Differenz P2 - P1 (", label, ") - Top 15:\n", sep = "")
  print(as.data.frame(head(band_sum, 15)), digits = 3, row.names = FALSE)
  
  invisible(list(test = test, band_sum = band_sum))
}

##### Beide Ebenen rechnen ####################################################
res_features <- run_period_analysis(X_all,       ID, id_info, futter_beide, "Features")
res_vollspek <- run_period_analysis(spectra_snv, ID, id_info, futter_beide, "Vollspektrum")

##### Gegenüberstellung #######################################################
cat("\n### Gegenueberstellung beider Ebenen ###\n")
print(as.data.frame(bind_rows(res_features$test, res_vollspek$test)),
      digits = 3, row.names = FALSE)