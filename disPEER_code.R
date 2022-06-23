library(MASS)
library(mdpeer)
library(stats)
library(igraph)
library(data.table)
library(readxl)
library(readr)
library(dplyr)
library(ggplot2)
library(Matrix)
library(gridExtra)
library(rms)
library(lme4)
library(textshape)
library(nloptr)
library(miscFuncs)
library(sjstats)
library(boot)
library(corrplot)
library(grid)
library(stringr)

### INFO ###
### all needed data files are stored in "~/Desktop/Scholar_Indiana/my_work/codes"/data_files" ###

setwd("~/Desktop/Scholar_Indiana/my_work/codes")

## centroids of cortical regions in 3D, all available regions ###
centroids_data = as.data.table(read_xls("data_files/Table_HCP_20171011.xls",
                                        col_names = TRUE,
                                        range = cell_cols("A:G"),
                                        sheet = 1))[c(1:33, 50:82), ]
n_regions = nrow(centroids_data) # number of all regions
centroids_data_left = centroids_data[1:(n_regions / 2), ] # centroids left hemi
centroids_data_right = centroids_data[(n_regions / 2 + 1):n_regions, ] # centroids right hemi

### geodesic matrices unordered ###
D_geodesic_left_true_unord = as.data.frame(fread("data_files/regeodesicdistancehcpdata/aparc_dist_matrix_lh.csv"))[-34, -34]
colnames(D_geodesic_left_true_unord) = rownames(D_geodesic_left_true_unord) = centroids_data_left$`aparc_aseg_region names`
D_geodesic_right_true_unord = as.data.frame(fread("data_files/regeodesicdistancehcpdata/aparc_dist_matrix_rh.csv"))[-34, -34]
colnames(D_geodesic_right_true_unord) = rownames(D_geodesic_right_true_unord) = centroids_data_right$`aparc_aseg_region names`

### ordering both hemi data according connectivity blocks ###
ord_blocks_l = order(centroids_data_left$`aparc_aseg_RSN 7`)
centroids_data_left = centroids_data_left[ord_blocks_l, ]
ord_blocks_r = order(centroids_data_right$`aparc_aseg_RSN 7`)
centroids_data_right = centroids_data_right[ord_blocks_r, ]
centroids = c("centroid X", "centroid Y", "centroid Z")

### removing the "islands" causing high variances in b covariance matrix ###
'%notin%' = Negate('%in%')
centroids_data_left = centroids_data_left[`aparc_aseg_RSN 7` %notin% c(3, 4, 6, 8)]
# centroids_data_right = centroids_data_right[`aparc_aseg_RSN 7` %notin% c(2, 3, 8)]
centroids_data_right = centroids_data_right[`aparc_aseg_RSN 7` %notin% c(3, 4, 8)]

### ordering geodesic matrices ###
D_geodesic_left_true = D_geodesic_left_true_unord[centroids_data_left$`aparc_aseg_region names`, 
                                                  centroids_data_left$`aparc_aseg_region names`]
D_geodesic_right_true = D_geodesic_right_true_unord[centroids_data_right$`aparc_aseg_region names`, 
                                                    centroids_data_right$`aparc_aseg_region names`]
### number of regions in each hemi ###
p_l = nrow(centroids_data_left)
p_r = nrow(centroids_data_right)

### normalized distances matrices for both hemispheres ###
Dnorm_r_true_geo = calc_norm_matrix_d(as.matrix(D_geodesic_right_true))
Dnorm_l_true_geo = calc_norm_matrix_d(as.matrix(D_geodesic_left_true))

### unrelated subjects Ids ###
unrel_ids = as.data.table(read_xlsx("data_files/IDs_100_1113_unrelated_all_tasks.xlsx", col_names = "Id"))

### demographics data ###
dmgrph_data_rel = as.data.table(read_csv("data_files/unrestricted_asteiner_2_2_2021_14_19_44.csv", col_names = TRUE))
dmgrph_data = subset(dmgrph_data_rel, Subject %in% unrel_ids$Id)

### restricted demographics data ###
dmgrph_data_rel_restr = as.data.table(read_csv("data_files/RESTRICTED_asteiner_10_15_2021_11_20_30.csv", col_names = TRUE))
dmgrph_data_restr = subset(dmgrph_data_rel_restr, Subject %in% unrel_ids$Id)

### FreeSurfer data ###
FS_data_rel = as.data.table(read_csv("data_files/unrestricted_hcp_freesurfer.csv", col_names = TRUE))
FS_data = subset(FS_data_rel, Subject %in% unrel_ids$Id)

### euclidean distance matrices ###
D_euclidean_left_true = calc_eucl_dist(centroids_data_left[, ..centroids])
colnames(D_euclidean_left_true) = centroids_data_left$`aparc_aseg_region names`
rownames(D_euclidean_left_true) = centroids_data_left$`aparc_aseg_region names`
D_euclidean_right_true = calc_eucl_dist(centroids_data_right[, ..centroids])
colnames(D_euclidean_right_true) = centroids_data_right$`aparc_aseg_region names`
rownames(D_euclidean_right_true) = centroids_data_right$`aparc_aseg_region names`
Dnorm_r_true_euc = calc_norm_matrix_d(D_euclidean_right_true)
Dnorm_l_true_euc = calc_norm_matrix_d(D_euclidean_left_true)


# ##########################################
# centroids_data_left = centroids_data[1:(n_regions / 2), ] # centroids left hemi
# centroids_data_right = centroids_data[(n_regions / 2 + 1):n_regions, ] # centroids right hemi
# 
# ### ordering both hemi data according connectivity blocks ###
# ord_blocks_l = order(centroids_data_left$`aparc_aseg_RSN 7`)
# centroids_data_left = centroids_data_left[ord_blocks_l, ]
# ord_blocks_r = order(centroids_data_right$`aparc_aseg_RSN 7`)
# centroids_data_right = centroids_data_right[ord_blocks_r, ]
# centroids = c("centroid X", "centroid Y", "centroid Z")
# 
# ### number of regions in each hemi ###
# p_l = nrow(centroids_data_left)
# p_r = nrow(centroids_data_right)
# 
# excluded_left = centroids_data_left$`aparc_aseg_region names`[centroids_data_left$`aparc_aseg_RSN 7` %in% c(3, 4, 6, 8)]
# excluded_right = centroids_data_right$`aparc_aseg_region names`[centroids_data_right$`aparc_aseg_RSN 7` %in% c(2, 3, 8)]
# 
# excluded_left = centroids_data_left[centroids_data_left$`aparc_aseg_RSN 7` %in% c(3, 4, 6, 8), c(1, 3)]
# excluded_right = centroids_data_right[centroids_data_right$`aparc_aseg_RSN 7` %in% c(2, 3, 8), c(1, 3)]
# 
# 
# ### connectivity matrices ###
# Ac_l_true = construct_adj_matrix_c(p_l, 
#                                    centroids_data_left, 
#                                    "aparc_aseg_RSN 7")
# Ac_r_true = construct_adj_matrix_c(p_r, 
#                                    centroids_data_right, 
#                                    "aparc_aseg_RSN 7")
# grid.arrange(plot_matrix(Ac_l_true),
#              plot_matrix(Ac_r_true))
# 
# ### removing the "islands" causing high variances in b covariance matrix ###
# '%notin%' = Negate('%in%')
# centroids_data_left = centroids_data_left[`aparc_aseg_RSN 7` %notin% c(3, 4, 6, 8)]
# centroids_data_right = centroids_data_right[`aparc_aseg_RSN 7` %notin% c(2, 3, 8)]
# 
# ### ordering both hemi data according connectivity blocks ###
# ord_blocks_l = order(centroids_data_left$`aparc_aseg_RSN 7`)
# centroids_data_left = centroids_data_left[ord_blocks_l, ]
# ord_blocks_r = order(centroids_data_right$`aparc_aseg_RSN 7`)
# centroids_data_right = centroids_data_right[ord_blocks_r, ]
# centroids = c("centroid X", "centroid Y", "centroid Z")
# 
# ### number of regions in each hemi ###
# p_l = nrow(centroids_data_left)
# p_r = nrow(centroids_data_right)
# 
# Ac_l_true = construct_adj_matrix_c(p_l, 
#                                    centroids_data_left, 
#                                    "aparc_aseg_RSN 7")
# Ac_r_true = construct_adj_matrix_c(p_r, 
#                                    centroids_data_right, 
#                                    "aparc_aseg_RSN 7")
# grid.arrange(plot_matrix(Ac_l_true),
#              plot_matrix(Ac_r_true))
# 
# ###############################

### connectivity matrices ###
Ac_l_true = construct_adj_matrix_c(p_l,
                                   centroids_data_left,
                                   "aparc_aseg_RSN 7")
Ac_r_true = construct_adj_matrix_c(p_r,
                                   centroids_data_right,
                                   "aparc_aseg_RSN 7")


set.seed(2)
Ac_l_obs = calc_Aobs(Atrue = Ac_l_true, diss = 0.2)
set.seed(4)
Ac_r_obs = calc_Aobs(Atrue = Ac_r_true, diss = 0.2)
colnames(Ac_l_true) = rownames(Ac_l_true) = colnames(D_euclidean_left_true) 
colnames(Ac_l_obs) = rownames(Ac_l_obs) = colnames(D_euclidean_left_true)
colnames(Ac_r_true) = rownames(Ac_r_true) = colnames(D_euclidean_right_true) 
colnames(Ac_r_obs) = rownames(Ac_r_obs) = colnames(D_euclidean_right_true)

### true and observed normalized connectivity Laplacians ###
Lc_l_true = calc_laplacian_c(Ac_l_true)
Qc_l_true = calc_norm_laplacian_c(Lc_l_true)
Lc_r_true = calc_laplacian_c(Ac_r_true)
Qc_r_true = calc_norm_laplacian_c(Lc_r_true)
Lc_l_obs = calc_laplacian_c(Ac_l_obs)
Qc_l_obs = calc_norm_laplacian_c(Lc_l_obs)
Lc_r_obs = calc_laplacian_c(Ac_r_obs)
Qc_r_obs = calc_norm_laplacian_c(Lc_r_obs)
colnames(Qc_l_true) = rownames(Qc_l_true) = colnames(D_euclidean_left_true)
colnames(Qc_l_obs) = rownames(Qc_l_obs) = colnames(D_euclidean_left_true)
colnames(Qc_r_true) = rownames(Qc_r_true) = colnames(D_euclidean_right_true)
colnames(Qc_r_obs) = rownames(Qc_r_obs) = colnames(D_euclidean_right_true)

### proximity matrices ###
### EUCLIDEAN ###
Ad_r_true_euc = calc_adj_matrix_d(Dnorm_r_true_euc, h = 5)
Ad_l_true_euc = calc_adj_matrix_d(Dnorm_l_true_euc, h = 5)
range_r_true = summary(as.numeric(Ad_r_true_euc[lower.tri(Ad_r_true_euc, diag = FALSE)]))
as.numeric(log10(range_r_true[6] / range_r_true[1]))
range_l_true = summary(as.numeric(Ad_l_true_euc[lower.tri(Ad_l_true_euc, diag = FALSE)]))
as.numeric(log10(range_l_true[6] / range_l_true[1]))
### GEODESIC ###
Ad_r_true_geo = calc_adj_matrix_d(Dnorm_r_true_geo, h = 100)
Ad_l_true_geo = calc_adj_matrix_d(Dnorm_l_true_geo, h = 100)
range_r_true = summary(as.numeric(Ad_r_true_geo[lower.tri(Ad_r_true_geo, diag = FALSE)]))
as.numeric(log10(range_r_true[6] / range_r_true[1]))
range_l_true = summary(as.numeric(Ad_l_true_geo[lower.tri(Ad_l_true_geo, diag = FALSE)]))
as.numeric(log10(range_l_true[6] / range_l_true[1]))

### true normalized proximity laplacians ###
Qd_r_true_euc = calc_norm_laplacian_d(Ad_r_true_euc)
Qd_l_true_euc = calc_norm_laplacian_d(Ad_l_true_euc)
Qd_r_true_geo = calc_norm_laplacian_d(Ad_r_true_geo)
Qd_l_true_geo = calc_norm_laplacian_d(Ad_l_true_geo)

###### FOR THE PLOTS ######
# Ad_r_true_euc = calc_adj_matrix_d(Dnorm_r_true_euc, h = 5)
# Ad_r_true_geo = calc_adj_matrix_d(Dnorm_r_true_geo, h = 100)
# Ad_l_true_euc = calc_adj_matrix_d(Dnorm_l_true_euc, h = 5)
# Ad_l_true_geo = calc_adj_matrix_d(Dnorm_l_true_geo, h = 100)
# 
# plot_matrix(Ad_r_true_euc,
#             which_mat = " ",
#             what_mat = " ",
#             which_hemi = " ")
# plot_matrix(Ad_r_true_geo,
#             which_mat = " ",
#             what_mat = " ",
#             which_hemi = " ")
# grid.arrange(h1_euc, h1_geo, ncol = 2)
# grid.arrange(h5_euc, h5_geo, ncol = 2)
# grid.arrange(h25_euc, h25_geo, ncol = 2)
# grid.arrange(h100_euc, h100_geo, ncol = 2)

# norm(Qd_r_true_euc)
# norm(Qd_l_true_euc)
# norm(Qd_r_true_geo)
# norm(Qd_l_true_geo)
# norm(Qc_l_true)
# norm(Qc_l_obs)
# norm(Qc_r_true)
# norm(Qc_r_obs)

### brain data preparation ###
### choosing thicknesses of brain regions from FS_data: for Z matrix ###
L_regs_ord = c("FS_L_Cuneus_Thck",
               "FS_L_Fusiform_Thck",
               "FS_L_Lateraloccipital_Thck",
               "FS_L_Lingual_Thck",
               "FS_L_Pericalcarine_Thck",
               "FS_L_Paracentral_Thck",
               "FS_L_Postcentral_Thck",
               "FS_L_Precentral_Thck",
               "FS_L_Superiortemporal_Thck",
               "FS_L_Transversetemporal_Thck",
               "FS_L_Entorhinal_Thck",
               "FS_L_Inferiortemporal_Thck",
               "FS_L_Lateralorbitofrontal_Thck",
               "FS_L_Medialorbitofrontal_Thck",
               "FS_L_Frontalpole_Thck",
               "FS_L_Temporalpole_Thck",
               "FS_L_Bankssts_Thck",
               "FS_L_Caudalmiddlefrontal_Thck",
               "FS_L_Inferiorparietal_Thck",
               "FS_L_Isthmuscingulate_Thck",
               "FS_L_Middletemporal_Thck",
               "FS_L_Parsorbitalis_Thck",
               "FS_L_Parstriangularis_Thck" ,
               "FS_L_Posteriorcingulate_Thck",
               "FS_L_Precuneus_Thck",
               "FS_L_Rostralanteriorcingulate_Thck",
               "FS_L_Superiorfrontal_Thck")

# R_regs_ord = c("FS_R_Cuneus_Thck",
#                "FS_R_Fusiform_Thck",
#                "FS_R_Lateraloccipital_Thck",
#                "FS_R_Lingual_Thck",
#                "FS_R_Pericalcarine_Thck",
#                "FS_R_Caudalanteriorcingulate_Thck",
#                "FS_R_Supramarginal_Thck" ,
#                "FS_R_Entorhinal_Thck",
#                "FS_R_Inferiortemporal_Thck",
#                "FS_R_Lateralorbitofrontal_Thck",
#                "FS_R_Medialorbitofrontal_Thck",
#                "FS_R_Frontalpole_Thck",
#                "FS_R_Temporalpole_Thck",
#                "FS_R_Caudalmiddlefrontal_Thck",
#                "FS_R_Parsopercularis_Thck",
#                "FS_R_Parstriangularis_Thck",
#                "FS_R_Rostralmiddlefrontal_Thck",
#                "FS_R_Bankssts_Thck",
#                "FS_R_Inferiorparietal_Thck",
#                "FS_R_Isthmuscingulate_Thck",
#                "FS_R_Middletemporal_Thck",
#                "FS_R_Parsorbitalis_Thck",
#                "FS_R_Precuneus_Thck",
#                "FS_R_Rostralanteriorcingulate_Thck",
#                "FS_R_Superiorfrontal_Thck")

R_regs_ord = c("FS_R_Cuneus_Thck",
               "FS_R_Fusiform_Thck",
               "FS_R_Lateraloccipital_Thck",
               "FS_R_Lingual_Thck",
               "FS_R_Pericalcarine_Thck",
               "FS_R_Paracentral_Thck",
               "FS_R_Postcentral_Thck",
               "FS_R_Posteriorcingulate_Thck",
               "FS_R_Precentral_Thck",
               "FS_R_Superiortemporal_Thck",
               "FS_R_Transversetemporal_Thck",
               "FS_R_Entorhinal_Thck",
               "FS_R_Inferiortemporal_Thck",
               "FS_R_Lateralorbitofrontal_Thck",
               "FS_R_Medialorbitofrontal_Thck",
               "FS_R_Frontalpole_Thck",
               "FS_R_Temporalpole_Thck",
               "FS_R_Caudalmiddlefrontal_Thck",
               "FS_R_Parsopercularis_Thck",
               "FS_R_Parstriangularis_Thck",
               "FS_R_Rostralmiddlefrontal_Thck",
               "FS_R_Bankssts_Thck",
               "FS_R_Inferiorparietal_Thck",
               "FS_R_Isthmuscingulate_Thck",
               "FS_R_Middletemporal_Thck",
               "FS_R_Parsorbitalis_Thck",
               "FS_R_Precuneus_Thck",
               "FS_R_Rostralanteriorcingulate_Thck",
               "FS_R_Superiorfrontal_Thck")

X = data.table(
  gender = 1 * (dmgrph_data$Gender == "F"),
  # edu = subset(dmgrph_data_rel_restr, Subject %in% dmgrph_data$Subject)$SSAGA_Educ,
  age = subset(dmgrph_data_rel_restr, Subject %in% dmgrph_data$Subject)$Age)

Z_l = as.matrix(FS_data[, ..L_regs_ord])
Z_r = as.matrix(FS_data[, ..R_regs_ord])

ThckArea_left = FS_data[, c(327:359, 395:427)] ### without insula
ThckArea_right = FS_data[, c(361:393, 429:461)] ### without insula

### Cortical regions Volumes ###
Vol_left = ThckArea_left[, 1:33] * ThckArea_left[, 34:66]
Vol_right = ThckArea_right[, 1:33] * ThckArea_right[, 34:66]
colnames(Vol_left) = gsub("Thck", "Vol", 
                          colnames(ThckArea_left)[1:33])
colnames(Vol_right) = gsub("Thck", "Vol", 
                           colnames(ThckArea_right)[1:33])
ThckArea_left = cbind(ThckArea_left, Vol_left)
ThckArea_right = cbind(ThckArea_right, Vol_right)

corrMat_left = cor(ThckArea_left[, 1:33], ThckArea_left[, 34:66])
corrMat_right = cor(ThckArea_right[, 1:33], ThckArea_right[, 34:66])
corrTA_left = diag(corrMat_left)
corrTA_right = diag(corrMat_right)
par(mfrow = c(1, 2))

CorHistL = ggplot() +
  geom_histogram(aes(x = corrTA_left),
                 alpha = 0.2,
                 bins = 50,
                 fill = "royalblue",
                 col = "black") +
  ggtitle("Histogram of corresponding thickness and area correlations (left hem.)") +
  xlab("Correlation") +
  ylab("Count") +
  ylim(0, 8) +
  xlim(-1, 1) +
  theme_bw()

CorHistR = ggplot() +
  geom_histogram(aes(x = corrTA_right),
                 alpha = 0.2,
                 bins = 50,
                 fill = "royalblue",
                 col = "black") +
  ggtitle("Histogram of corresponding thickness and area correlations (right hem.)") +
  xlab("Correlation") +
  ylab("Count") +
  ylim(0, 8) +
  xlim(-1, 1) +
  theme_bw()

grid.arrange(CorHistL, CorHistR)

# noticed one high negative corr within left hem.
regions_names_left = colnames(ThckArea_left[, 1:33])
regions_names_right = colnames(ThckArea_right[, 1:33])
high_corr_number_left = which.max(abs(corrTA_left))
high_corr_number_right = which.max(abs(corrTA_right))
regions_names_left[high_corr_number_left] # FS_L_Isthmuscingulate_Thck
regions_names_right[high_corr_number_right] # FS_R_Caudalanteriorcingulate_Thck

### Thickness - Area high correlated regions ###
par(mfrow = c(1, 3))
plot(ThckArea_left$FS_L_Isthmuscingulate_Thck, 
     ThckArea_left$FS_L_Isthmuscingulate_Area,
     main = "Scatterplot of area and thickness measurements")

plot(ThckArea_right$FS_R_Caudalanteriorcingulate_Thck, 
     ThckArea_right$FS_R_Caudalanteriorcingulate_Area,
     main = "Scatterplot of area and thickness measurements")

plot(ThckArea_right$FS_R_Posteriorcingulate_Thck, 
     ThckArea_right$FS_R_Posteriorcingulate_Area,
     main = "Scatterplot of area and thickness measurements")

# loop for all the regions
# left hemisphere
Thck_left = ThckArea_left[, 1:33]
Area_left = ThckArea_left[, 34:66]
Volume_left = ThckArea_left[, 67:99]

# left hemisphere
Thck_right = ThckArea_right[, 1:33]
Area_right = ThckArea_right[, 34:66]
Volume_right = ThckArea_right[, 67:99]

### LEFT
### FS_L_Isthmuscingulate_Thck
reg_ThckArea_isth = data.table(thck = ThckArea_left$FS_L_Isthmuscingulate_Thck,
                               area = ThckArea_left$FS_L_Isthmuscingulate_Area)
reg_ThckArea_isth[, outlier := as.factor(ifelse((thck < 1.75 & area > 2000), 1, 0))]
reg_ThckArea_isth[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
outlier_subjects_1_l = which(reg_ThckArea_isth$outlier == 1)

# FS_L_Bankssts_Thck
reg_ThckArea_bank = data.table(thck = ThckArea_left$FS_L_Bankssts_Thck,
                               area = ThckArea_left$FS_L_Bankssts_Area,
                               vol = ThckArea_left$FS_L_Bankssts_Thck)
reg_ThckArea_bank[, outlier := as.factor(ifelse(thck < 2, 1, 0))]
reg_ThckArea_bank[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
outlier_subjects_2_l = which(reg_ThckArea_bank$outlier == 1)

### RIGHT
### FS_R_Caudalanteriorcingulate_Thck
reg_ThckArea_caudal = data.table(thck = ThckArea_right$FS_R_Caudalanteriorcingulate_Thck,
                                 area = ThckArea_right$FS_R_Caudalanteriorcingulate_Area)
reg_ThckArea_caudal[, outlier := as.factor(ifelse((thck < 2.25 & area > 1250), 1, 0))]
reg_ThckArea_caudal[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
outlier_subjects_1_r = which(reg_ThckArea_caudal$outlier == 1)

# FS_R_Posteriorcingulate_Thck
reg_ThckArea_postcing = data.table(thck = ThckArea_right$FS_R_Posteriorcingulate_Thck,
                                   area = ThckArea_right$FS_R_Posteriorcingulate_Area,
                                   vol = ThckArea_right$FS_R_Posteriorcingulate_Vol)
reg_ThckArea_postcing[, outlier := as.factor(ifelse(thck < 2.2, 1, 0))]
reg_ThckArea_postcing[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
outlier_subjects_2_r = which(reg_ThckArea_postcing$outlier == 1)
# outlier_subjects_1_r == outlier_subjects_2_r

setwd("~/Desktop/Scholar_Indiana/my_work/codes")
### REPLACEMENT ###
### replacement of the outlier subjects (left) ###
all_outliers_rows = c(outlier_subjects_1_l, outlier_subjects_2_l) # outlier rows from unrelated FS_data
outliers_ids = FS_data[all_outliers_rows, ]$Subject # IDs of outlier subjects
outlier_fam_ids = subset(dmgrph_data_rel_restr, Subject %in% outliers_ids)$Family_ID # IDs of outlier subjects' families
outlier_fam_data = subset(dmgrph_data_rel_restr, Family_ID %in% outlier_fam_ids) # demographic data of outlier subjects' families
outlier_fam_mem_ids = outlier_fam_data$Subject # all subjects' IDs from outliers' families
FS_data_rel = as.data.table(read_csv("data_files/unrestricted_hcp_freesurfer.csv", col_names = TRUE))
outlier_families_FSdata = subset(FS_data_rel, Subject %in% outlier_fam_mem_ids) # demographic data of outlier subjects' families

outlier_families_AreaThck = data.table(Subject = outlier_families_FSdata$Subject,
                                       Family_ID = as.character(NA),
                                       FS_L_Isthmuscingulate_Area = outlier_families_FSdata$FS_L_Isthmuscingulate_Area,
                                       FS_L_Isthmuscingulate_Thck = outlier_families_FSdata$FS_L_Isthmuscingulate_Thck)
for (i in 1:nrow(outlier_families_AreaThck)) {
  outlier_families_AreaThck[i, Family_ID := dmgrph_data_rel_restr[Subject == outlier_families_AreaThck[i, ]$Subject]$Family_ID]
}
outlier_families_AreaThck$is_outlier = outlier_families_AreaThck$Subject %in% outliers_ids
possible_replacements_data = outlier_families_AreaThck[is_outlier == FALSE & FS_L_Isthmuscingulate_Area < 2000 & FS_L_Isthmuscingulate_Thck > 1.75]
possible_replacements_data = possible_replacements_data[order(possible_replacements_data$Family_ID), ]
replacements_data = possible_replacements_data[c(1, 3, 4, 6, 
                                                 8, 10, 12, 
                                                 15, 16, 18, 
                                                 20, 22, 23, 
                                                 25, 26, 27, 28)]
replacements_ids = replacements_data$Subject

dmgrph_data_replacements_l = subset(dmgrph_data_rel, Subject %in% replacements_ids)
FS_data_replacements_l = subset(FS_data_rel, Subject %in% replacements_ids)

FS_data_l = subset(FS_data, Subject %notin% outliers_ids)
dmgrph_data_l = subset(dmgrph_data, Subject %notin% outliers_ids)

FS_data_l = rbind(FS_data_l, FS_data_replacements_l)
dmgrph_data_l = rbind(dmgrph_data_l, dmgrph_data_replacements_l)

posteriorcingulate_indices = FS_data_l$FS_L_Posteriorcingulate_Area < 2000 # 2 outliers
FS_data_l = FS_data_l[posteriorcingulate_indices, ]
dmgrph_data_l = dmgrph_data_l[posteriorcingulate_indices, ]

precentral_indices = FS_data_l$FS_L_Precentral_Area < 8000 # 1 outlier
FS_data_l = FS_data_l[precentral_indices, ]
dmgrph_data_l = dmgrph_data_l[precentral_indices, ]

### replacement of the outlier subjects (right) ###
all_outliers_rows = outlier_subjects_1_r # outlier rows from unrelated FS_data
outliers_ids = FS_data[all_outliers_rows, ]$Subject # IDs of outlier subjects
outlier_fam_ids = subset(dmgrph_data_rel_restr, Subject %in% outliers_ids)$Family_ID # IDs of outlier subjects' families
outlier_fam_data = subset(dmgrph_data_rel_restr, Family_ID %in% outlier_fam_ids) # demographic data of outlier subjects' families
outlier_fam_mem_ids = outlier_fam_data$Subject # all subjects' IDs from outliers' families
FS_data_rel = as.data.table(read_csv("data_files/unrestricted_hcp_freesurfer.csv", col_names = TRUE))
outlier_families_FSdata = subset(FS_data_rel, Subject %in% outlier_fam_mem_ids) # demographic data of outlier subjects' families

outlier_families_AreaThck = data.table(Subject = outlier_families_FSdata$Subject,
                                       Family_ID = as.character(NA),
                                       FS_R_Caudalanteriorcingulate_Area = outlier_families_FSdata$FS_R_Caudalanteriorcingulate_Area,
                                       FS_R_Caudalanteriorcingulate_Thck = outlier_families_FSdata$FS_R_Caudalanteriorcingulate_Thck)
for (i in 1:nrow(outlier_families_AreaThck)) {
  outlier_families_AreaThck[i, Family_ID := dmgrph_data_rel_restr[Subject == outlier_families_AreaThck[i, ]$Subject]$Family_ID]
}
outlier_families_AreaThck$is_outlier = outlier_families_AreaThck$Subject %in% outliers_ids
possible_replacements_data = outlier_families_AreaThck[is_outlier == FALSE & FS_R_Caudalanteriorcingulate_Area < 1750 & FS_R_Caudalanteriorcingulate_Thck > 2]
possible_replacements_data = possible_replacements_data[order(possible_replacements_data$Family_ID), ]
replacements_data = possible_replacements_data[c(1, 3, 5, 
                                                 6, 9, 10, 
                                                 11, 14, 15)]
replacements_ids = replacements_data$Subject

dmgrph_data_replacements_r = subset(dmgrph_data_rel, Subject %in% replacements_ids)
FS_data_replacements_r = subset(FS_data_rel, Subject %in% replacements_ids)

FS_data_r = subset(FS_data, Subject %notin% outliers_ids)
dmgrph_data_r = subset(dmgrph_data, Subject %notin% outliers_ids)

FS_data_r = rbind(FS_data_r, FS_data_replacements_r)
dmgrph_data_r = rbind(dmgrph_data_r, dmgrph_data_replacements_r)

# parsorbitalis?

######## ONCE AGAIN without outliers ########
Z_l = as.matrix(FS_data_l[, ..L_regs_ord])
Z_r = as.matrix(FS_data_r[, ..R_regs_ord])

ThckArea_left = FS_data_l[, c(327:359, 395:427)] ### without insula
ThckArea_right = FS_data_r[, c(361:393, 429:461)] ### without insula

### Cortical regions Volumes ###
Vol_left = ThckArea_left[, 1:33] * ThckArea_left[, 34:66]
Vol_right = ThckArea_right[, 1:33] * ThckArea_right[, 34:66]
colnames(Vol_left) = gsub("Thck", "Vol", 
                          colnames(ThckArea_left)[1:33])
colnames(Vol_right) = gsub("Thck", "Vol", 
                           colnames(ThckArea_right)[1:33])
ThckArea_left = cbind(ThckArea_left, Vol_left)
ThckArea_right = cbind(ThckArea_right, Vol_right)

# loop for all the regions
# left hemisphere
Thck_left = ThckArea_left[, 1:33]
Area_left = ThckArea_left[, 34:66]
Volume_left = ThckArea_left[, 67:99]

# left hemisphere
Thck_right = ThckArea_right[, 1:33]
Area_right = ThckArea_right[, 34:66]
Volume_right = ThckArea_right[, 67:99]

###################### SIMULATION STUDY #########################

### signal strength in the data ###
n = 400
m = 2
rC = 0
rD = 1 - rC
k = 0.01 
p = p_l
Qc_true = Qc_l_true
Qc_obs = Qc_l_obs
Qd_true = Qd_l_true_euc
Qd_obs = Qd_l_true_euc
n_loops = 100
sigma2b = 0.1
sigma2eps = 5

b = mvrnorm(1, 
            mu = rep(0, p), 
            Sigma = sigma2b * solve(0.1 * diag(p) + 
                                      rC * Qc_true + 
                                      rD * Qd_true))

Z = matrix(rnorm(n * p), nrow = n, ncol = p)
X = matrix(rnorm(n * m), nrow = n, ncol = m)
beta = rep(0, m)
eta = Z %*% b + X %*% beta
epsilon = rnorm(n, mean = 0, sd = sqrt(sigma2eps))
y = eta + epsilon

response_data_sb001_se1 = data.frame(x = 1:400,
                                     y1 = eta,
                                     y2 = y)

response_data_sb001_se5 = data.frame(x = 1:400,
                                     y1 = eta,
                                     y2 = y)

response_data_sb01_se1 = data.frame(x = 1:400,
                                    y1 = eta,
                                    y2 = y)

response_data_sb01_se5 = data.frame(x = 1:400,
                                    y1 = eta,
                                    y2 = y)

#################################################################
sb001_se1 = ggplot(data = response_data_sb001_se1) +
  geom_point(aes(x = y1, y = y2 - y1)) +
  xlab("response without the error") + 
  ylab("the error") +
  geom_abline(intercept = 0, slope = 0, colour = "darkred") +
  ylim(-3, 3) +
  xlim(-3, 3) +
  theme_bw()

sb001_se5 = ggplot(data = response_data_sb001_se5) +
  geom_point(aes(x = y1, y = y2 - y1)) +
  xlab("response without the error") + 
  ylab("the error") +
  geom_abline(intercept = 0, slope = 0, colour = "darkred") +
  ylim(-10, 10) +
  xlim(-10, 10) +
  theme_bw()

sb01_se1 = ggplot(data = response_data_sb01_se1) +
  geom_point(aes(x = y1, y = y2 - y1)) +
  xlab("response without the error") + 
  ylab("the error") +
  geom_abline(intercept = 0, slope = 0, colour = "darkred") +
  ylim(-10, 10) +
  xlim(-10, 10) +
  theme_bw()

sb01_se5 = ggplot(data = response_data_sb01_se5) +
  geom_point(aes(x = y1, y = y2 - y1)) +
  xlab("response without the error") + 
  ylab("the error") +
  geom_abline(intercept = 0, slope = 0, colour = "darkred") +
  ylim(-15, 15) +
  xlim(-15, 15) +
  theme_bw()

grid.arrange(sb001_se1, sb001_se5, ncol = 2)
grid.arrange(sb01_se1, sb01_se5, ncol = 2)

### penalty parameters behaviour in disPEER ###
rC = c(0, 0.1, 0.5, 0.9, 1)
lambdas_data_all = data.table()

for (i in 1:length(rC)) {
  set.seed(1970)
  MSE_Bias = calc_est_MSE(n = 400, m = 2, 
                          sigma2b = 0.01,
                          sigma2eps = 1,
                          rC = rC[i], rD = 1 - rC[i],
                          k = 0.01, p = p_r,
                          Qc_true = Qc_r_true,
                          Qc_obs = Qc_r_obs,
                          Qd_true = Qd_r_true_euc,
                          Qd_obs = Qd_r_true_euc,
                          n_loops = 100)
  
  ### penalty parameters ###
  lambdas_rC = as.data.table(MSE_Bias[[6]])
  colnames(lambdas_rC) = c("lambda_C",
                           "lambda_D",
                           "lambda_R")
  lambdas_data = melt(lambdas_rC)
  colnames(lambdas_data) = c("lambda",
                             "value")
  lambdas_data$`simulation run` = rep(1:100, times = 3)
  lambdas_data$rC = rep(paste("rC =", rC[i]), times = 300)
  
  lambdas_data_all = rbind(lambdas_data_all,
                           lambdas_data)
}

ggplot(data = lambdas_data_all[lambda %in% c("lambda_C",
                                             "lambda_D")],
       aes(x = `simulation run`,
           y = value,
           colour = lambda)) +
  facet_wrap(. ~ rC) +
  geom_point() +
  geom_line() +
  theme_bw()

ggplot(data = lambdas_data_all[lambda %in% c("lambda_C",
                                             "lambda_D")],
       aes(x = lambda,
           y = value,
           colour = lambda)) +
  facet_wrap(. ~ rC) +
  geom_boxplot() +
  geom_boxplot() +
  theme_bw()

mean(lambdas_data_all[rC == 0.5 & lambda == "lambda_C"]) 
lambdas_data_all[rC == 0.5 & lambda == "lambda_D"]$value)

############################################################
rC = 0.5
set.seed(1970)
MSE_Bias = calc_est_MSE(n = 400, m = 2, 
                        sigma2b = 0.01,
                        sigma2eps = 1,
                        rC = rC, rD = 1 - rC,
                        k = 0.01, p = p_l,
                        Qc_true = Qc_l_true,
                        Qc_obs = Qc_l_obs,
                        Qd_true = Qd_l_true_euc,
                        Qd_obs = Qd_l_true_euc,
                        n_loops = 100)

### Variance ###
E_b_hat_ridge = apply(MSE_Bias[[1]], 1, mean)
E_b_hat_ols = apply(MSE_Bias[[2]], 1, mean)
E_b_hat_riPEER = apply(MSE_Bias[[3]], 1, mean)
E_b_hat_disPEER = apply(MSE_Bias[[4]], 1, mean)
b = as.vector(MSE_Bias[[5]])

Var_ridge = apply((MSE_Bias[[1]] - E_b_hat_ridge)^2, 1, mean)
Var_ols = apply((MSE_Bias[[2]] - E_b_hat_ols)^2, 1, mean)
Var_riPEER = apply((MSE_Bias[[3]] - E_b_hat_riPEER)^2, 1, mean)
Var_disPEER = apply((MSE_Bias[[4]] - E_b_hat_disPEER)^2, 1, mean)

### Bias^2 ###
Bias2_ridge = (E_b_hat_ridge - b)^2
Bias2_ols = (E_b_hat_ols - b)^2
Bias2_riPEER = (E_b_hat_riPEER - b)^2
Bias2_disPEER = (E_b_hat_disPEER - b)^2

### MSE ###
MSE_ridge = apply((MSE_Bias[[1]] - b)^2, 1, mean)
MSE_ols = apply((MSE_Bias[[2]] - b)^2, 1, mean)
MSE_riPEER = apply((MSE_Bias[[3]] - b)^2, 1, mean)
MSE_disPEER = apply((MSE_Bias[[4]] - b)^2, 1, mean)

### relative MSE ###
relMSE_ridge = apply(((MSE_Bias[[1]] - b)^2) / (b^2), 1, mean)
relMSE_ols = apply(((MSE_Bias[[2]] - b)^2) / (b^2), 1, mean)
relMSE_riPEER = apply(((MSE_Bias[[3]] - b)^2) / (b^2), 1, mean)
relMSE_disPEER = apply(((MSE_Bias[[4]] - b)^2) / (b^2), 1, mean)

############################################################

Var = data.frame(OLS = Var_ols,
                 Ridge = Var_ridge,
                 riPEER = Var_riPEER,
                 disPEER = Var_disPEER)
Bias = data.frame(OLS = Bias2_ols,
                  Ridge = Bias2_ridge,
                  riPEER = Bias2_riPEER,
                  disPEER = Bias2_disPEER)
MSE = data.frame(OLS = MSE_ols,
                 Ridge = MSE_ridge,
                 riPEER = MSE_riPEER,
                 disPEER = MSE_disPEER)
relMSE = data.frame(OLS = relMSE_ols,
                    Ridge = relMSE_ridge,
                    riPEER = relMSE_riPEER,
                    disPEER = relMSE_disPEER)

### checking MSE = BIAS^2 + VAR ###
# all(round(Var_ridge + Bias2_ridge, 4) == round(MSE_ridge, 4)) ### OK!

ggplot(data = Var, aes(x = 1:27)) +
  geom_line(aes(y = Ridge, colour = "Ridge")) +
  geom_point(aes(y = Ridge, colour = "Ridge")) +
  geom_line(aes(y = OLS, colour = "OLS")) +
  geom_point(aes(y = OLS, colour = "OLS")) +
  geom_line(aes(y = riPEER, colour = "riPEER")) +
  geom_line(aes(y = disPEER, colour = "disPEER")) +
  geom_point(aes(y = riPEER, colour = "riPEER")) +
  geom_point(aes(y = disPEER, colour = "disPEER")) +
  ylab("Variance") +
  xlab("b index") +
  ylim(c(0, 0.00405)) +
  scale_colour_manual("",
                      breaks = c("Ridge", 
                                 "OLS",
                                 "riPEER", 
                                 "disPEER"),
                      values = c("black", 
                                 "#E7B800",
                                 "firebrick2",
                                 "#00AFBB")) +
  theme_bw()
  # ggtitle(paste("Variance for b coefficients with rC equal to", rC, "with"))

 b=ggplot(data = Bias, aes(x = 1:25)) +
  geom_line(aes(y = Ridge, colour = "Ridge")) +
  geom_point(aes(y = Ridge, colour = "Ridge")) +
  geom_line(aes(y = OLS, colour = "OLS")) +
  geom_point(aes(y = OLS, colour = "OLS")) +
  geom_line(aes(y = riPEER, colour = "riPEER")) +
  geom_line(aes(y = disPEER, colour = "disPEER")) +
  geom_point(aes(y = riPEER, colour = "riPEER")) +
  geom_point(aes(y = disPEER, colour = "disPEER")) +
  ylab("Squared bias") +
  xlab("b index") +
  scale_colour_manual("",
                      breaks = c("Ridge", 
                                 "OLS",
                                 "riPEER", 
                                 "disPEER"),
                      values = c("black", 
                                 "#E7B800",
                                 "firebrick2",
                                 "#00AFBB")) +
   theme_bw()
  # ggtitle(paste("Bias for b coefficients with rC equal to", rC))

m=ggplot(data = MSE, aes(x = 1:25)) +
  geom_line(aes(y = Ridge, colour = "Ridge")) +
  geom_line(aes(y = OLS, colour = "OLS")) +
  geom_line(aes(y = riPEER, colour = "riPEER")) +
  geom_line(aes(y = disPEER, colour = "disPEER")) +
  geom_point(aes(y = Ridge, colour = "Ridge")) +
  geom_point(aes(y = OLS, colour = "OLS")) +
  geom_point(aes(y = riPEER, colour = "riPEER")) +
  geom_point(aes(y = disPEER, colour = "disPEER")) +
  ylab("MSE") +
  xlab("b index") +
  scale_colour_manual("",
                      breaks = c("Ridge", 
                                 "OLS",
                                 "riPEER", 
                                 "disPEER"),
                      values = c("black", 
                                 "#E7B800",
                                 "firebrick2",
                                 "#00AFBB")) +
  theme_bw()
  # ggtitle(paste("MSE for b coefficients with rC equal to", rC))

v
b
m
grid.arrange(v, b, m, ncol = 1)

ggplot(data = relMSE, aes(x = 1:25)) +
  geom_line(aes(y = Ridge, colour = "Ridge")) +
  geom_line(aes(y = OLS, colour = "OLS")) +
  geom_line(aes(y = riPEER, colour = "riPEER")) +
  geom_line(aes(y = disPEER, colour = "disPEER")) +
  geom_point(aes(y = Ridge, colour = "Ridge")) +
  geom_point(aes(y = OLS, colour = "OLS")) +
  geom_point(aes(y = riPEER, colour = "riPEER")) +
  geom_point(aes(y = disPEER, colour = "disPEER")) +
  ylab("relative MSE") +
  xlab("b index") +
  scale_colour_manual("",
                      breaks = c("Ridge", 
                                 "OLS",
                                 "riPEER", 
                                 "disPEER"),
                      values = c("black", 
                                 "#E7B800",
                                 "red1",
                                 "#00AFBB")) +
  theme_bw()

# boxplot_rMSE_01 = melt(relMSE)
# boxplot_rMSE_001 = melt(relMSE)
# boxplot_rMSE = rbind(boxplot_rMSE_01,
#                      boxplot_rMSE_001)
# boxplot_rMSE$signal_strength = factor(rep(c(0.1, 0.01), 
#                                           each = 25*4))
# rMSE_data = as.data.table(boxplot_rMSE)

###
maxval <- 0.8
boxplot_rMSE_001 = melt(relMSE)
boxplot_rMSE_001 = boxplot_rMSE_001[boxplot_rMSE_001$variable != "OLS", ]

p = ggplot(data = boxplot_rMSE_001, aes(y = value,
                                        x = variable,
                                        fill = variable)) 

dd_data = boxplot_rMSE_001
dd_data[dd_data$value > maxval,]$value = round(dd_data[dd_data$value > maxval,]$value, 1)

dd <- dd_data %>% filter(value > maxval) %>%
  group_by(variable) %>%
  summarise(outlier_txt = paste(value, collapse = ", "))

p2 <- p + 
  geom_violin(scale = "width") +
  geom_boxplot(colour = "white",
               width = 0.1) +
  coord_cartesian(ylim = c(0, maxval),
                  clip = "off") + 
  # coord_cartesian(ylim = c(min(boxplot_rMSE_001$value), 
  #                          maxval * 1.2),
  #                 clip = "off") +
  geom_text(data = dd, aes(y = maxval, 
                           label = str_wrap(outlier_txt, 
                                            width = 50)),
            size = 3, vjust = 0, hjust = -0.2,
            position = position_dodge(width = 1)) +
  geom_segment(data = dd, 
               aes(y = maxval * 0.9,
                   x = as.numeric(variable) - 0.85,
                   yend = maxval * 0.98,
                   xend = as.numeric(variable) - 0.85),
               arrow = arrow(length = unit(0.3, 
                                           "cm"))) +
  ylab("relative MSE") +
  xlab(" ") +
  scale_fill_manual(name = "method", values = c(
    # "#E7B800",
    "black", 
    "red1",
    "#00AFBB")) +
  theme_bw()
p2


b001 = ggplot(data = rMSE_data[signal_strength == 0.01]) +
  geom_violin(aes(y = value,
                   x = variable,
                   fill = variable),
               outlier.shape = NA) +
  ylab("relative MSE") +
  xlab(" ") +
  # coord_cartesian(ylim = c(0, 5)) +
  scale_fill_manual(name = "method", values = c("#E7B800",
                                                "black", 
                                                "red1",
                                                "#00AFBB")) +
  # facet_wrap(. ~ signal_strength, scales = "free_y") +
  theme_bw()


b01 = ggplot(data = rMSE_data[signal_strength == 0.1]) +
  geom_violin(aes(y = value,
                   x = variable,
                   fill = variable),
               outlier.shape = NA) +
  ylab("relative MSE") +
  xlab(" ") +
  # coord_cartesian(ylim = c(0, 5)) +
  scale_fill_manual(name = "method", values = c("#E7B800",
                                                "black", 
                                                "red1",
                                                "#00AFBB")) +
  # facet_wrap(. ~ signal_strength, scales = "free_y") +
  theme_bw()

grid.arrange(b001, b01, ncol = 2)

no_ols_id = boxplot_rMSE$variable != "OLS"
no_ols = c(1:216)[no_ols_id]
rMSE_data = boxplot_rMSE[no_ols, ]

ggplot(data = rMSE_data) +
  geom_boxplot(aes(y = value, 
                   x = signal_strength,
                   fill = variable)) +
  ylab("rMSE") +
  xlab("signal stregth") +
  # coord_cartesian(ylim = c(0, 2.5)) +
  scale_fill_manual(values = c("#CC0000", "#009900", "#3399FF")) +
  facet_wrap(. ~ signal_strength, scales = "free_y") +
  theme_bw()

grid.arrange(boxplot_0001, boxplot_01001)

# ggplot() +
#   geom_line(aes(x = 1:27, y = b), col = "dodgerblue") +
#   geom_point(aes(x = 1:27, y = b), col = "dodgerblue") +
#   geom_vline(xintercept = c(12, 14, 15), lty = 2, col = "red") +
#   ylab("MSE") +
#   xlab("b index") +
#   geom_hline(yintercept = 0, lty = 2) +
#   theme_bw()

c(mean(relMSE_ridge),
  mean(relMSE_riPEER),
  mean(relMSE_disPEER))
c(median(relMSE_ridge),
  median(relMSE_riPEER),
  median(relMSE_disPEER))
mean(relMSE_riPEER) / mean(relMSE_disPEER)
median(relMSE_riPEER) / median(relMSE_disPEER)

round(c(mean(MSE_ridge),
        mean(MSE_riPEER),
        mean(MSE_disPEER)), 4)
round(c(median(MSE_ridge),
        median(MSE_riPEER),
        median(MSE_disPEER)), 4)

mean(Var_riPEER > Var_disPEER)

###################### b ESTIMATES #########################
### LEFT HEMISPHERE ###
set.seed(1970)
est_results = calc_estimations(n = 400, m = 2, 
                               sigma2b = 0.01,
                               sigma2eps = 1,
                               rC = 0, rD = 1,
                               k = 0.004, p = p_r,
                               Qc_true = Qc_r_true,
                               Qc_obs = Qc_r_obs,
                               Qd_true = Qd_r_true_euc,
                               Qd_obs = Qd_r_true_euc,
                               which_hemi = "left")
# est_results$lambdas_disPEER
est_results$lambdas_disPEER
# lambda = (lambda_C, lambda_D, lambda_R)
# rC = 1, rD = 0, lambdas: 79.722053 1.513885 9.769178
# rC = 0.9, rD = 0.1, lambdas: 35.13153 39.10471 15.46462
# rC = 0.5, rD = 0.5, lambdas: 0.00001 119.221695 6.867136
# rC = 0.1, rD = 0.9, lambdas: 0.00001 103.35648 7.70194
# rC = 0, rD = 1, lambdas: 0.00001 129.27293  11.93026

rC_vec = c(0.1, 0.5, 0.9)
sigma2b_vec = c(0.1, 0.01)
sigma2eps_vec = 1

Results = matrix(NA,
                 length(rC_vec) * length(sigma2b_vec) * length(sigma2eps_vec),
                 6)

Results[, 1] = rep(rC_vec, each = 2)
Results[, 2] = rep(sigma2b_vec, times = 3)
Results[, 3] = rep(sigma2eps_vec, each = 6)

l = 0
for (i in 1:length(rC_vec)) {
  for (j in 1:length(sigma2b_vec)) {
    set.seed(1970)
    MSE_Bias = calc_est_MSE(n = 400, m = 2, 
                            sigma2b = sigma2b_vec[j],
                            sigma2eps = sigma2eps_vec,
                            rC = rC_vec[i], rD = 1 - rC_vec[i],
                            k = 0.01, p = p_l,
                            Qc_true = Qc_l_true,
                            Qc_obs = Qc_l_obs,
                            Qd_true = Qd_l_true_euc,
                            Qd_obs = Qd_l_true_euc,
                            n_loops = 100)
    
    ### relative MSE ###
    b = as.vector(MSE_Bias[[5]])
    rMSE_ridge = apply((MSE_Bias[[1]] - b)^2, 1, mean) / b^2
    rMSE_riPEER = apply((MSE_Bias[[3]] - b)^2, 1, mean) / b^2
    rMSE_disPEER = apply((MSE_Bias[[4]] - b)^2, 1, mean) / b^2
    result = c(median(rMSE_ridge),
               median(rMSE_riPEER),
               median(rMSE_disPEER))
    l = l + 1
    Results[l, 4:6] = round(result, 4)
  }
}
latextable(Results, digits = 4, scientific = F)
Results = as.data.table(Results)
colnames(Results) = c("rC", "sb", "se", 
                      "ridge", "riPEER", "disPEER")
ggplot_results = cbind(rbind(Results[, 1:3], 
                             Results[, 1:3], 
                             Results[, 1:3]),
                       melt(Results[, 4:6]))
colnames(ggplot_results) = c("rC", "sb", "se", 
                             "method", "rMSE")
ggplot_results$rC = rep(c("rC = 0.1",
                          "rC = 0.5",
                          "rC = 0.9"), 
                        each = 2,
                        times = 3)

glimpse(ggplot_results)

ggplot(data = ggplot_results, 
       aes(x = sb, y = rMSE, colour = method)) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ rC) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_x_continuous(labels = c(0.01, 0.1),
                     breaks = c(0.01, 0.1)) +
  scale_colour_manual(values = c("black",
                                 "red1",
                                 "#00AFBB"))

### MULTIPLE RESULTS ###
rC_vec = c(0.1)
sigma2b_vec = c(0.1, 0.01)
h_vec = c(1, 5, 25, 100)
sigma2eps_vec = 1

Results = matrix(NA,
                 length(rC_vec) * length(sigma2b_vec) * length(sigma2eps_vec) * length(h_vec),
                 7)

Results[, 1] = rep(h_vec, each = 2)
Results[, 2] = rep(rC_vec, times = 8)
Results[, 3] = rep(sigma2b_vec, times = 4)
Results[, 4] = rep(sigma2eps_vec, times = 8)

l = 0
for (k in 1:length(h_vec)) {
  for (i in 1:length(rC_vec)) {
    for (j in 1:length(sigma2b_vec)) {
      Ad_r_true_geo = calc_adj_matrix_d(Dnorm_r_true_geo, h = h_vec[k])
      Ad_l_true_geo = calc_adj_matrix_d(Dnorm_l_true_geo, h = h_vec[k])
      Qd_r_true_geo = calc_norm_laplacian_d(Ad_r_true_geo)
      Qd_l_true_geo = calc_norm_laplacian_d(Ad_l_true_geo)
      set.seed(1970)
      MSE_Bias = calc_est_MSE(n = 400, m = 2, 
                              sigma2b = sigma2b_vec[j],
                              sigma2eps = sigma2eps_vec,
                              rC = rC_vec[i], rD = 1 - rC_vec[i],
                              k = 0.01, p = p_r,
                              Qc_true = Qc_r_true,
                              Qc_obs = Qc_r_obs,
                              Qd_true = Qd_r_true_geo,
                              Qd_obs = Qd_r_true_geo,
                              n_loops = 100)
      
      ### relative MSE ###
      b = as.vector(MSE_Bias[[5]])
      rMSE_ridge = apply((MSE_Bias[[1]] - b)^2, 1, mean) / b^2
      rMSE_riPEER = apply((MSE_Bias[[3]] - b)^2, 1, mean) / b^2
      rMSE_disPEER = apply((MSE_Bias[[4]] - b)^2, 1, mean) / b^2
      result = c(median(rMSE_ridge),
                 median(rMSE_riPEER),
                 median(rMSE_disPEER))
      l = l + 1
      Results[l, 5:7] = round(result, 4)
    }
  }
}
latextable(Results, digits = 4, scientific = F)
colnames(Results) = c("rC", "sb", "se", 
                      "h_hat", "ridge", "riPEER", "disPEER")

Results = as.data.table(Results)
ggplot_results = cbind(rbind(Results[, 1:4], Results[, 1:4], Results[, 1:4]),
                       melt(Results[, 5:7]))
colnames(ggplot_results) = c("rC", "sb", "se", 
                             "h_hat", "method", "rMSE")
glimpse(ggplot_results)
ggplot_results$rC = rep(c("rC=0.1", "rC=0.5", "rC=0.9"), 
                        each = 2, times = 9)
ggplot_results$h_hat = rep(c("h=1", "h=10", "h=25"), 
                        each = 6, times = 3)


ggplot(data = ggplot_results, 
       aes(x = sb, y = rMSE, colour = method)) +
  geom_line() +
  geom_point() +
  facet_grid(rC ~ h_hat) +
  theme_bw() +
  scale_x_continuous(labels = c(0.01, 0.1),
                     breaks = c(0.01, 0.1)) +
  scale_colour_manual(values = c("black",
                                 "red1",
                                 "#00AFBB"), )

ggplot_results_trueAd = ggplot_results[c(1:6,
                                         19:24,
                                         37:42), ]
ggplot_results_trueAd$h_hat = rep(5, 18)
ggplot_results_trueAd$rMSE[13:18] = c(0.0326, 0.3936,
                                      0.0372, 0.2343,
                                      0.0357, 0.2854)
ggplot(data = ggplot_results_trueAd, 
       aes(x = sb, y = rMSE, colour = method)) +
  geom_line() +
  geom_point() +
  facet_grid(. ~ rC) +
  theme_bw()

results_bothAd = rbind(ggplot_results, ggplot_results_trueAd)
ggplot(data = results_bothAd, 
       aes(x = sb, y = rMSE, colour = method)) +
  geom_line() +
  geom_point() +
  facet_grid(rC ~ h_hat) +
  theme_bw()

############################################################


######################################################################
### MULTIPLE RESULTS ###
# set.seed(1997)
# rC = c(0.9, 0.5, 0.1)
# sigma2b = c(0.1, 0.01)
# sigma2eps = c(1, 5, 10)
# 
# Results = matrix(NA,
#                  length(rC) * length(sigma2b) * length(sigma2eps),
#                  6)
# Results[, 1] = rep(rC, each = 6)
# Results[, 2] = rep(sigma2b, each = 3, times = 3)
# Results[, 3] = rep(sigma2eps, 6)
# l = 0
# 
# for (i in 1:length(rC)) {
#   for (j in 1:length(sigma2b)) {
#     for (k in 1:length(sigma2eps)) {
#       mult_est_results = calc_multiple_est(n = 400, m = 2, 
#                                            sigma2b = sigma2b[j],
#                                            sigma2eps = sigma2eps[k],
#                                            rC = rC[i], rD = 1 - rC[i],
#                                            k = 0.01, p = p_l,
#                                            Qc_true = Qc_l_true,
#                                            Qc_obs = Qc_l_obs,
#                                            Qd_true = Qd_l_true_geo,
#                                            Qd_obs = Qd_l_true_geo,
#                                            n_loops = 100)
#       l = l + 1
#       Results[l, 4:6] = round(mult_est_results$mean, 4)
#     }
#   }
# }
# latextable(Results, digits = 4, scientific = F)
# 
# MSEr_result = list(Results_euc_l, Results_geo_l, Results_euc_r, Results_geo_r)
# which.max(MSEr_result[[1]][, 5] - MSEr_result[[1]][, 6])
# which.max(MSEr_result[[2]][, 5] - MSEr_result[[1]][, 6])
# which.max(MSEr_result[[3]][, 5] - MSEr_result[[1]][, 6])
# which.max(MSEr_result[[4]][, 5] - MSEr_result[[1]][, 6])
# 
# MSE1 = MSEr_result[[1]][MSEr_result[[1]][, 5] < 0.5 & MSEr_result[[1]][, 6] < 0.5, ]
# MSE2 = MSEr_result[[2]][MSEr_result[[2]][, 5] < 0.5 & MSEr_result[[2]][, 6] < 0.5, ]
# MSE3 = MSEr_result[[3]][MSEr_result[[3]][, 5] < 0.5 & MSEr_result[[3]][, 6] < 0.5, ]
# MSE4 = MSEr_result[[4]][MSEr_result[[4]][, 5] < 0.5 & MSEr_result[[4]][, 6] < 0.5, ]
# 
# which.max(MSE1[, 5] - MSE1[, 6])
# which.max(MSE2[, 5] - MSE2[, 6])
# which.max(MSE3[, 5] - MSE3[, 6])
# which.max(MSE4[, 5] - MSE4[, 6])

#### PLOTS ####
# # particular coef estimate distribution with real value marked
# ridge = ggplot() +
#   geom_histogram(aes(x = MSE_Bias[[1]][13, ]),
#                  fill = "grey",
#                  col = "black") +
#   geom_vline(xintercept = b[13], 
#              col = "red", 
#              lty = 3,
#              lwd = 1.5) +
#   xlab("ridge estimates") +
#   theme_bw()
# ols = ggplot() +
#   geom_histogram(aes(x = MSE_Bias[[2]][13, ]),
#                  fill = "grey",
#                  col = "black") +
#   geom_vline(xintercept = b[13], 
#              col = "red", 
#              lty = 3,
#              lwd = 1.5) +
#   xlab("ols estimates") +
#   xlim(-0.3, 0.3) +
#   theme_bw()
# ripeer = ggplot() +
#   geom_histogram(aes(x = MSE_Bias[[3]][13, ]),
#                  fill = "grey",
#                  col = "black") +
#   geom_vline(xintercept = b[13], 
#              col = "red", 
#              lty = 3,
#              lwd = 1.5) +
#   xlab("riPEER estimates") +
#   xlim(-0.3, 0.3) +
#   theme_bw()
# dispeer = ggplot() +
#   geom_histogram(aes(x = MSE_Bias[[4]][13, ]),
#                  fill = "grey",
#                  col = "black") +
#   geom_vline(xintercept = b[13], 
#              col = "red", 
#              lty = 3,
#              lwd = 1.5) +
#   xlab("disPEER estimates") +
#   xlim(-0.3, 0.3) +
#   theme_bw()
# 
# grid.arrange(ridge, ols, 
#              ripeer, dispeer)
# 
# ############################################################
# 
# whole b vector distribution
est_distribution = t(MSE_Bias[[1]])
colnames(est_distribution) = 1:length(b)
est_distribution = melt(est_distribution)
est_distribution = est_distribution[, -1]
colnames(est_distribution) = c("index", "estimation")
est_distribution$index = as.numeric(est_distribution$index)

ridge = ggplot() +
  geom_line(aes(x = 1:length(b), y = b),
            lwd = 1.5, col = "dodgerblue") +
  geom_point(aes(x = 1:length(b), y = b),
             lwd = 3, col = "dodgerblue") +
  geom_jitter(data = est_distribution,
              aes(x = index,
                  y = estimation),
              alpha = 0.3,
              width = 0.1,) +
  ylim(-0.5, 0.5) +
  xlab("ridge estimates") +
  geom_vline(xintercept = c(12, 14, 15), lty = 2) +
  theme_bw()

est_distribution = t(MSE_Bias[[2]])
colnames(est_distribution) = 1:length(b)
est_distribution = melt(est_distribution)
est_distribution = est_distribution[, -1]
colnames(est_distribution) = c("index", "estimation")
est_distribution$index = as.numeric(est_distribution$index)

ols = ggplot() +
  geom_line(aes(x = 1:length(b), y = b),
            lwd = 1.5, col = "dodgerblue") +
  geom_point(aes(x = 1:length(b), y = b),
             lwd = 3, col = "dodgerblue") +
  geom_jitter(data = est_distribution,
              aes(x = index,
                  y = estimation),
              alpha = 0.3,
              width = 0.1) +
  ylim(-0.5, 0.5) +
  xlab("ols estimates") +
  geom_vline(xintercept = c(12, 14, 15), lty = 2) +
  theme_bw()

est_distribution = t(MSE_Bias[[3]])
colnames(est_distribution) = 1:length(b)
est_distribution = melt(est_distribution)
est_distribution = est_distribution[, -1]
colnames(est_distribution) = c("index", "estimation")
est_distribution$index = as.numeric(est_distribution$index)

ripeer = ggplot() +
  geom_line(aes(x = 1:length(b), y = b),
            lwd = 1.5, col = "dodgerblue") +
  geom_point(aes(x = 1:length(b), y = b),
             lwd = 3, col = "dodgerblue") +
  geom_jitter(data = est_distribution,
              aes(x = index,
                  y = estimation),
              alpha = 0.3,
              width = 0.1) +
  ylim(-0.5, 0.5) +
  xlab("riPEER estimates") +
  geom_vline(xintercept = c(12, 14, 15), lty = 2) +
  theme_bw()

est_distribution = t(MSE_Bias[[4]])
colnames(est_distribution) = 1:length(b)
est_distribution = melt(est_distribution)
est_distribution = est_distribution[, -1]
colnames(est_distribution) = c("index", "estimation")
est_distribution$index = as.numeric(est_distribution$index)

dispeer = ggplot() +
  geom_line(aes(x = 1:length(b), y = b),
            lwd = 1.5, col = "dodgerblue") +
  geom_point(aes(x = 1:length(b), y = b),
             lwd = 3, col = "dodgerblue") +
  geom_jitter(data = est_distribution,
              aes(x = index,
                  y = estimation),
              alpha = 0.3,
              width = 0.1) +
  ylim(-0.5, 0.5) +
  xlab("disPEER estimates") +
  geom_vline(xintercept = c(12, 14, 15), lty = 2) +
  theme_bw()

grid.arrange(ridge, ols,
             ripeer, dispeer)



D_G_l = apply(D_geodesic_left_true, 1, as.numeric)
D_E_l = apply(D_euclidean_left_true, 1, as.numeric)
D_G_r = apply(D_geodesic_right_true, 1, as.numeric)
D_E_r = apply(D_euclidean_right_true, 1, as.numeric)
Ad_E_l = apply(Ad_l_true_euc, 1, as.numeric)
Ad_E_r = apply(Ad_r_true_euc, 1, as.numeric)
Ad_G_l = apply(Ad_l_true_geo, 1, as.numeric)
Ad_G_r = apply(Ad_r_true_geo, 1, as.numeric)
Ac_l = apply(Ac_l_true, 1, as.numeric)
Ac_r = apply(Ac_r_true, 1, as.numeric)



correlations_A_l = sapply(1:nrow(Ac_l), function(i) {
  cor(Ad_E_l[i, ], Ac_l[i, ])
})

correlations_A_r = sapply(1:nrow(Ac_r), function(i) {
  cor(Ad_E_r[i, ], Ac_r[i, ])
})

par(mfrow = c(1, 2))
plot(correlations_A_l, ylim = c(-1, 1))
plot(correlations_A_r, ylim = c(-1, 1))


correlations_l = sapply(1:nrow(D_G_l), function(i) {
  cor(D_G_l[i, ], D_E_l[i, ])
})

correlations_r = sapply(1:nrow(D_G_r), function(i) {
  cor(D_G_r[i, ], D_E_r[i, ])
})

mean_cor_l = mean(sapply(1:nrow(D_G_l), function(i) {
  cor(D_G_l[i, ], D_E_l[i, ])
}))

mean_cor_r = mean(sapply(1:nrow(D_G_r), function(i) {
  cor(D_G_r[i, ], D_E_r[i, ])
}))

regions_l = rownames(Qc_l_true)
regions_r = rownames(Qc_r_true)

cor_data_l = data.table(indexes = 1:p_l,
                        correlations_l = correlations_l)
cor_data_r = data.table(indexes = 1:p_r,
                        correlations_r = correlations_r)

ggplot(data = cor_data_l, aes(x = indexes, y = correlations_l)) +
  geom_point() +
  geom_line(lty = 2) +
  geom_hline(yintercept = mean_cor_l, col = "red", lty = 2) +
  scale_x_continuous(breaks = 1:p_l,
                     labels = rownames(Qc_l_true)) +
  ylab("Correlation") +
  xlab("Cortical regions") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   vjust = 0.3)) 


ggplot(data = cor_data_r, aes(x = indexes, y = correlations_r)) +
  geom_point() +
  geom_line(lty = 2) +
  geom_hline(yintercept = mean_cor_r, col = "red", lty = 2) +
  scale_x_continuous(breaks = 1:p_r,
                     labels = rownames(Qc_r_true)) +
  ylab("Correlation") +
  xlab("Cortical regions") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 9,
                                   angle = 90,
                                   vjust = 0.3)) 

############################################################
L_regs_ord_vol = c("FS_L_Cuneus_Vol",
                   "FS_L_Fusiform_Vol",
                   "FS_L_Lateraloccipital_Vol",
                   "FS_L_Lingual_Vol",
                   "FS_L_Pericalcarine_Vol",
                   "FS_L_Paracentral_Vol",
                   "FS_L_Postcentral_Vol",
                   "FS_L_Precentral_Vol",
                   "FS_L_Superiortemporal_Vol",
                   "FS_L_Transversetemporal_Vol",
                   "FS_L_Entorhinal_Vol",
                   "FS_L_Inferiortemporal_Vol",
                   "FS_L_Lateralorbitofrontal_Vol",
                   "FS_L_Medialorbitofrontal_Vol",
                   "FS_L_Frontalpole_Vol",
                   "FS_L_Temporalpole_Vol",
                   "FS_L_Bankssts_Vol",
                   "FS_L_Caudalmiddlefrontal_Vol",
                   "FS_L_Inferiorparietal_Vol",
                   "FS_L_Isthmuscingulate_Vol",
                   "FS_L_Middletemporal_Vol",
                   "FS_L_Parsorbitalis_Vol",
                   "FS_L_Parstriangularis_Vol" ,
                   "FS_L_Posteriorcingulate_Vol",
                   "FS_L_Precuneus_Vol",
                   "FS_L_Rostralanteriorcingulate_Vol",
                   "FS_L_Superiorfrontal_Vol")

# R_regs_ord_vol = c("FS_R_Cuneus_Vol",
#                    "FS_R_Fusiform_Vol",
#                    "FS_R_Lateraloccipital_Vol",
#                    "FS_R_Lingual_Vol",
#                    "FS_R_Pericalcarine_Vol",
#                    "FS_R_Caudalanteriorcingulate_Vol",
#                    "FS_R_Supramarginal_Vol" ,
#                    "FS_R_Entorhinal_Vol",
#                    "FS_R_Inferiortemporal_Vol",
#                    "FS_R_Lateralorbitofrontal_Vol",
#                    "FS_R_Medialorbitofrontal_Vol",
#                    "FS_R_Frontalpole_Vol",
#                    "FS_R_Temporalpole_Vol",
#                    "FS_R_Caudalmiddlefrontal_Vol",
#                    "FS_R_Parsopercularis_Vol",
#                    "FS_R_Parstriangularis_Vol",
#                    "FS_R_Rostralmiddlefrontal_Vol",
#                    "FS_R_Bankssts_Vol",
#                    "FS_R_Inferiorparietal_Vol",
#                    "FS_R_Isthmuscingulate_Vol",
#                    "FS_R_Middletemporal_Vol",
#                    "FS_R_Parsorbitalis_Vol",
#                    "FS_R_Precuneus_Vol",
#                    "FS_R_Rostralanteriorcingulate_Vol",
#                    "FS_R_Superiorfrontal_Vol")

R_regs_ord_vol = c("FS_R_Cuneus_Vol",
                   "FS_R_Fusiform_Vol",
                   "FS_R_Lateraloccipital_Vol",
                   "FS_R_Lingual_Vol",
                   "FS_R_Pericalcarine_Vol",
                   "FS_R_Paracentral_Vol",
                   "FS_R_Postcentral_Vol",
                   "FS_R_Posteriorcingulate_Vol",
                   "FS_R_Precentral_Vol",
                   "FS_R_Superiortemporal_Vol",
                   "FS_R_Transversetemporal_Vol",
                   "FS_R_Entorhinal_Vol",
                   "FS_R_Inferiortemporal_Vol",
                   "FS_R_Lateralorbitofrontal_Vol",
                   "FS_R_Medialorbitofrontal_Vol",
                   "FS_R_Frontalpole_Vol",
                   "FS_R_Temporalpole_Vol",
                   "FS_R_Caudalmiddlefrontal_Vol",
                   "FS_R_Parsopercularis_Vol",
                   "FS_R_Parstriangularis_Vol",
                   "FS_R_Rostralmiddlefrontal_Vol",
                   "FS_R_Bankssts_Vol",
                   "FS_R_Inferiorparietal_Vol",
                   "FS_R_Isthmuscingulate_Vol",
                   "FS_R_Middletemporal_Vol",
                   "FS_R_Parsorbitalis_Vol",
                   "FS_R_Precuneus_Vol",
                   "FS_R_Rostralanteriorcingulate_Vol",
                   "FS_R_Superiorfrontal_Vol")


Z_vol_l = Volume_left[, ..L_regs_ord_vol]
Z_vol_r = Volume_right[, ..R_regs_ord_vol]

y_l = as.matrix(dmgrph_data_l$Flanker_Unadj) 
y_l = as.matrix(dmgrph_data_l$PicSeq_Unadj) # $1 \cdot 10^{-5}$ & 1.000039 \cdot 10^{-5}
y_l = as.matrix(dmgrph_data_l$ReadEng_Unadj) # $1 \cdot 10^{-5}$ & $89.2627$
y_l = as.matrix(dmgrph_data_l$ListSort_Unadj) # 2.877 \cdot 10^2 & 1.000011 \cdot 10^{-5}
y_l = as.matrix(dmgrph_data_l$ProcSpeed_Unadj)
y_l = as.matrix(dmgrph_data_l$CardSort_Unadj)
y_l = as.matrix(dmgrph_data_l$PicVocab_Unadj) # $1 \cdot 10^{-5}$ & $1048.97553$


y_r = as.matrix(dmgrph_data_r$Flanker_Unadj) 
y_r = as.matrix(dmgrph_data_r$PicSeq_Unadj) # $1 \cdot 10^{-5}$ & 1.000039 \cdot 10^{-5}
y_r = as.matrix(dmgrph_data_r$ReadEng_Unadj) # $1 \cdot 10^{-5}$ & $89.2627$
y_r = as.matrix(dmgrph_data_r$ListSort_Unadj) # 2.877 \cdot 10^2 & 1.000011 \cdot 10^{-5}
y_r = as.matrix(dmgrph_data_r$ProcSpeed_Unadj)
y_r = as.matrix(dmgrph_data_r$CardSort_Unadj)
y_r = as.matrix(dmgrph_data_r$PicVocab_Unadj) # $1 \cdot 10^{-5}$ & $1048.97553$

y = y_r

# ggplot(data = dmgrph_data_l, aes(x = Age, y = y_l)) +
#   geom_boxplot()
# ggplot(data = dmgrph_data_r, aes(x = Age, y = y_r)) +
#   geom_boxplot()

Z_vol_norm_l = Z_vol_l / FS_data_l$FS_LCort_GM_Vol
Z_vol_norm_r = Z_vol_r / FS_data_r$FS_RCort_GM_Vol

Z_vol_norm = Z_vol_norm_r

X_l = data.table(
  gender = 1 * (dmgrph_data_l$Gender == "F"),
  # edu = subset(dmgrph_data_rel_restr, Subject %in% dmgrph_data$Subject)$SSAGA_Educ,
  age = subset(dmgrph_data_rel_restr, 
               Subject %in% dmgrph_data_l$Subject)$Age)

X_r = data.table(
  gender = 1 * (dmgrph_data_r$Gender == "F"),
  # edu = subset(dmgrph_data_rel_restr, Subject %in% dmgrph_data$Subject)$SSAGA_Educ,
  age = subset(dmgrph_data_rel_restr, 
               Subject %in% dmgrph_data_r$Subject)$Age)

X = X_r

### Estimation ###
### left ###
n = nrow(X)
m = ncol(X)
X = as.matrix(X)
X_norm = scale(X)
Z = scale(as.matrix(Z_vol_norm))
y = scale(y)

Pcx = diag(n) - X_norm %*% solve(t(X_norm) %*% X_norm) %*% t(X_norm)
y_P = Pcx %*% y
Z_P = Pcx %*% Z
Qc_obs = Qc_r_obs
Qc_true = Qc_r_true
Qd_true = Qd_r_true_euc
p = nrow(Qd_true)

set.seed(1971)
### disPEER estimation ###
l_tilde_function = function(lambdas) {
  lambdaC = lambdas[1]
  lambdaD = lambdas[2]
  lambdaR = lambdas[3]
  B_lambda =
    lambdaC * Qc_obs +
    lambdaD * Qd_true +
    lambdaR * diag(p)
  
  n * log(sum(y_P^2) - t(y_P) %*% Z_P %*% solve(B_lambda + t(Z_P) %*% Z_P) %*% t(Z_P) %*% y_P) +
    log(det(B_lambda + t(Z_P) %*% Z_P)) - 
    log(det(B_lambda))
}

lambdas_hat = sbplx(x0 = c(1, 1, 1),
                    fn = l_tilde_function,
                    lower = c(10^(-5), 10^(-5), 10^(-5)),
                    upper = c(1e6, 1e6, 1e6))
lambda_C = lambdas_hat$par[1]
lambda_D = lambdas_hat$par[2]
lambda_R = lambdas_hat$par[3]
lambdas_hat$par

b_disPEER_geo_r = as.vector(solve(t(Z_P) %*% Z_P + 
                              lambda_C * Qc_obs +
                              lambda_D * Qd_true +
                              lambda_R * diag(p)) %*% t(Z_P) %*% y_P)

### CI for disPEER ###
Blambda_opt = lambda_C * Qc_obs +
  lambda_D * Qd_true +
  lambda_R * diag(p)

b.est.boot <- function(data, indices) {
  data.boot = data[indices, ]
  y_P_boot = as.matrix(data.boot[, 1])
  Z_P_boot = as.matrix(data.boot[, 2:ncol(data)])
  b.est.dis.boot <- as.vector(solve(t(Z_P_boot) %*% Z_P_boot +
                                      Blambda_opt) %*% t(Z_P_boot) %*% y_P_boot)
  return(as.vector(b.est.dis.boot))
}
boot.out <- boot::boot(data = cbind(y_P, Z_P), statistic = b.est.boot,
                       R = 500, parallel = "multicore", ncpus = 7)
boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx) {
  boot.ci.out <- (boot::boot.ci(boot.out, type = "bca", index = idx, conf = 0.95))$bca
  return(boot.ci.out[c(4, 5)])
})
ConfBand_disPEER_euc_r <- do.call(rbind.data.frame, boot.CI.l)
names(ConfBand_disPEER_euc_r) <- c("lower", "upper")

significant_disPEER_euc_r = (1:p)[!(0 >= ConfBand_disPEER_euc_r$lower & 0 <= ConfBand_disPEER_euc_r$upper)]
significant_disPEER_euc_r

###############################################################
###############################################################

### Multiple linear regression ###
lm_model = lm(y ~ cbind(X_norm, Z))
b_lm = as.vector(coef(lm_model)[(m + 2):(p + m + 1)])
beta_lm = as.vector(coef(lm_model)[2:(m + 1)])

### Ridge estimation: without any additional information ###
### left ###
ridge_cv = cv.glmnet(x = cbind(X_norm, Z), y = y, 
                     alpha = 0, intercept = FALSE)
b_ridge = as.vector(coef(ridge_cv)[(m + 2):(p + m + 1)]) # exclude intercept and covs
beta_ridge = as.vector(coef(ridge_cv)[2:(m + 1)])
            
### disPEER estimation ###
l_tilde_function = function(lambdas) {
  lambdaC = lambdas[1]
  lambdaD = lambdas[2]
  lambdaR = lambdas[3]
  B_lambda =
    lambdaC * Qc_obs +
    lambdaD * Qd_true +
    lambdaR * diag(p)
  
  n * log(sum(y_P^2) - t(y_P) %*% Z_P %*% solve(B_lambda + t(Z_P) %*% Z_P) %*% t(Z_P) %*% y_P) +
    log(det(B_lambda + t(Z_P) %*% Z_P)) - 
    log(det(B_lambda))
}

lambdas_hat = sbplx(x0 = c(1, 1, 1),
                    fn = l_tilde_function,
                    lower = c(10^(-5), 10^(-5), 10^(-5)),
                    upper = c(1e6, 1e6, 1e6))
lambda_C = lambdas_hat$par[1]
lambda_D = lambdas_hat$par[2]
lambda_R = lambdas_hat$par[3]
lambdas_hat$par

b_disPEER = as.vector(solve(t(Z_P) %*% Z_P + 
                              lambda_C * Qc_obs +
                              lambda_D * Qd_true +
                              lambda_R * diag(p)) %*% t(Z_P) %*% y_P)
beta_disPEER = as.vector(solve(t(X_norm) %*% X_norm) %*% t(X_norm) %*% (y - Z %*% b_disPEER))

### INFORMATION FROM THE CONNECTIONS ONLY ###
### We have here y, Z, Qc and Qd_r and Qd_l ###
riPEER_1 = riPEER(Q = Qc_obs, 
                  y = y,
                  Z = Z,
                  X = X_norm,
                  compute.boot.CI = TRUE,
                  boot.R = 500)
b_riPEER = riPEER_1$b.est
beta_riPEER = riPEER_1$beta.est
ConfBand_riPEER <- riPEER_1$boot.CI

ggplot() +
  geom_line(aes(x = 1:p, y = b_riPEER, colour = "riPEER")) +
  geom_line(aes(x = 1:p, y = b_ridge, colour = "ridge")) +
  geom_line(aes(x = 1:p, y = b_lm, colour = "OLS")) +
  geom_line(aes(x = 1:p, y = b_disPEER, colour = "disPEER")) +
  geom_point(aes(x = 1:p, y = b_riPEER, colour = "riPEER")) +
  geom_point(aes(x = 1:p, y = b_ridge, colour = "ridge")) +
  geom_point(aes(x = 1:p, y = b_lm, colour = "OLS")) +
  geom_point(aes(x = 1:p, y = b_disPEER, colour = "disPEER")) +
  scale_colour_manual(" ",
                      breaks = c("riPEER", "ridge", "OLS", "disPEER"),
                      values = c("red", "black", "orange", "dodgerblue")) +
  scale_x_continuous(breaks = 1:p,
                     labels = colnames(Ad_r_true_geo)) +
  ylab("Estimates") +
  xlab("b index") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "bold", angle = 45)) 

ggplot() +
  geom_line(aes(x = 1:p, y = b_riPEER, colour = "riPEER")) +
  geom_line(aes(x = 1:p, y = b_ridge, colour = "ridge")) +
  geom_line(aes(x = 1:p, y = b_disPEER, colour = "disPEER")) +
  geom_point(aes(x = 1:p, y = b_riPEER, colour = "riPEER")) +
  geom_point(aes(x = 1:p, y = b_ridge, colour = "ridge")) +
  geom_point(aes(x = 1:p, y = b_disPEER, colour = "disPEER")) +
  scale_colour_manual(" ",
                      breaks = c("riPEER", "ridge", "OLS", "disPEER"),
                      values = c("red", "black", "orange", "dodgerblue")) +
  ylab("Estimates") +
  xlab("b index") +
  theme_bw()

latextable(data.table(b_lm,
           b_ridge,
           b_riPEER,
           b_disPEER), digits = 4)

set.seed(1970)
### CI for OLS ###
b.est.boot <- function(data, indices) {
  data.boot = data[indices, ]
  y_boot = as.matrix(data.boot[, 1])
  Z_boot = as.matrix(data.boot[, 2:(p + 1)])
  X_boot = as.matrix(data.boot[, (p + 2):ncol(data)])
  lm_model = lm(y_boot ~ cbind(X_boot, Z_boot))
  b.est.gR.boot = as.vector(coef(lm_model)[(m + 2):(p + m + 1)]) 
  return(as.vector(b.est.gR.boot))
}
boot.out <- boot::boot(data = cbind(y, Z, X_norm), statistic = b.est.boot,
                       R = 500, parallel = "multicore", ncpus = 7)
boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx) {
  boot.ci.out <- (boot::boot.ci(boot.out, type = "bca", index = idx, conf = 0.95))$bca
  return(boot.ci.out[c(4, 5)])
})
ConfBand_OLS <- do.call(rbind.data.frame, boot.CI.l)
names(ConfBand_OLS) <- c("lower", "upper")

### CI for ridge ###
b.est.boot <- function(data, indices) {
  data.boot = data[indices, ]
  y_boot = as.matrix(data.boot[, 1])
  Z_boot = as.matrix(data.boot[, 2:(p + 1)])
  X_boot = as.matrix(data.boot[, (p + 2):ncol(data)])
  ridge_cv = cv.glmnet(x = cbind(X_boot, Z_boot), 
                       y = y_boot, 
                       alpha = 0, intercept = FALSE)
  b.est.gR.boot = as.vector(coef(ridge_cv)[(m + 2):(p + m + 1)]) 
  return(as.vector(b.est.gR.boot))
}
boot.out <- boot::boot(data = cbind(y, Z, X_norm), statistic = b.est.boot,
                       R = 500, parallel = "multicore", ncpus = 7)
boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx) {
  boot.ci.out <- (boot::boot.ci(boot.out, type = "bca", index = idx, conf = 0.95))$bca
  return(boot.ci.out[c(4, 5)])
})
ConfBand_ridge <- do.call(rbind.data.frame, boot.CI.l)
names(ConfBand_ridge) <- c("lower", "upper")

### CI for disPEER ###
Blambda_opt = lambda_C * Qc_obs +
  lambda_D * Qd_true +
  lambda_R * diag(p)

b.est.boot <- function(data, indices) {
  data.boot = data[indices, ]
  y_P_boot = as.matrix(data.boot[, 1])
  Z_P_boot = as.matrix(data.boot[, 2:ncol(data)])
  b.est.dis.boot <- as.vector(solve(t(Z_P_boot) %*% Z_P_boot +
                                      Blambda_opt) %*% t(Z_P_boot) %*% y_P_boot)
  return(as.vector(b.est.dis.boot))
}
boot.out <- boot::boot(data = cbind(y_P, Z_P), statistic = b.est.boot,
                       R = 500, parallel = "multicore", ncpus = 7)
boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx) {
  boot.ci.out <- (boot::boot.ci(boot.out, type = "bca", index = idx, conf = 0.95))$bca
  return(boot.ci.out[c(4, 5)])
})
ConfBand_disPEER_euc_r <- do.call(rbind.data.frame, boot.CI.l)
names(ConfBand_disPEER) <- c("lower", "upper")

significant_riPEER = (1:p)[!(0 >= ConfBand_riPEER$lower & 0 <= ConfBand_riPEER$upper)]
significant_ridge = (1:p)[!(0 >= ConfBand_ridge$lower & 0 <= ConfBand_ridge$upper)]
significant_OLS = (1:p)[!(0 >= ConfBand_OLS$lower & 0 <= ConfBand_OLS$upper)]
significant_disPEER_geo_r = (1:p)[!(0 >= ConfBand_disPEER_geo_r$lower & 0 <= ConfBand_disPEER_geo_r$upper)]
significant_disPEER_euc_r = (1:p)[!(0 >= ConfBand_disPEER_euc_r$lower & 0 <= ConfBand_disPEER_euc_r$upper)]


# estimations_data_b = data.table(indexes = 1:p,
#                                 b_ridge = b_ridge,
#                                 b_riPEER = riPEER_1$b.est,
#                                 b_OLS = b_lm,
#                                 lower_riPEER = ConfBand_riPEER$lower,
#                                 upper_riPEER = ConfBand_riPEER$upper,
#                                 lower_ridge = ConfBand_ridge$lower,
#                                 upper_ridge = ConfBand_ridge$upper,
#                                 lower_OLS = ConfBand_OLS$lower,
#                                 upper_OLS = ConfBand_OLS$upper,
#                                 lower_disPEER_geo_l = ConfBand_disPEER_geo_r$lower,
#                                 upper_disPEER_geo_l = ConfBand_disPEER_geo_r$upper,
#                                 lower_disPEER_euc_l = ConfBand_disPEER_euc_r$lower,
#                                 upper_disPEER_euc_l = ConfBand_disPEER_euc_r$upper)

estimations_data_b_r = data.table(indexes = 1:p,
                                  b_ridge = b_ridge,
                                  b_riPEER = riPEER_1$b.est,
                                  b_OLS = b_lm,
                                  lower_riPEER = ConfBand_riPEER$lower,
                                  upper_riPEER = ConfBand_riPEER$upper,
                                  lower_ridge = ConfBand_ridge$lower,
                                  upper_ridge = ConfBand_ridge$upper,
                                  lower_OLS = ConfBand_OLS$lower,
                                  upper_OLS = ConfBand_OLS$upper,
                                  lower_disPEER_euc_r = ConfBand_disPEER_euc_r$lower,
                                  upper_disPEER_euc_r = ConfBand_disPEER_euc_r$upper,
                                  lower_disPEER_geo_r = ConfBand_disPEER_geo_r$lower,
                                  upper_disPEER_geo_r = ConfBand_disPEER_geo_r$upper)

# estimations_data_b_r = data.table(indexes = 1:p,
#                                   b_ridge = b_ridge,
#                                   b_riPEER = riPEER_1$b.est,
#                                   b_OLS = b_lm,
#                                   lower_riPEER = ConfBand_riPEER$lower,
#                                   upper_riPEER = ConfBand_riPEER$upper,  
#                                   lower_ridge = ConfBand_ridge$lower,
#                                   upper_ridge = ConfBand_ridge$upper,
#                                   lower_OLS = ConfBand_OLS$lower,
#                                   upper_OLS = ConfBand_OLS$upper,
#                                   lower_disPEER_geo_r = ConfBand_disPEER$lower,
#                                   upper_disPEER_geo_r = ConfBand_disPEER$upper)


estimation_plot_riPEER = ggplot(data = estimations_data_b_r, 
                                aes(x = indexes)) +
  geom_line(aes(y = b_riPEER), lty = 2) + 
  geom_point(aes(y = b_riPEER), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_riPEER, ymax = upper_riPEER), 
              alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_riPEER, col = "red", lty = 2) +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("riPEER \n ") +
  xlab(" ") +
  ylim(-0.35, 0.35) +
  theme_bw()

estimation_plot_ridge = ggplot(data = estimations_data_b_r, 
                               aes(x = indexes)) +
  geom_line(aes(y = b_ridge), lty = 2) +
  geom_point(aes(y = b_ridge), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_ridge, ymax = upper_ridge), 
              alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_ridge, col = "red", lty = 2) +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("ridge \n ") +
  xlab(" ") +
  ylim(-0.35, 0.35) +
  theme_bw()

estimation_plot_OLS = ggplot(data = estimations_data_b_r, 
                             aes(x = indexes)) +
  geom_line(aes(y = b_OLS), lty = 2) +
  geom_point(aes(y = b_OLS), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_OLS, ymax = upper_OLS), 
              alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_OLS, col = "red", lty = 2) +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("OLS \n ") +
  xlab(" ") +
  ylim(-0.35, 0.35) +
  theme_bw()

estimation_plot_disPEER_euc_r = ggplot(data = estimations_data_b_r, aes(x = indexes)) +
  geom_line(aes(y = b_disPEER_euc_r), lty = 2) +
  geom_point(aes(y = b_disPEER_euc_r), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_disPEER_euc_r, ymax = upper_disPEER_euc_r), 
              alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_disPEER_euc_r, col = "red", lty = 2) +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("disPEER \n (Euclidean right)") +
  xlab(" ") +
  ylim(-0.35, 0.35) +
  theme_bw() 

estimation_plot_disPEER_geo_r = ggplot(data = estimations_data_b_r, aes(x = indexes)) +
  geom_line(aes(y = b_disPEER_geo_r), lty = 2) +
  geom_point(aes(y = b_disPEER_geo_r), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_disPEER_geo_r, ymax = upper_disPEER_geo_r), 
              alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_disPEER_geo_r, col = "red", lty = 2) +
  ylab("disPEER \n (geodesic right)") +
  scale_x_continuous(breaks = 1:p,
                     labels = colnames(Ad_r_true_geo)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.6,
                                   face = "bold")) +
  ylim(-0.35, 0.35)

estimation_plot_disPEER_euc_l = ggplot(data = estimations_data_b, aes(x = indexes)) +
  geom_line(aes(y = b_disPEER_euc_l), lty = 2) +
  geom_point(aes(y = b_disPEER_euc_l), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_disPEER_euc_l, ymax = upper_disPEER_euc_l), 
              alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_disPEER_euc_l, col = "red", lty = 2) +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("disPEER \n (Euclidean left)") +
  xlab(" ") +
  ylim(-0.35, 0.35) +
  theme_bw() 

estimation_plot_disPEER_geo_l = ggplot(data = estimations_data_b, aes(x = indexes)) +
  geom_line(aes(y = b_disPEER_geo_l), lty = 2) +
  geom_point(aes(y = b_disPEER_geo_l), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_disPEER_geo_l, ymax = upper_disPEER_geo_l), 
              alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_disPEER_geo_l, col = "red", lty = 2) +
  ylab("disPEER \n (geodesic left)") +
  xlab("cortical regions") +
  scale_x_continuous(breaks = 1:27,
                     labels = colnames(Ad_l_true_geo)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.6,
                                   face = "bold")) +
  ylim(-0.35, 0.35)


grid.arrange(estimation_plot_OLS, 
             estimation_plot_ridge, 
             estimation_plot_riPEER, 
             estimation_plot_disPEER_euc_l,
             estimation_plot_disPEER_geo_l,
             ncol = 1,
             heights = c(2, 2, 2, 2, 5)) 

grid.arrange(estimation_plot_OLS, 
             estimation_plot_ridge, 
             estimation_plot_riPEER, 
             estimation_plot_disPEER_euc_r,
             estimation_plot_disPEER_geo_r,
             ncol = 1,
             heights = c(2, 2, 2, 2, 5)) 





estimations_data_y = data.table(indexes = 1:30,
                                y = y[1:30],
                                y_ridge = y_ridge[1:30],
                                y_riPEER = y_riPEER[1:30],
                                y_disPEER = y_disPEER[1:30])

ggplot(data = estimations_data_y, aes(x = indexes)) +
  geom_line(aes(y = y, colour = "y real values"), lty = 2) +
  geom_point(aes(y = y, colour = "y real values"), lwd = 1.7) +
  geom_line(aes(y = y_ridge, colour = "Ridge estimates"), lty = 2) +
  geom_point(aes(y = y_ridge, colour = "Ridge estimates"), lwd = 1.7) +
  geom_line(aes(y = y_riPEER, colour = "riPEER estimates"), lty = 2) + 
  geom_point(aes(y = y_riPEER, colour = "riPEER estimates"), lwd = 1.7) +
  geom_line(aes(y = y_disPEER, colour = "disPEER estimates"), lty = 2) +
  geom_point(aes(y = y_disPEER, colour = "disPEER estimates"), lwd = 1.7) +
  scale_colour_manual("",
                      breaks = c("y real values",
                                 "Ridge estimates",
                                 "riPEER estimates", 
                                 "disPEER estimates"),
                      values = c("black",
                                 "#E7B800",
                                 "#00AFBB", 
                                 "red1")) +
  ylab("y") +
  scale_x_continuous(breaks = 1:length(y)) +
  ggtitle("Estimations of y for the left hemisphere")

### Training & Test ###
n = length(y)
training_ids = sample(n, n / 2)
test_ids = (1:n)[-training_ids]

y_train = y[training_ids]
y_test = y[test_ids]
X_train = X[training_ids, ]
X_test = X[test_ids, ]

n_1 = length(y_train)
Z_train = scale(Z_l[training_ids, ])
Z_test = scale(Z_l[test_ids, ])
Qc_obs = Qc_l_obs
Qc_true = Qc_l_true
Qd_true = Qd_l_true_geo
p = nrow(Qd_true)
Pcx_train = diag(n_1) - X_train %*% solve(t(X_train) %*% X_train) %*% t(X_train)
y_P_train = Pcx_train %*% y_train
Z_P_train = Pcx_train %*% Z_train

### Multiple linear regression ###
lm_model_train = lm(y_train ~ cbind(X_train, Z_train))
b_lm_train = as.vector(coef(lm_model_train)[(m + 2):(p + m + 1)])
beta_lm_train = as.vector(coef(lm_model_train)[2:(m + 1)])
y_lm_train = coef(lm_model_train)[1] + X_train %*% beta_lm_train + Z_train %*% b_lm_train
y_lm_test = coef(lm_model_train)[1] + X_test %*% beta_lm_train + Z_test %*% b_lm_train

### Ridge estimation: without any additional information ###
### left ###
ridge_cv_train = cv.glmnet(x = cbind(X_train, Z_train), y = y_train, 
                           alpha = 0, intercept = FALSE)
b_ridge_train = coef(ridge_cv_train)[(m + 2):(p + m + 1)] # exclude intercept and covs
beta_ridge_train = coef(ridge_cv_train)[2:(m + 1)]
y_ridge_train = X_train %*% beta_ridge_train + Z_train %*% b_ridge_train
y_ridge_test = X_test %*% beta_ridge_train + Z_test %*% b_ridge_train
calc_MSEr(y_ridge_test, y_test)

l_tilde_function = function(lambdas) {
  lambdaC = lambdas[1]
  lambdaD = lambdas[2]
  lambdaR = lambdas[3]
  B_lambda =
    lambdaC * Qc_obs +
    lambdaD * Qd_true +
    lambdaR * diag(p)
  n * log(sum(y_P_train^2) - t(y_P_train) %*% Z_P_train %*% solve(B_lambda + t(Z_P_train) %*% Z_P_train) %*% t(Z_P_train) %*% y_P_train) + 
    log(det(B_lambda + t(Z_P_train) %*% Z_P_train)) - log(det(B_lambda))
}

lambdas_hat_train = sbplx(x0 = c(1, 1, 1),
                          fn = l_tilde_function,
                          lower = c(10^(-5), 10^(-5), 10^(-5)),
                          upper = c(1e6, 1e6, 1e6))
lambdas_hat_train$par
lambda_C_train = lambdas_hat_train$par[1]
lambda_D_train = lambdas_hat_train$par[2]
lambda_R_train = lambdas_hat_train$par[3]

b_disPEER_train = as.vector(solve(t(Z_P_train) %*% Z_P_train + 
                                    lambda_C_train * Qc_obs +
                                    lambda_D_train * Qd_true +
                                    lambda_R_train * diag(p)) %*% t(Z_P_train) %*% y_P_train )
beta_disPEER_train = as.vector(solve(t(X_train) %*% X_train) %*% t(X_train) %*% (y_train - Z_train %*% b_disPEER_train))
y_disPEER_train = X_train %*% beta_disPEER_train + Z_train %*% b_disPEER_train
y_disPEER_test = X_test %*% beta_disPEER_train + Z_test %*% b_disPEER_train
calc_MSEr(y_disPEER_test, y_test)

### INFORMATION FROM THE CONNECTIONS ONLY ###
### We have here y, Z, Qc and Qd_r and Qd_l ###
riPEER_1_train = riPEER(Q = Qc_obs, 
                        y = as.matrix(y_train),
                        Z = Z_train,
                        X = X_train)
riPEER_1_train$lambda.Q
riPEER_1_train$lambda.R
riPEER_1_train$lambda.2
b_riPEER_train = riPEER_1_train$b.est
beta_riPEER_train = riPEER_1_train$beta.est
y_riPEER_train = X_train %*% beta_riPEER_train[-1] + Z_train %*% b_riPEER_train + beta_riPEER_train[1]
y_riPEER_test = X_test %*% beta_riPEER_train[-1] + Z_test %*% b_riPEER_train + beta_riPEER_train[1]
calc_MSEr(y_riPEER_test, y_test)


estimations_data_y_test = data.frame(y = y_test,
                                     y_ridge = y_ridge_test,
                                     y_riPEER = y_riPEER_test,
                                     y_disPEER = y_disPEER_test,
                                     y_lm = y_lm_test)

ggplot(data = estimations_data_y_test, aes(x = y_test)) +
  geom_line(aes(y = y_test), lty = 2) +
  geom_point(aes(y = y_ridge, colour = "Ridge test estimates"), lwd = 1.7) +
  geom_point(aes(y = y_ridge, colour = "Ridge test estimates"), lwd = 1.7) +
  geom_point(aes(y = y_riPEER, colour = "riPEER test estimates"), lwd = 1.7) +
  geom_point(aes(y = y_disPEER, colour = "disPEER test estimates"), lwd = 1.7) +
  scale_colour_manual("",
                      breaks = c("Ridge test estimates",
                                 "riPEER test estimates", 
                                 "disPEER test estimates"),
                      values = c("#E7B800",
                                 "#00AFBB", 
                                 "red1")) +
  ylab("Estimates") +
  xlab("True values") 

differences_test = data.frame(b_ridge = y_ridge_test - y_test,
                              b_riPEER = y_riPEER_test - y_test,
                              b_disPEER = y_disPEER_test - y_test)

ggplot(data = differences_test) +
  geom_violin(aes(x = "1", y = b_ridge, fill = "Ridge")) +
  geom_violin(aes(x = "2", y = b_riPEER, fill = "riPEER")) + 
  geom_violin(aes(x = "3", y = b_disPEER, fill = "disPEER")) +
  geom_boxplot(aes(x = "1", y = b_ridge), width = 0.2) +
  geom_boxplot(aes(x = "2", y = b_riPEER), width = 0.2) + 
  geom_boxplot(aes(x = "3", y = b_disPEER), width = 0.2) +
  geom_jitter(aes(x = "1", y = b_ridge)) +
  geom_jitter(aes(x = "2", y = b_riPEER)) + 
  geom_jitter(aes(x = "3", y = b_disPEER)) +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(" ",
                    breaks = c("Ridge",
                               "riPEER", 
                               "disPEER"),
                    values = c("#E7B800",
                               "#00AFBB", 
                               "red1")) +
  ylab("Real and estimate values differences") +
  xlab("Method")

y_test



### CI in test ###
### INFORMATION FROM THE CONNECTIONS ONLY ###
### We have here y, Z, Qc and Qd_r and Qd_l ###
riPEER_1_train_ci = riPEER(Q = Qc_obs, 
                           y = as.matrix(y_train),
                           Z = Z_train,
                           X = X_train,
                           compute.boot.CI = TRUE,
                           boot.R = 500)
b_riPEER_train_ci = riPEER_1$b.est
beta_riPEER_train_ci = riPEER_1$beta.est
ConfBand_riPEER_train_ci <- riPEER_1_train_ci$boot.CI

### CI for OLS ###
b.est.boot <- function(data, indices) {
  data.boot = data[indices, ]
  y_boot = as.matrix(data.boot[, 1])
  Z_boot = as.matrix(data.boot[, 2:(p + 1)])
  X_boot = as.matrix(data.boot[, (p + 2):ncol(data)])
  lm_model = lm(y_boot ~ cbind(X_boot, Z_boot))
  b.est.gR.boot = as.vector(coef(lm_model)[(m + 2):(p + m + 1)]) 
  return(as.vector(b.est.gR.boot))
}
boot.out <- boot::boot(data = cbind(y_train, Z_train, X_train), statistic = b.est.boot,
                       R = 500, parallel = "multicore", ncpus = 7)
boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx) {
  boot.ci.out <- (boot::boot.ci(boot.out, type = "bca", index = idx, conf = 0.95))$bca
  return(boot.ci.out[c(4, 5)])
})
ConfBand_OLS_train_ci <- do.call(rbind.data.frame, boot.CI.l)
names(ConfBand_OLS_train_ci) <- c("lower", "upper")

### CI for ridge ###
b.est.boot <- function(data, indices) {
  data.boot = data[indices, ]
  y_boot = as.matrix(data.boot[, 1])
  Z_boot = as.matrix(data.boot[, 2:(p + 1)])
  X_boot = as.matrix(data.boot[, (p + 2):ncol(data)])
  ridge_cv = cv.glmnet(x = cbind(X_boot, Z_boot), 
                       y = y_boot, 
                       alpha = 0, intercept = FALSE)
  b.est.gR.boot = as.vector(coef(ridge_cv)[(m + 2):(p + m + 1)]) 
  return(as.vector(b.est.gR.boot))
}
boot.out <- boot::boot(data = cbind(y_train, Z_train, X_train), statistic = b.est.boot,
                       R = 500, parallel = "multicore", ncpus = 7)
boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx) {
  boot.ci.out <- (boot::boot.ci(boot.out, type = "bca", index = idx, conf = 0.95))$bca
  return(boot.ci.out[c(4, 5)])
})
ConfBand_ridge_train_ci <- do.call(rbind.data.frame, boot.CI.l)
names(ConfBand_ridge_train_ci) <- c("lower", "upper")

### CI for disPEER ###
Blambda_opt = lambda_C * Qc_obs +
  lambda_D * Qd_true +
  lambda_R * diag(p)

b.est.boot <- function(data, indices) {
  data.boot = data[indices, ]
  y_P_boot = as.matrix(data.boot[, 1])
  Z_P_boot = as.matrix(data.boot[, 2:ncol(data)])
  b.est.dis.boot <- as.vector(solve(t(Z_P_boot) %*% Z_P_boot +
                                      Blambda_opt) %*% t(Z_P_boot) %*% y_P_boot)
  return(as.vector(b.est.dis.boot))
}
boot.out <- boot::boot(data = cbind(y_P_train, Z_P_train), statistic = b.est.boot,
                       R = 500, parallel = "multicore", ncpus = 7)
boot.CI.l <- lapply(1:(dim(boot.out$t)[2]), function(idx) {
  boot.ci.out <- (boot::boot.ci(boot.out, type = "bca", index = idx, conf = 0.95))$bca
  return(boot.ci.out[c(4, 5)])
})
ConfBand_disPEER_train_ci <- do.call(rbind.data.frame, boot.CI.l)
names(ConfBand_disPEER_train_ci) <- c("lower", "upper")

Pval_disPEER_train_ci = boot_p(boot.out$t)
Pval_disPEER_train_ci$term = colnames(Qd_true)
Pval_disPEER_train_ci$estimates = b_disPEER_train
ConfBand_disPEER_pval_train_ci = ConfBand_disPEER_train_ci[order(Pval_disPEER_train_ci$p.value, decreasing = FALSE), ]
rownames(ConfBand_disPEER_pval_train_ci) = 1:p
Pval_disPEER_train_ci = Pval_disPEER_train_ci[order(Pval_disPEER_train_ci$p.value, decreasing = FALSE), ]
rownames(Pval_disPEER_train_ci) = 1:p
disPEER_significant_pval_train_ci = (1:p)[Pval_disPEER_train_ci$p.value < 0.04]



ggplot(data = Pval_disPEER_train_ci, 
       aes(x = 1:p)) +
  geom_line(aes(y = estimates), lty = 2) +
  geom_point(aes(y = estimates), lwd = 1.7) +
  geom_ribbon(aes(ymin = ConfBand_disPEER_pval_train_ci$lower, ymax = ConfBand_disPEER_pval_train_ci$upper), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = disPEER_significant_pval_train_ci, col = "red") +
  ylab("disPEER estimations") +
  xlab("cortical regions") +
  ggtitle("Estimates of b coefficients by disPEER with confidence intervals and significance denoted") +
  scale_x_continuous(breaks = 1:p,
                     labels = Pval_disPEER_train_ci$term) +
  theme(axis.text.x = element_text(size = 9, 
                                   angle = 90,
                                   vjust = 0.3)) +
  ylim(-30, 30) 




plot_matrix(Ad_r_true_euc,
            which_mat = " ",
            what_mat = " ",
            which_hemi = " ",
            based_on = " ")

Ac_reordered = plot_matrix(Ac_l_obs[Pval_disPEER_train_ci$term, Pval_disPEER_train_ci$term],
                           which_mat = "Observed",
                           what_mat = "adjacency",
                           which_hemi = "left",
                           based_on = "structural connectivity")

grid.arrange(Ad_reordered, Ac_reordered, ncol = 2)



significant_riPEER_train_ci = (1:p)[!(0 >= ConfBand_riPEER_train_ci$lower & 0 <= ConfBand_riPEER_train_ci$upper)]
significant_disPEER_train_ci = (1:p)[!(0 >= ConfBand_disPEER_train_ci$lower & 0 <= ConfBand_disPEER_train_ci$upper)]
significant_ridge_train_ci = (1:p)[!(0 >= ConfBand_ridge_train_ci$lower & 0 <= ConfBand_ridge_train_ci$upper)]
significant_OLS_train_ci = (1:p)[!(0 >= ConfBand_OLS_train_ci$lower & 0 <= ConfBand_OLS_train_ci$upper)]

estimations_data_b_train_ci = data.table(indexes = 1:p,
                                         b_ridge = b_ridge_train,
                                         b_riPEER = riPEER_1_train_ci$b.est,
                                         b_disPEER = b_disPEER_train,
                                         b_OLS = b_lm,
                                         lower_riPEER = ConfBand_riPEER_train_ci$lower,
                                         upper_riPEER = ConfBand_riPEER_train_ci$upper,
                                         lower_disPEER = ConfBand_disPEER_train_ci$lower,
                                         upper_disPEER = ConfBand_disPEER_train_ci$upper,
                                         lower_ridge = ConfBand_ridge_train_ci$lower,
                                         upper_ridge = ConfBand_ridge_train_ci$upper,
                                         lower_OLS = ConfBand_OLS_train_ci$lower,
                                         upper_OLS = ConfBand_OLS_train_ci$upper)

estimation_plot_riPEER = ggplot(data = estimations_data_b_train_ci, 
                                aes(x = indexes)) +
  geom_line(aes(y = b_riPEER), lty = 2) + 
  geom_point(aes(y = b_riPEER), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_riPEER, ymax = upper_riPEER), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_riPEER, col = "red") +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("riPEER") +
  xlab(" ") +
  ylim(-20, 20) 

estimation_plot_ridge = ggplot(data = estimations_data_b_train_ci, 
                               aes(x = indexes)) +
  geom_line(aes(y = b_ridge), lty = 2) +
  geom_point(aes(y = b_ridge), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_ridge, ymax = upper_ridge), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_ridge, col = "red") +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("ridge") +
  xlab(" ") +
  ylim(-20, 20) 

estimation_plot_OLS = ggplot(data = estimations_data_b_train_ci, 
                             aes(x = indexes)) +
  geom_line(aes(y = b_OLS), lty = 2) +
  geom_point(aes(y = b_OLS), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_OLS, ymax = upper_OLS), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_OLS, col = "red") +
  scale_x_continuous(breaks = 1:p,
                     labels = NULL) +
  ylab("OLS") +
  xlab(" ") +
  ylim(-120, 120) +
  ggtitle(paste0("Estimations of b for the left hemisphere")) 



estimation_plot_disPEER = ggplot(data = estimations_data_b_train_ci, aes(x = indexes)) +
  geom_line(aes(y = b_disPEER), lty = 2) +
  geom_point(aes(y = b_disPEER), lwd = 1.7) +
  geom_ribbon(aes(ymin = lower_disPEER, ymax = upper_disPEER), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2, col = "red") +
  geom_vline(xintercept = significant_disPEER, col = "red") +
  ylab("disPEER") +
  xlab("cortical regions") +
  scale_x_continuous(breaks = 1:p,
                     labels = colnames(Qd_true)) +
  theme(axis.text.x = element_text(size = 9, 
                                   angle = 90,
                                   vjust = 0.3)) +
  ylim(-20, 20) 

grid.arrange(estimation_plot_OLS, 
             estimation_plot_ridge, 
             estimation_plot_riPEER, 
             estimation_plot_disPEER, 
             heights = c(1.2, 1, 1, 2.2))

### VARIATION OF POSSIBLE RESPONSES (y) ###
# ggplot(data = dmgrph_data, aes(x = Age, y = MMSE_Score)) +
#   geom_boxplot() 
# ### practically no dementia observed (one case) ###
# ### Penn Progressive Matrices: number of correct responses ###
# ggplot(data = dmgrph_data, aes(x = Age, y = PMAT24_A_CR)) +
#   ylim(c(0, 25)) +
#   geom_boxplot() 
# ### Pattern Completion Processing Speed ###
# ggplot(data = dmgrph_data, aes(x = Age, y = ProcSpeed_Unadj)) +
#   geom_boxplot() 
# ###
# ggplot(data = dmgrph_data, aes(x = Age, y = ReadEng_Unadj)) +
#   geom_boxplot() 
# ###
# ggplot(data = dmgrph_data, aes(x = Age, y = Flanker_Unadj)) +
#   geom_boxplot() 
# ###
# ggplot(data = dmgrph_data, aes(x = Age, y = PicVocab_Unadj)) +
#   geom_boxplot() 





## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##  

# ### SOLVER PERFORMANCE ###
# ### CHECKING HOM MANY TIMES SOLVER CANNOT SOLVE THE OPTIMIZATION ###
# calc_solver_errors(n = 400, m = 2, 
#                    sigma2b = 0.1,
#                    rC = 0.1, k = 0.01,
#                    p = p_l,
#                    Qc_true = Qc_l_true,
#                    Qc_obs = Qc_l_obs,
#                    Qd_true = Qd_l_true,
#                    Qd_obs = Qd_l_obs,
#                    n_loops = 1000)
# 
# calc_solver_errors(n = 400, m = 2, 
#                    sigma2b = 0.1,
#                    rC = 0.5, k = 0.01,
#                    p = p_l,
#                    Qc_true = Qc_r_true,
#                    Qc_obs = Qc_r_obs,
#                    Qd_true = Qd_r_true,
#                    Qd_obs = Qd_r_obs,
#                    n_loops = 100)
# 
# 

# ### THICKNESS - AREA (original data) ###
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/corrected_AdjMat_results/article_plots")
# 
# p = lapply(1:ncol(Thck_right), function(i) {
#   thck_colnames = colnames(Thck_right)
#   area_colnames = colnames(Area_right)
#   reg_ThckArea = data.table(thck = Thck_right[, thck_colnames[i], 
#                                              with = FALSE],
#                             area = Area_right[, area_colnames[i],
#                                              with = FALSE])
#   reg_ThckArea[, gender := as.factor(dmgrph_data$Gender)]
#   area_name = unlist(strsplit(unlist(strsplit(thck_colnames[i], "FS_R_"))[2], "_Thck"))
#   ThckArea_data = reg_ThckArea
#   ggplot() +
#     geom_point(aes(x = ThckArea_data$thck,
#                    y = ThckArea_data$area,
#                    col = ThckArea_data$gender)) +
#     scale_color_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "thickness measurements")) +
#     ylab(paste(area_name, "area measurements")) +
#     labs(col = "Gender") +
#     ggtitle(paste("Scatterplot of", 
#                   area_name, 
#                   "thickness and area measurements")) +
#     theme_bw()
# })
# ggsave(
#   filename = "scatter_plots_TA_r.pdf", 
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )

# ### THICKNESS - AREA (without outliers) ###
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/corrected_AdjMat_results/article_plots")
# 
# p = lapply(1:ncol(Thck_right), function(i) {
#   thck_colnames = colnames(Thck_right)
#   area_colnames = colnames(Area_right)
#   reg_ThckArea = data.table(thck = Thck_right[, thck_colnames[i], 
#                                               with = FALSE],
#                             area = Area_right[, area_colnames[i],
#                                               with = FALSE])
#   reg_ThckArea[, gender := as.factor(dmgrph_data_r$Gender)]
#   area_name = unlist(strsplit(unlist(strsplit(thck_colnames[i], "FS_R_"))[2], "_Thck"))
#   ThckArea_data = reg_ThckArea
#   ggplot() +
#     geom_point(aes(x = ThckArea_data$thck,
#                    y = ThckArea_data$area,
#                    col = ThckArea_data$gender)) +
#     scale_color_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "thickness measurements")) +
#     ylab(paste(area_name, "area measurements")) +
#     labs(col = "Gender") +
#     ggtitle(paste("Scatterplot of", 
#                   area_name, 
#                   "thickness and area measurements")) +
#     theme_bw()
# })
# ggsave(
#   filename = "scatter_plots_TA_nooutliers_r.pdf", 
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )

# ################################################################
# ### THICKNESS - AREA (without outliers) ###
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/article_plots")
# 
# p = lapply(1:ncol(Thck_left), function(i) {
#   thck_colnames = colnames(Thck_left)
#   area_colnames = colnames(Area_left)
#   reg_ThckArea = data.table(thck = Thck_left[, thck_colnames[i], 
#                                              with = FALSE],
#                             area = Area_left[, area_colnames[i],
#                                              with = FALSE])
#   reg_ThckArea[, gender := as.factor(dmgrph_data$Gender)]
#   area_name = unlist(strsplit(unlist(strsplit(thck_colnames[i], "FS_L_"))[2], "_Thck"))
#   ThckArea_data = reg_ThckArea
#   ggplot() +
#     geom_point(aes(x = ThckArea_data$thck,
#                    y = ThckArea_data$area,
#                    col = ThckArea_data$gender)) +
#     scale_color_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "thickness measurements")) +
#     ylab(paste(area_name, "area measurements")) +
#     labs(col = "Gender") +
#     ggtitle(paste("Scatterplot of", 
#                   area_name, 
#                   "thickness and area measurements")) +
#     theme_bw()
# })
# ggsave(
#   filename = "scatter_plots_TA_nooutliers_l.pdf", 
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )
# 
# ### THICKNESS - AREA (with outliers) ###
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/article_plots")
# 
# p = lapply(1:ncol(Thck_right), function(i) {
#   thck_colnames = colnames(Thck_right)
#   area_colnames = colnames(Area_right)
#   reg_ThckArea = data.table(thck = Thck_right[, thck_colnames[i], 
#                                              with = FALSE],
#                             area = Area_right[, area_colnames[i],
#                                              with = FALSE])
#   reg_ThckArea[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
#   ThckArea_data = reg_ThckArea[-c(outlier_subjects_1, outlier_subjects_2), ]
#   ThckArea_data_out_1 = reg_ThckArea[outlier_subjects_1, ]
#   ThckArea_data_out_2 = reg_ThckArea[outlier_subjects_2, ]
#   area_name = unlist(strsplit(unlist(strsplit(thck_colnames[i], "FS_L_"))[2], "_Thck"))
#   
#   ggplot() +
#     geom_point(aes(x = ThckArea_data$thck,
#                    y = ThckArea_data$area,
#                    col = ThckArea_data$gender)) +
#     geom_point(aes(x = ThckArea_data_out_1$thck,
#                    y = ThckArea_data_out_1$area,
#                    col = ThckArea_data_out_1$gender),
#                shape = 15,
#                size = 4) +
#     geom_point(aes(x = ThckArea_data_out_2$thck,
#                    y = ThckArea_data_out_2$area,
#                    col = ThckArea_data_out_2$gender),
#                shape = 17,
#                size = 4) +
#     scale_color_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "thickness measurements")) +
#     ylab(paste(area_name, "area measurements")) +
#     labs(col = "Gender") +
#     ggtitle(paste("Scatterplot of", 
#                   area_name, 
#                   "thickness and area measurements")) +
#     theme_bw()
# })
# ggsave(
#   filename = "scatter_plots_TA_r.pdf", 
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )
# 
# ### VOLUME - AREA ###
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/article_plots")
# 
# p = lapply(1:ncol(Thck_left), function(i) {
#   vol_colnames = colnames(Vol_left)
#   area_colnames = colnames(Area_left)
#   reg_ThckArea = data.table(vol = Vol_left[, vol_colnames[i],
#                                            with = FALSE],
#                             area = Area_left[, area_colnames[i],
#                                              with = FALSE])
#   reg_ThckArea[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
#   ThckArea_data = reg_ThckArea[-c(outlier_subjects_1, outlier_subjects_2), ]
#   ThckArea_data_out_1 = reg_ThckArea[outlier_subjects_1, ]
#   ThckArea_data_out_2 = reg_ThckArea[outlier_subjects_2, ]
#   area_name = unlist(strsplit(unlist(strsplit(vol_colnames[i], "FS_L_"))[2], "_Vol"))
# 
#   ggplot() +
#     geom_point(aes(x = ThckArea_data$vol,
#                    y = ThckArea_data$area,
#                    col = ThckArea_data$gender)) +
#     geom_point(aes(x = ThckArea_data_out_1$vol,
#                    y = ThckArea_data_out_1$area,
#                    col = ThckArea_data_out_1$gender),
#                shape = 15,
#                size = 4) +
#     geom_point(aes(x = ThckArea_data_out_2$vol,
#                    y = ThckArea_data_out_2$area,
#                    col = ThckArea_data_out_2$gender),
#                shape = 17,
#                size = 4) +
#     scale_color_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "volume measurements")) +
#     ylab(paste(area_name, "area measurements")) +
#     labs(col = "Gender") +
#     ggtitle(paste("Scatterplot of",
#                   area_name,
#                   "volume and area measurements")) +
#     theme_bw()
# })
# ggsave(
#   filename = "scatter_plots_AV_nooutliers.pdf",
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )
# 
# ### VOLUME - THICKNESS ###
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/article_plots")
# 
# p = lapply(1:ncol(Thck_left), function(i) {
#   vol_colnames = colnames(Vol_left)
#   thck_colnames = colnames(Thck_left)
#   reg_ThckArea = data.table(vol = Vol_left[, vol_colnames[i],
#                                            with = FALSE],
#                             thck = Thck_left[, thck_colnames[i],
#                                              with = FALSE])
#   reg_ThckArea[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
#   ThckArea_data = reg_ThckArea[-c(outlier_subjects_1, outlier_subjects_2), ]
#   ThckArea_data_out_1 = reg_ThckArea[outlier_subjects_1, ]
#   ThckArea_data_out_2 = reg_ThckArea[outlier_subjects_2, ]
#   area_name = unlist(strsplit(unlist(strsplit(vol_colnames[i], "FS_L_"))[2], "_Vol"))
# 
#   ggplot() +                           
#     geom_point(aes(x = ThckArea_data$vol,
#                    y = ThckArea_data$thck,
#                    col = ThckArea_data$gender)) +
#     geom_point(aes(x = ThckArea_data_out_1$vol,
#                    y = ThckArea_data_out_1$thck,
#                    col = ThckArea_data_out_1$gender),
#                shape = 15,
#                size = 4) +
#     geom_point(aes(x = ThckArea_data_out_2$vol,
#                    y = ThckArea_data_out_2$thck,
#                    col = ThckArea_data_out_2$gender),
#                shape = 17,
#                size = 4) +
#     scale_color_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "volume measurements")) +
#     ylab(paste(area_name, "thickness measurements")) +
#     labs(col = "Gender") +
#     ggtitle(paste("Scatterplot of",
#                   area_name,
#                   "volume and thickness measurements")) +
#     theme_bw()
# })
# ggsave(
#   filename = "scatter_plots_VT_nooutliers.pdf",
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )
# 
# ############################################################
# ### VOLUME WITHIN GENDER (not normalized) ###
# 
# # Volume_left = Volume_left[-c(outlier_subjects_1, outlier_subjects_2), ]
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/article_plots")
# p = lapply(1:ncol(Volume_left), function(i) {
#   vol_colnames = colnames(Volume_left)
#   reg_ThckArea = data.table(Volume_left[, vol_colnames[i], 
#                                         with = FALSE])
#   reg_ThckArea[, gender := as.factor(dmgrph_data$Gender)]
#   # reg_ThckArea[, gender := ifelse(X[-c(outlier_subjects_1, outlier_subjects_2), ]$gender == 1, 
#   #                                 "female", "male")]
#   colnames(reg_ThckArea) = c("vol", "gender")
#   area_name = unlist(strsplit(unlist(strsplit(vol_colnames[i], "FS_L_"))[2], "_Vol"))
#   
#   ggplot(data = reg_ThckArea) +
#     geom_violin(aes(y = vol,
#                     x = gender,
#                     fill = gender),
#                 lwd = 2, 
#                 trim = FALSE) +
#     geom_jitter(aes(y = vol,
#                     x = gender),
#                 lwd = 2) +
#     scale_fill_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "volume measurements")) +
#     ylab("Volume") +
#     labs(col = "Gender") +
#     ggtitle(paste("Violin plots of", 
#                   area_name, 
#                   "volume measurements for males and females")) +
#     theme_bw()
# })
# ggsave(
#   filename = "NotNormVol_violins_nooutliers.pdf", 
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )
# 
# ### VOLUME WITHIN GENDER (normalized) ###
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/corrected_AdjMat_results/article_plots")
# Volume_left_norm = Volume_left / FS_data$FS_LCort_GM_Vol
# p = lapply(1:ncol(Volume_left_norm), function(i) {
#   vol_colnames = colnames(Volume_left_norm)
#   reg_ThckArea = data.table(Volume_left_norm[, vol_colnames[i], 
#                                              with = FALSE])
#   reg_ThckArea[, gender := as.factor(dmgrph_data$Gender)]
#   # reg_ThckArea[, gender := ifelse(X[-c(outlier_subjects_1, outlier_subjects_2), ]$gender == 1, "female", "male")]
#   colnames(reg_ThckArea) = c("vol", "gender")
#   area_name = unlist(strsplit(unlist(strsplit(vol_colnames[i], "FS_L_"))[2], "_Vol"))
#   
#   ggplot(data = reg_ThckArea) +
#     geom_violin(aes(y = vol,
#                     x = gender,
#                     fill = gender),
#                 lwd = 2, 
#                 trim = FALSE) +
#     geom_jitter(aes(y = vol,
#                     x = gender),
#                 lwd = 2) +
#     scale_fill_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "volume measurements")) +
#     ylab("Volume") +
#     labs(col = "Gender") +
#     ggtitle(paste("Violin plots of", 
#                   area_name, 
#                   "volume measurements for males and females")) +
#     theme_bw()
# })
# ggsave(
#   filename = "NormVolume_violins_nooutliers.pdf", 
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )


# 
# ### REGION VOLUME - TOTAL VOLUME ###
# Vols_left = cbind(Vol_left,
#                   TotVol = FS_data$FS_LCort_GM_Vol)
# 
# setwd("~/Desktop/Scholar_Indiana/my_work/Article/article_plots")
# 
# p = lapply(1:ncol(Thck_left), function(i) {
#   vol_colnames = colnames(Vol_left)
#   reg_ThckArea = data.table(Vols_left[, vol_colnames[i], 
#                                             with = FALSE],
#                             total_vol = Vols_left$TotVol)
#   colnames(reg_ThckArea) = c("vol", "total_vol")
#   reg_ThckArea[, gender := as.factor(ifelse(X$gender == 1, "female", "male"))]
#   ThckArea_data = reg_ThckArea[-c(outlier_subjects_1, outlier_subjects_2), ]
#   ThckArea_data_out_1 = reg_ThckArea[outlier_subjects_1, ]
#   ThckArea_data_out_2 = reg_ThckArea[outlier_subjects_2, ]
#   area_name = unlist(strsplit(unlist(strsplit(vol_colnames[i], "FS_L_"))[2], "_Vol"))
#   
#   ggplot() +
#     geom_point(aes(x = ThckArea_data$vol,
#                    y = ThckArea_data$total_vol,
#                    col = ThckArea_data$gender)) +
#     geom_point(aes(x = ThckArea_data_out_1$vol,
#                    y = ThckArea_data_out_1$total_vol,
#                    col = ThckArea_data_out_1$gender),
#                shape = 15,
#                size = 4) +
#     geom_point(aes(x = ThckArea_data_out_2$vol,
#                    y = ThckArea_data_out_2$total_vol,
#                    col = ThckArea_data_out_2$gender),
#                shape = 17,
#                size = 4) +
#     scale_color_manual(values = c("red", "dodgerblue")) +
#     guides(colour = guide_legend(override.aes = list(shape = 16))) +
#     xlab(paste(area_name, "volume measurements")) +
#     ylab("Left hemisphere total gray matter cortical volume measurements") +
#     labs(col = "Gender") +
#     ggtitle(paste("Scatterplot of", 
#                   area_name, 
#                   "volume and total left hemisphere volume measurements")) +
#     theme_bw()
# })
# ggsave(
#   filename = "scatter_plots_Vols.pdf", 
#   plot = marrangeGrob(p, nrow = 1, ncol = 1)
# )

# VolNorm_left = Vol_left / FS_data$FS_LCort_GM_Vol


# corrplot(cor(ThckArea_left[, 1:33], ThckArea_left[, 34:66]), 
#          method = "shade",
#          title = "Correlation plot between thickness and area for left hemisphere",
#          mar = c(0, 0, 1, 0))
# corrplot(cor(ThckArea_right[, 1:33], ThckArea_right[, 34:66]),
#          method = "shade",
#          title = "Correlation plot between thickness and area for right hemisphere",
#          mar = c(0, 0, 1, 0))

# AvgThck_left = sapply(1:nrow(ThckArea_left), function(i) {
#   weighted.mean(ThckArea_left[i, 1:33], ThckArea_left[i, 34:66])
# })
# 
# AvgThck_right = sapply(1:nrow(ThckArea_right), function(i) {
#   weighted.mean(ThckArea_right[i, 1:33], ThckArea_right[i, 34:66])
# })
# 
# Z_l = Z_l / AvgThck_left
# Z_r = Z_r / AvgThck_right