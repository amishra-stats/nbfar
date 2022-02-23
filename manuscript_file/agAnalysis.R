
###########################################################
## R code for simulation examples and applications       ##
##     Negative binomial factor regression model         ##
##    with application to microbiome data anaysis        ##
###########################################################



##########################################################
########     American Gut Data Analysis     ##############
##########################################################

rm(list = ls())

# Load relevant package
library(nbfar)
library(MASS)
require(ggplot2)
require(reshape2)
library(magrittr)
library(mpath)
library(RcppParallel)
library(tictoc)
library(phyloseq)




# set seed
SD = 1234; set.seed(SD)


## Generate predictor matrix:
## load("final_ag_vio.rda")
data_url = "https://github.com/amishra-stats/nbfar/raw/master/manuscript_file/final_ag_vio.rda"
download.file(data_url, destfile='temp.rds')
load("temp.rds")

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 59 taxa and 626 samples ]
# sample_data() Sample Data:       [ 626 samples by 360 sample variables ]
# tax_table()   Taxonomy Table:    [ 59 taxa by 6 taxonomic ranks ]




X <- sample_data(agdata_vio) %>% as.matrix() %>%
  apply(2,as.numeric)
# X <- X %>% scale();
X <- X[,apply(X, 2, sd) > 0];
num_index <- apply(X,2, function(x)length(unique(x))) > 2
X[,num_index] <- (X[,num_index] %>% scale(center = T, scale = T))/2
Y <- otu_table(agdata_vio)@.Data %>% t()
Y <- Y + 1
# preprocess Y and add a pseudo count for the zeros
n <- nrow(X); q = ncol(Y)
# Offset term is log sum of the abundance
offset <- matrix(rowMeans(log(Y)),n,q)


cIndex <- c('sex', 'bmi', "age_years")
cIndex <- match(cIndex, colnames(X))
X <- cbind(X[,cIndex],X[,-cIndex]); cIndex <- 1:length(cIndex)
## fit the model with or without control variable [medical condition] in the model
colnames(X)[cIndex]
example_seed = SD




## ----- ----- ----- ----- ----- Full data model ----- ----- ----- -----
## Define relevant output files
fname <- 'agdata-nbfar2.rda'  # out file name

# NBFAR: Model with offset
tic()
nthread = 1; rank.est = 10
RcppParallel::setThreadOptions(numThreads = nthread)
set.seed(example_seed)
nlam = 40; eps = 1e-5; sp = 0.5
control_nbfar <- nbfar_control(gamma0 = 1, spU = sp, spV = 0.95, maxit = 3000,
                               lamMinFac = 1e-8, lamMaxFac = 1, epsilon = eps, objI = 0,
                               initmaxit = 10000, initepsilon = 1e-8)
nbfar_test_of <- nbfar(Y, X, maxrank = rank.est, nlambda = nlam, cIndex = cIndex,
                       ofset = offset, control = control_nbfar, nfold = 10,
                       PATH = F,nthread = nthread, trace = F)
nbfar_test_of$time = toc()


# NBRRR
tic()
set.seed(example_seed)
rank.est = 5
control_r3 <- nbfar_control(initmaxit = 10000,
                            initepsilon = 1e-6, objI = 1)
nbrrr_test_of <- nbrrr(Y, X, maxrank = rank.est, ofset = offset,
                       control = control_r3, nfold = 10,
                       cIndex = cIndex,trace = F)
nbrrr_test_of$time = toc()
save(list = ls(),file = fname)


## We have saved the output file
data_url = "https://github.com/amishra-stats/nbfar/raw/master/manuscript_file/agdata-nbfar2.rda"
download.file(data_url, destfile='temp.rds')
load("temp.rds")




## -------------------------------------------
## Inference based on NB-FAR output

## Download annotation of the covariates from GitHub
data_url = "https://github.com/amishra-stats/nbfar/raw/master/manuscript_file/agdata_anot.rda"
download.file(data_url, destfile='temp.rds')
load("temp.rds")



## Process V estimate from the NB-FAR output
library(phyloseq)
Phylum <- tax_table(agdata_vio)@.Data[,'Phylum']
species_effect = data.frame(Phylum, Mean_abundance = colMeans(Y), nbfar_test_of$V)
rownames(species_effect) = colnames(Y)
species_effect = species_effect[with(species_effect, order(Mean_abundance,decreasing = T)),]
#species_effect
Vest <- as.matrix(species_effect[,c(1,3,4,5)])
cutv <- quantile(abs(as.numeric(as.vector(Vest[,2:4]))),0.5 )
temp <- apply(Vest[,2:4],2,as.numeric)
Vest[,2:4] <- temp <- temp*(abs(temp)>cutv)
temp <- sign(temp)
temp <- data.frame(Phylum = Vest[,'Phylum'], temp)
stoOrd <- with(temp, order(Phylum,X1,X2))
Vest <- Vest[stoOrd, ]
species_effect <- species_effect[stoOrd,]
v_diff <- which(diff(as.numeric(species_effect$Phylum))>0)
v_diff  ## Index of chnage in the Phylum level
## Number of non-zero entries in columns of estimated V: [11,9 ,23]
colSums(Vest!=0)




## Process U estimate from the NB-FAR output
covariate_effect = data.frame(nbfar_test_of$U)
rownames(covariate_effect) = colnames(X)[-cIndex]
covariate_effect$Type <- covarname[rownames(covariate_effect),'Category']
Uest <- as.matrix(covariate_effect)
cutu <- quantile(abs(as.numeric(as.vector(Uest[,1:3]))),0.5 )

temp <- apply(Uest[,1:3],2,as.numeric)
Uest[,1:3] <- temp <- temp*(abs(temp)>cutu)
temp <- sign(temp)
temp <- data.frame(Type = Uest[,'Type'], temp)
stoOrd <- with(temp, order(Type,X1,X2,X3))
Uest <- Uest[stoOrd, ]
covariate_effect <- covariate_effect[stoOrd,]
u_diff <- which(diff(as.numeric(as.factor(Uest[,'Type'])))>0)
## Number of non-zero entries in columns of estimated U: [57, 61, 72 ]
colSums(covariate_effect[,1:3]!=0)



# -----------------------------------------------------
## Unit rank components of the coefficient matrix
ind <- nbfar_test_of$D > 0
Uest_nbfar <- apply(as.matrix(Uest[,1:3]),2,as.numeric)
Vest_nbfar <- apply(as.matrix(Vest[,2:4]),2,as.numeric)
Cvar <- Uest_nbfar%*%(nbfar_test_of$D*t(Vest_nbfar))
plC <- vector("list",sum(ind) + 1)    # list of full coefficient matrix and its unit rank components
mxcvar <- max(abs(Cvar))
plC[[1]] <-  Cvar/mxcvar     # full estimate of the low-rank coefficient matrix
for (i in 1:sum(ind)) {
  # estimate of the unit-rank component of the low-rank coefficient matrix
  plC[[i + 1]] <- (nbfar_test_of$D[i]*as.matrix(Uest_nbfar[,i]) %*% t(Vest_nbfar[,i]))/mxcvar
}
pmat <- plC[[1]]
ncl <- 500
brk <- seq(-1,1,length.out = ncl)
tem1 <- max(abs(pmat))
color_palette <- colorRampPalette(c( "blue","grey98","red"))(length(brk) - 1)




# load packages
library(gridExtra)
library(tidyverse)
library(dplyr)
library(repr)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
suppressPackageStartupMessages(library(dendextend))
options(repr.plot.width=20, repr.plot.height=24)


# -----------------------------------------------------
# Figure 3 in the main manuscript


unique(species_effect$Phylum)
options(repr.plot.width=30, repr.plot.height=10)
col_fun = colorRamp2(c(-1, 0, 1), c( "blue","white","red"))
Uname <- rownames(covariate_effect)
Vname <- rownames(species_effect)
i <- 1
pmat <-  1.5*plC[[i]]
r_sel <- rowSums(pmat)!=0
c_sel <- colSums(pmat)!=0
u_diff_sub <- which(diff(as.numeric(as.factor(covariate_effect$Type[r_sel])))>0)
v_diff_sub <- which(diff(as.numeric(as.factor(species_effect$Phylum[c_sel])))>0)
#
Phylum <- unique(species_effect$Phylum[c_sel])
Features <- unique(covariate_effect$Type[r_sel])
annot_col <- list(Phylum = brewer.pal(length(Phylum), "Set1"),
                  Features = brewer.pal(length(Features), "Dark2"))
names(annot_col$Phylum) <- Phylum
names(annot_col$Features) <- Features
## -------------------
i <- 1
pmat <-  1.5*plC[[i]]
pmat <- pmat[r_sel,c_sel]
pmat2 <- pmat
rownames(pmat2) <- Uname[r_sel]
colnames(pmat2) <- Vname[c_sel]
pmat2 <- t(pmat2)


i <- 2
pmat3 <-  1.5*plC[[i]][r_sel,c_sel]
rownames(pmat3) <- Uname[r_sel]
colnames(pmat3) <- Vname[c_sel]
pmat3 <- t(pmat3)

i <- 3
pmat4 <-  1.5*plC[[i]][r_sel,c_sel]
rownames(pmat4) <- Uname[r_sel]
colnames(pmat4) <- Vname[c_sel]
pmat4 <- t(pmat4)

i <- 4
pmat5 <-  1.5*plC[[i]][r_sel,c_sel]
rownames(pmat5) <- Uname[r_sel]
colnames(pmat5) <- Vname[c_sel]
pmat5 <- t(pmat5)

## -------------------
spanot <- data.frame(Phylum = species_effect$Phylum[c_sel])
covaranot <- data.frame(Features = covariate_effect$Type[r_sel])
rownames(covaranot) <- Uname[r_sel]
rownames(spanot) <- Vname[c_sel]
ha_col = HeatmapAnnotation(Features = covaranot$Features, col = annot_col,
                           show_annotation_name = F,
                           show_legend = FALSE)
ha_row = rowAnnotation(df = spanot, col = annot_col, show_annotation_name = F,
                       show_legend = FALSE)
ht1 = Heatmap(pmat2, cluster_rows = FALSE, cluster_columns =  FALSE, name = "C",
              col = col_fun, border = TRUE,
              column_title = "C", column_title_gp = gpar(fontsize = 14, fontface = 'bold'),
              show_row_names = F, show_column_names = FALSE, show_heatmap_legend = T,
              row_names_side = 'right', row_names_gp = gpar(fontsize = 12,fontface = 'bold'),
              top_annotation = ha_col, left_annotation = ha_row,
              #               rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          legend_height = unit(15, "cm"),
                                          labels_gp = gpar(fontsize = 18),
                                          at = c(-1, 0,  1),
                                          title_position = 'topcenter',
                                          title = "")
)
ht2 = Heatmap(pmat3, cluster_rows = FALSE, cluster_columns =  FALSE,  name = "C1",
              col = col_fun,border = TRUE,
              column_title = expression(bold('C'[1])),
              column_title_gp = gpar(fontsize = 14, fontface = 'bold'),
              show_row_names = F, show_column_names = FALSE, show_heatmap_legend = F,
              row_names_side = 'right', row_names_gp = gpar(fontsize = 12,fontface = 'bold'),
              top_annotation = ha_col, #left_annotation = ha_row,
              #               rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          at = c(-1, -0.5,0, 0.5, 1),
                                          title_position = 'topcenter',
                                          title = "")
)
ht3 = Heatmap(pmat4, cluster_rows = FALSE, cluster_columns =  FALSE,  name = "C2",
              column_title = expression(bold('C'[2])),
              col = col_fun,border = TRUE,
              column_title_gp = gpar(fontsize = 14, fontface = 'bold'),
              show_row_names = F, show_column_names = FALSE, show_heatmap_legend = F,
              row_names_side = 'right', row_names_gp = gpar(fontsize = 12,fontface = 'bold'),
              top_annotation = ha_col, #left_annotation = ha_row,
              #               rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          at = c(-1, -0.5,0, 0.5, 1),
                                          title_position = 'topcenter',
                                          title = "")
)
ht4 = Heatmap(pmat5, cluster_rows = FALSE, cluster_columns =  FALSE,name = "C3",
              column_title = expression(bold('C'[3])),
              col = col_fun,border = TRUE,
              column_title_gp = gpar(fontsize = 14, fontface = 'bold'),
              show_row_names = T, show_column_names = FALSE, show_heatmap_legend = F,
              row_names_side = 'right', row_names_gp = gpar(fontsize = 14,fontface = 'bold'),
              top_annotation = ha_col, #left_annotation = ha_row,
              #               rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          at = c(-1, 0, 1),
                                          title_position = 'topcenter',
                                          title = "")
)


col_fun = colorRamp2(c(-1, 0, 1), c( "blue","grey98","red"))
lgdp = Legend(labels =  Phylum, title = "Phylum[R]", direction = "horizontal", nrow = 2,
              title_position = 'topcenter', gap = unit(50, "mm"),
              grid_height = unit(8, "mm"),
              title_gp = gpar(fontsize = 20, fontface = "bold",col = 'blue'),
              labels_gp = gpar(fontsize = 20, fontface = 'bold'),
              title_gap = unit(3, "mm"),
              legend_gp = gpar(fill = brewer.pal(length(Phylum), "Set1") ))
lgdc = Legend(labels =  Features, title = "Features[C]", direction = "horizontal",
              title_position = 'topcenter', gap = unit(5, "mm"),   nrow = 2,
              title_gp = gpar(fontsize = 20, fontface = "bold",col = 'blue'),
              title_gap = unit(3, "mm"), grid_height = unit(8, "mm"),
              labels_gp = gpar(fontsize = 20, fontface = 'bold'),
              legend_gp = gpar(fill = brewer.pal(length(Phylum), "Dark2")) )


pd = packLegend(lgdp, lgdc, direction = "horizontal", column_gap = unit(70, "mm"))

# setEPS()
# postscript(file = file.path(figfol, 'agdata_fig1.eps'), height = 8 , width = 20 )

draw(ht1 + ht2+ ht3+ ht4, heatmap_legend_side = "left",  annotation_legend_side = "top", gap= unit(5, "mm"),
     annotation_legend_list = pd)
decorate_heatmap_body("C", {
  for(i in u_diff_sub)
    grid.lines(c(i,i)/sum(r_sel),c(0, 1), gp = gpar(col = "black",lty = 2,lwd = 2))
  for(i in v_diff_sub)
    grid.lines(c(0, 1), 1-c(i,i)/sum(c_sel), gp = gpar(col = "black",lty = 2,lwd = 2))

})
decorate_heatmap_body("C1", {
  for(i in u_diff_sub)
    grid.lines(c(i,i)/sum(r_sel),c(0, 1), gp = gpar(col = "black",lty = 2,lwd = 2))
  for(i in v_diff_sub)
    grid.lines(c(0, 1), 1-c(i,i)/sum(c_sel), gp = gpar(col = "black",lty = 2,lwd = 2))

})
decorate_heatmap_body("C2", {
  for(i in u_diff_sub)
    grid.lines(c(i,i)/sum(r_sel),c(0, 1), gp = gpar(col = "black",lty = 2,lwd = 2))
  for(i in v_diff_sub)
    grid.lines(c(0, 1), 1-c(i,i)/sum(c_sel), gp = gpar(col = "black",lty = 2,lwd = 2))

})
decorate_heatmap_body("C3", {
  for(i in u_diff_sub)
    grid.lines(c(i,i)/sum(r_sel),c(0, 1), gp = gpar(col = "black",lty = 2,lwd = 2))
  for(i in v_diff_sub)
    grid.lines(c(0, 1), 1-c(i,i)/sum(c_sel), gp = gpar(col = "black",lty = 2,lwd = 2))

})

# dev.off()







# -----------------------------------------------------
# Figure 4 in the main manuscript



library(circlize)
col_fun = colorRamp2(c(-1, 0, 1), c( "blue","white","red"))
unique(species_effect$Phylum)
options(repr.plot.width=30, repr.plot.height=11)
Uname <- rownames(covariate_effect)
Vname <- rownames(species_effect)
i <- 1
pmat <-  plC[[i]]
rownames(pmat) <- Uname
colnames(pmat) <- Vname
pmat <- t(pmat)

spanot <- data.frame(Phylum = species_effect$Phylum)
covaranot <- data.frame(Features = covariate_effect$Type)
rownames(covaranot) <- Uname
rownames(spanot) <- Vname


Phylum <- unique(species_effect$Phylum)
Features <- unique(covariate_effect$Type)
annot_col <- list(Phylum = brewer.pal(length(Phylum), "Set1"),
                  Features = brewer.pal(length(Features), "Dark2"))
names(annot_col$Phylum) <- Phylum
names(annot_col$Features) <- Features

# covariate selection
ord_index <- order(colSums(abs(pmat)), decreasing = T)
# colSums(abs(pmat))[ord_index]
pmat <- pmat[,ord_index[1:25]]
pmat <- pmat[rowSums(abs(pmat))!=0,]
sp_sel <- rownames(pmat)
covar_sel <- colnames(pmat)



ha_col = HeatmapAnnotation(Features = covaranot[covar_sel,'Features'], col = annot_col,
                           show_annotation_name = F,  show_legend = F)
ha_row = rowAnnotation(Phylum = spanot[sp_sel,'Phylum'], col = annot_col, show_annotation_name = F,
                       show_legend = F)
ht1 = Heatmap(pmat, cluster_rows = FALSE, cluster_columns =  FALSE, name = "C1",
              col = col_fun,
              #               column_title = "C", column_title_gp = gpar(fontsize = 20, fontface = 'bold'),
              show_row_names = F, show_column_names = T, show_heatmap_legend = T,
              row_names_side = 'left', row_names_gp = gpar(fontsize = 14,fontface = 'bold'),
              column_names_gp = gpar(fontsize = 16,fontface = 'bold'),
              top_annotation = ha_col, left_annotation = ha_row,
              border = TRUE, column_names_max_height = unit(10, "cm"),
              rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          legend_height = unit(10, "cm"),
                                          labels_gp = gpar(fontsize = 18),
                                          at = c(-1, 0,  1),
                                          title_position = 'topcenter',
                                          title = "")
)


# covariate selection
pmat2 <-  plC[[1]]
rownames(pmat2) <- Uname
colnames(pmat2) <- Vname
pmat2 <- t(pmat2)
sel_index_sp <- apply(pmat2, 1,function(x) {
  order(abs(x), decreasing = T)[1:7]
})
cov_index_sel <- sort(unique(as.vector(sel_index_sp)))
pmat2 <- pmat2[,cov_index_sel]
sp_sel
pmat2 <- pmat2[sp_sel,]
sp_sel <- rownames(pmat2)
covar_sel <- colnames(pmat2)

ha_col2 = HeatmapAnnotation(Features = covaranot[covar_sel,'Features'], col = annot_col,
                            show_annotation_name = F,  show_legend = F)
ha_row2 = rowAnnotation(Phylum = spanot[sp_sel,'Phylum'], col = annot_col, show_annotation_name = F,
                        show_legend = F)
ht2 = Heatmap(pmat2, cluster_rows = FALSE, cluster_columns =  FALSE, name = "C2",
              col = col_fun,
              #               column_title = "C", column_title_gp = gpar(fontsize = 20, fontface = 'bold'),
              show_row_names = T, show_column_names = T, show_heatmap_legend = F,
              row_names_side = 'right', row_names_gp = gpar(fontsize = 14,fontface = 'bold'),
              column_names_gp = gpar(fontsize = 16,fontface = 'bold'),
              top_annotation = ha_col2, left_annotation = ha_row2,
              border = TRUE, column_names_max_height = unit(10, "cm"),
              rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          legend_height = unit(10, "cm"),
                                          labels_gp = gpar(fontsize = 18),
                                          at = c(-1, 0,  1),
                                          title_position = 'topcenter',
                                          title = "")
)



lgdp = Legend(labels =  Phylum, title = "Phylum[R]", direction = "horizontal", nrow = 2,
              title_position = 'topcenter', gap = unit(50, "mm"),
              grid_height = unit(8, "mm"),
              title_gp = gpar(fontsize = 20, fontface = "bold",col = 'blue'),
              labels_gp = gpar(fontsize = 20, fontface = 'bold'),
              title_gap = unit(5, "mm"),
              legend_gp = gpar(fill = brewer.pal(length(Phylum), "Set1") ))
lgdc = Legend(labels =  Features, title = "Features[C]", direction = "horizontal",
              title_position = 'topcenter', gap = unit(3, "mm"),   nrow = 2,
              title_gp = gpar(fontsize = 20, fontface = "bold",col = 'blue'),
              title_gap = unit(5, "mm"), grid_height = unit(8, "mm"),
              labels_gp = gpar(fontsize = 20, fontface = 'bold'),
              legend_gp = gpar(fill = brewer.pal(length(Phylum), "Dark2")) )


pd = packLegend(lgdp, lgdc, direction = "horizontal", column_gap = unit(90, "mm"))

# setEPS()
# postscript(file = file.path(figfol, 'agdata_fig2.eps'), height = 11 , width = 20 )

draw(ht1 + ht2, heatmap_legend_side = "left",  annotation_legend_side = "top", gap= unit(10, "mm"),
     annotation_legend_list = pd)
# dev.off()





# -----------------------------------------------------
# Figure 5 in the main manuscript



unique(species_effect$Phylum)
Uname <- rownames(covariate_effect)
Vname <- rownames(species_effect)

## Color anotdation
spanot <- data.frame(Phylum = species_effect$Phylum)
covaranot <- data.frame(Features = covariate_effect$Type)
rownames(covaranot) <- Uname
rownames(spanot) <- Vname
Phylum <- unique(species_effect$Phylum)
Features <- unique(covariate_effect$Type)
annot_col <- list(Phylum = brewer.pal(length(Phylum), "Set1"),
                  Features = brewer.pal(length(Features), "Dark2"))
names(annot_col$Phylum) <- Phylum
names(annot_col$Features) <- Features


## Prepare important covariates matrix
pmat2 <-  plC[[1]]
range(pmat2)
rownames(pmat2) <- Uname
colnames(pmat2) <- Vname
imp_cov_mat <- apply(pmat2, 2,function(x) {
  out <- rep(0,length(x))
  out[order(abs(x), decreasing = T)[1:7]] <- 1
  out
})
imp_cov_mat <- t(imp_cov_mat)
sum(imp_cov_mat)/39 # sanity check

imp_entries = function(j, i, x, y, width, height, fill) {
  if(imp_covar[i, j] == 1){
    grid.points(x, y, pch = 8)
  }
}
set.seed(1)
dist_method = "euclidean" # maximum, binary
clust_method = "average"
anot_name <- c()
var_name_T <- c()




## ************************************ Uncomment it to save it in file  -----
# setEPS()
# postscript(file = file.path(figfol, 'agdata_fig3.eps'), height = 16 , width = 38 )

pd = packLegend(lgdp, lgdc, direction = "vertical", column_gap = unit(90, "mm"))
grid.newpage()

## Draw Heatmap for the selected subset of the covariates
i <- 2
pmat <-  plC[[i]]
r_sel <- order(rowSums(abs(pmat)), decreasing = 2)[1:20]
c_sel <- colSums(abs(pmat))> 0
pmat <- pmat[r_sel,c_sel]
rownames(pmat) <- Uname[r_sel]
colnames(pmat) <- Vname[c_sel]
pmat <- t(pmat)
imp_covar <- imp_cov_mat[c_sel,r_sel]
# Plot
sp_sel <- rownames(pmat)
covar_sel <- colnames(pmat)
ha_col = HeatmapAnnotation(Features = covaranot[covar_sel,'Features'], col = annot_col,
                           show_annotation_name = F,  show_legend = F)
ha_row = rowAnnotation(Phylum = spanot[sp_sel,'Phylum'], col = annot_col, show_annotation_name = F,
                       show_legend = F)
colnames(pmat) <- paste('A',1:ncol(pmat),sep = '')
anot_name <- c(anot_name, colnames(pmat)[colSums(imp_covar)>0])
var_name_T <- c(var_name_T,  covar_sel[colSums(imp_covar)>0])
ht1 = Heatmap(pmat, cluster_rows = T, cluster_columns =  T, name = "C1",
              col = col_fun, column_km = 2, row_km = 2,
              cell_fun = imp_entries,
              clustering_distance_rows = dist_method,
              clustering_method_rows = clust_method,
              clustering_distance_columns = dist_method,
              clustering_method_columns = clust_method,
              column_title = expression(bold('C')[1]),
              show_row_dend = F, show_column_dend = F,
              row_names_max_width = unit(9, "cm"),
              column_title_gp = gpar(fontsize = 30, fontface = 'bold'),
              show_row_names = T, show_column_names = T, show_heatmap_legend = F,
              row_names_side = 'left', row_names_gp = gpar(fontsize = 25,fontface = 'bold'),
              column_names_gp = gpar(fontsize = 30,fontface = 'bold'),
              top_annotation = ha_col, left_annotation = ha_row,
              border = TRUE, column_names_max_height = unit(10, "cm"),
              rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "horizontal",
                                          col_fun = col_fun,
                                          legend_height = unit(10, "cm"),
                                          labels_gp = gpar(fontsize = 18),
                                          at = c(-1, 0,  1),
                                          title_position = 'topcenter',
                                          title = "")
)



options(repr.plot.width=20, repr.plot.height=10)

pushViewport(viewport(layout = grid.layout(nr = 2, nc = 2, heights = c(1,1.6))))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(ht1, newpage = FALSE)
upViewport()


i <- 3
pmat <-  plC[[i]]
r_sel <- order(rowSums(abs(pmat)), decreasing = 2)[1:20]
c_sel <- colSums(abs(pmat))> 0
pmat <- pmat[r_sel,c_sel]
rownames(pmat) <- Uname[r_sel]
colnames(pmat) <- Vname[c_sel]
pmat <- t(pmat)
imp_covar <- imp_cov_mat[c_sel,r_sel]
# Plot
sp_sel <- rownames(pmat)
covar_sel <- colnames(pmat)
ha_col = HeatmapAnnotation(Features = covaranot[covar_sel,'Features'], col = annot_col,
                           show_annotation_name = F,  show_legend = F)
ha_row = rowAnnotation(Phylum = spanot[sp_sel,'Phylum'], col = annot_col, show_annotation_name = F,
                       show_legend = F)
colnames(pmat) <- paste('B',1:ncol(pmat),sep = '')
anot_name <- c(anot_name, colnames(pmat)[colSums(imp_covar)>0])
var_name_T <- c(var_name_T,  covar_sel[colSums(imp_covar)>0])
ht2 = Heatmap(pmat, cluster_rows = T, cluster_columns =  T, name = "C2",
              col = col_fun, column_km = 2, row_km = 2,
              cell_fun = imp_entries,
              clustering_distance_rows = dist_method,
              clustering_method_rows = clust_method,
              clustering_distance_columns = dist_method,
              clustering_method_columns = clust_method,
              column_title = expression(bold('C')[2]),
              show_row_dend = F, show_column_dend = F,
              row_names_max_width = unit(9, "cm"),
              column_title_gp = gpar(fontsize = 30, fontface = 'bold'),
              show_row_names = T, show_column_names = T, show_heatmap_legend = F,
              row_names_side = 'left', row_names_gp = gpar(fontsize = 25,fontface = 'bold'),
              column_names_gp = gpar(fontsize = 30,fontface = 'bold'),
              top_annotation = ha_col, left_annotation = ha_row,
              border = TRUE, column_names_max_height = unit(10, "cm"),
              rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          legend_height = unit(10, "cm"),
                                          labels_gp = gpar(fontsize = 18),
                                          at = c(-1, 0,  1),
                                          title_position = 'topcenter',
                                          title = "")
)


pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(ht2, newpage = FALSE)
upViewport()


i <- 4
pmat <-  plC[[i]]
r_sel <- order(rowSums(abs(pmat)), decreasing = 2)[1:20]
c_sel <- colSums(abs(pmat))> 0.1
pmat <- pmat[r_sel,c_sel]
rownames(pmat) <- Uname[r_sel]
colnames(pmat) <- Vname[c_sel]
pmat <- t(pmat)
imp_covar <- imp_cov_mat[c_sel,r_sel]
# Plot
sp_sel <- rownames(pmat)
covar_sel <- colnames(pmat)
ha_col = HeatmapAnnotation(Features = covaranot[covar_sel,'Features'], col = annot_col,
                           show_annotation_name = F,  show_legend = F)
ha_row = rowAnnotation(Phylum = spanot[sp_sel,'Phylum'], col = annot_col, show_annotation_name = F,
                       show_legend = F)
colnames(pmat) <- paste('C',1:ncol(pmat),sep = '')
anot_name <- c(anot_name, colnames(pmat)[colSums(imp_covar)>0])
var_name_T <- c(var_name_T,  covar_sel[colSums(imp_covar)>0])
ht3 = Heatmap(pmat, cluster_rows = T, cluster_columns =  T, name = "C3",
              col = col_fun, column_km = 2, row_km = 2,
              cell_fun = imp_entries,
              clustering_distance_rows = dist_method,
              clustering_method_rows = clust_method,
              clustering_distance_columns = dist_method,
              clustering_method_columns = clust_method,
              column_title = expression(bold('C')[3]),
              show_row_dend = F, show_column_dend = F,
              row_names_max_width = unit(9, "cm"),
              column_title_gp = gpar(fontsize = 30, fontface = 'bold'),
              show_row_names = T, show_column_names = T, show_heatmap_legend = F,
              row_names_side = 'left', row_names_gp = gpar(fontsize = 25,fontface = 'bold'),
              column_names_gp = gpar(fontsize = 30,fontface = 'bold'),
              top_annotation = ha_col, left_annotation = ha_row,
              border = TRUE, column_names_max_height = unit(11, "cm"),
              rect_gp = gpar(col = "white", lwd = 1),
              heatmap_legend_param = list(direction = "vertical",
                                          col_fun = col_fun,
                                          legend_height = unit(10, "cm"),
                                          labels_gp = gpar(fontsize = 18),
                                          at = c(-1, 0,  1),
                                          title_position = 'topcenter',
                                          title = "")
)

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(ht3, newpage = FALSE)
upViewport()


lgdp = Legend(labels =  Phylum, title = "Phylum[R]", direction = "vertical", ncol = 1,
              title_position = 'topcenter', gap = unit(50, "mm"),
              grid_height = unit(8, "mm"),
              title_gp = gpar(fontsize = 25, fontface = "bold",col = 'blue'),
              labels_gp = gpar(fontsize = 25, fontface = 'bold'),
              title_gap = unit(5, "mm"),
              legend_gp = gpar(fill = brewer.pal(length(Phylum), "Set1") ))
lgdc = Legend(labels =  Features, title = "Features[C]", direction = "vertical",
              title_position = 'topcenter', gap = unit(3, "mm"),   ncol = 1,
              title_gp = gpar(fontsize = 25, fontface = "bold",col = 'blue'),
              title_gap = unit(5, "mm"), grid_height = unit(8, "mm"),
              labels_gp = gpar(fontsize = 25, fontface = 'bold'),
              legend_gp = gpar(fill = brewer.pal(length(Phylum), "Dark2")) )
lgds = Legend(col_fun = col_fun, title = "", at = c(-1, 0,  1), legend_height = unit(5, "cm"),
              grid_width = unit(1, "cm"), labels_gp = gpar(fontface = "bold", fontsize = 15))


pd = packLegend(lgds, lgdp, lgdc, direction = "vertical", row_gap = unit(5, "mm"))
labelname = paste(anot_name, var_name_T, sep = ': ')
lgd = Legend(labels = labelname, legend_gp = gpar(fill = 'grey'), title = "Dictionary",
             title_position = "topcenter", ncol = 2, title_gap = unit(7, "mm"),
             gap = unit(50, "mm"),
             grid_height = unit(8, "mm"), grid_width = unit(6, "mm"),
             title_gp = gpar(fontsize = 25, fontface = "bold",col = 'blue'),
             labels_gp = gpar( fontsize = 25,fontface = 'bold'))


pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(pd, x = unit(0.1, "npc"), y = unit(0.5, "npc"))
draw(lgd, x = unit(0.6, "npc"), y = unit(0.5, "npc"))
upViewport()
upViewport()

# dev.off()

