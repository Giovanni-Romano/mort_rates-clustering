rm(list = ls())
load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/data/simdata_PPMx/simulated_data.RData")
library(ggplot2)
library(reshape2)
library(gridExtra)





### ### ### ### ### ### ### ### ##
##### ggplot2 utils functions ####
### ### ### ### ### ### ### ### ##
# theme_set(
#   theme_light() +
#     theme(text = element_text(size = 10),
#           strip.background = element_rect(color = 'gray', fill = 'white'),
#           strip.text.x = element_text(color = 'black'),
#           strip.text.y = element_text(color = 'black')
#     )
# )

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


### ### ### ### ### ###
#### Plot of labels ####
### ### ### ### ### ###
labels <- out$labs
labels_plots <- list()

for (j in seq_along(labels)){
  lab <- labels[[j]]
  rownames(lab) <- paste("unit", 1:nrow(lab), sep = "")
  colnames(lab) <- paste("t", 1:ncol(lab), sep = "")
  
  lab_melt <- melt(lab,
                   value.name = "cluster",
                   varnames = c("unit", "time"))
  
  lab_melt$cluster <- as.factor(lab_melt$cluster)
  
  tmp.plot <- ggplot(data = lab_melt,
                     mapping = aes(x = time, y = unit, 
                                   pch = cluster, col = cluster)) +
    geom_point(size = 3.5) + 
    scale_color_manual(values = 2:4) +
    ggtitle(paste0("Spline", j, sep = ""))
  labels_plots <- append(labels_plots, list(tmp.plot))
}

# Transform to have just one global legend (and not one per plot)
labels_lgnd <- get_legend(labels_plots[[1]])
for (j in seq_along(labels_plots)){
  labels_plots[[j]] <-labels_plots[[j]] + 
    theme(legend.position="none")
}

# Create grid of plots and save to pdf
pdf("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/sim_study_code/labels.pdf",
    width = 10, height = 5)
grid.arrange(grobs = append(labels_plots, list(labels_lgnd)), 
             ncol = 4, 
             layout_matrix = cbind(matrix(1:6, nrow = 2, byrow = T), 7),
             widths = c(rep(7, 3), 2))
dev.off()




### ### ### ### ### ##
#### Plot of beta ####
### ### ### ### ### ##
beta <- out$beta
beta_plots <- list()

for (b in beta){
  rownames(b) <- paste("cluster", 1:nrow(b), sep = "")
  colnames(b) <- paste("t", 1:ncol(b), sep = "")
  
  b_melt <- melt(b,
                 value.name = "beta",
                 varnames = c("cluster", "time"))
  
  tmp.plot <- ggplot(data = b_melt,
                     mapping = aes(x = time, y = beta, 
                                   group = cluster,
                                   col = cluster)) +
    geom_line()+ 
    scale_color_manual(values = 2:4)
    
  beta_plots <- append(beta_plots, list(tmp.plot))
}

# Transform to have just one global legend (and not one per plot)
beta_lgnd <- get_legend(beta_plots[[1]])
for (j in seq_along(beta_plots)){
  beta_plots[[j]] <-beta_plots[[j]] + 
    theme(legend.position="none")
}

# Create grid of plots and save to pdf
pdf("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/sim_study_code/beta.pdf",
    width = 10, height = 5)
grid.arrange(grobs = append(beta_plots, list(beta_lgnd)), 
             ncol = 4, 
             layout_matrix = cbind(matrix(1:6, nrow = 2, byrow = T), 7),
             widths = c(rep(7, 3), 2))
dev.off()




### ### ### ### ### ### ### 
#### Plot of covariates ####
### ### ### ### ### ### ###
x1 <- out$x1
x2 <- out$x2

rownames(x1) <- rownames(x2) <- paste("unit", 1:out$n, sep = "")
colnames(x1) <- colnames(x2) <- paste("t", 1:ncol(x1), sep = "")

x1_melt <- melt(x1,
                value.name = "x1",
                varnames = c("unit", "time"))
x2_melt <- melt(x2,
                value.name = "x2",
                varnames = c("unit", "time"))


x1.plot <- ggplot(data = x1_melt,
                  mapping = aes(x = time, y = x1, 
                                group = unit,
                                col = unit)) +
  geom_line() + 
  scale_color_manual(values = 5:7)
x2.plot <- ggplot(data = x2_melt,
                  mapping = aes(x = time, y = x2, 
                                group = unit,
                                col = unit)) +
  geom_line()+ 
  scale_color_manual(values = 5:7)

x_lgnd <- get_legend(x1.plot)
x1.plot <- x1.plot + theme(legend.position="none")
x2.plot <- x2.plot + theme(legend.position="none")

# Create grid of plots and save to pdf
pdf("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/sim_study_code/covariates.pdf",
    width = 10, height = 5)
grid.arrange(grobs = list(x1.plot, x2.plot, x_lgnd), 
             ncol = 2, 
             layout_matrix = cbind(matrix(1:2, nrow = 2, byrow = F), 3),
             widths = c(5, 0.5))
dev.off()
