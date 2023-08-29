rm(list = ls())
library(ggplot2)
library(reshape2)
library(dplyr)

theme_set(
  theme_light() +
    theme(text = element_text(size = 14),
      strip.background = element_rect(color = 'gray', fill = 'white'),
      strip.text.x = element_text(color = 'black'),
      strip.text.y = element_text(color = 'black')
    )
)


burn_in <- 1:1000





# How many clusters?
agesmax <- apply(S, 2, which.max) - 1
windows()
par(mfrow = c(5, 4))
for (j in 1:20){
  plot(apply(labels_res[[j]][,,-burn_in], 2, function(x) mean(apply(x, 2, function(y) max(y)))), 
       type = "l", 
       ylim = c(1, 2),
       xlab = "time",
       ylab = "avg # of cluster",
       main = paste("Coeff of ~", agesmax[j], " yrs old", sep = ""))
}
# savePlot("plots/avg_#_clusters", type = "jpg")





# Probability of co-clustering for coeff of 1st spline (only one active in 0) 
#   and 12th one (the one active between 46 and 64)
same_cluster0 <- same_cluster55 <- matrix(NA, 6, 87)

up3 <- upper.tri(matrix(1, 4, 4))

same_cluster0 <- sapply(1:T_final, 
                        function(t) rowSums(apply(labels_res[[1]][, t, -burn_in], 2, 
                                                 function(x) outer(x, x, "==")[up3]))
                       )/(n_iter-max(burn_in))
same_cluster55 <- sapply(1:T_final, 
                         function(t) rowSums(apply(labels_res[[12]][, t, -burn_in], 2, 
                                                   function(x) outer(x, x, "==")[up3]))
                         )/(n_iter-max(burn_in))

rownames(same_cluster0) <- rownames(same_cluster55) <- 
  c("It:Sw", "It:UK", "It:US",
    "Sw:UK", "Sw:US", "UK:US")
colnames(same_cluster0) <- colnames(same_cluster55) <- 1933:2019

sc_melt0 <- melt(same_cluster0)
sc_melt55 <- melt(same_cluster55)
colnames(sc_melt0) <- colnames(sc_melt55) <- c("Countries", "Year", "Value")

pl0 <- ggplot(sc_melt0,
              aes(x = Year, y = Value, 
                  group = Countries, 
                  colour = Countries, lty = Countries)) + 
  geom_line(lwd = 1) + 
  ylim(0.55, 1) +
  labs(y = "MonteCarlo prob. of co-clustering", x = "Year",
       title = "0 yrs old")
pl55 <- ggplot(sc_melt55,
               aes(x = Year, y = Value, 
                   group = Countries, 
                   colour = Countries, lty = Countries)) + 
  geom_line(lwd = 1) + 
  ylim(0.55, 1) +
  labs(y = "MonteCarlo prob. of co-clustering", x = "Year",
       title = "55-64 yrs old")
windows(); pl0
windows(); pl55
# ggsave("plots/co-clustering_prob_0.jpg", plot = pl0)
# ggsave("plots/co-clustering_prob_55.jpg", plot = pl55)





# Plot of betas
est_clust_infant <- est_clust[[1]]
colnames(est_clust_infant) <- (1933:2019)
rownames(est_clust_infant) <- names(Y)

beta1_mean <- rowMeans(beta_res[[1]][1, , ])
beta2_mean <- rowMeans(beta_res[[1]][2, , ], na.rm = T)
beta3_mean <- rowMeans(beta_res[[1]][3, , ], na.rm = T)
beta4_mean <- rowMeans(beta_res[[1]][4, , ], na.rm = T)

beta_mean <- cbind(beta1_mean, beta2_mean, beta3_mean, beta4_mean)

beta_it <- beta_mean[cbind(1:87, est_clust_infant[1, ])] 
beta_sw <- beta_mean[cbind(1:87, est_clust_infant[2, ])]
beta_uk <- beta_mean[cbind(1:87, est_clust_infant[3, ])]
beta_us <- beta_mean[cbind(1:87, est_clust_infant[4, ])]

windows()
plot(1933:2019, beta_it, type = 'b', pch = 0,
     xlab = "Year", ylab = "", 
     main = "Values of splines' coeff. for Infant")
lines(1933:2019, beta_sw, type = 'b', col = 2, pch = 1)
lines(1933:2019, beta_uk, type = 'l', col = 7)
lines(1933:2019, beta_us, type = 'b', col = 4, pch = 8)
legend("topright", legend = c("it", "sw", "uk", "us"),
       col = c(1, 2, 7, 4), pch = c(0, 1, 15, 8))

