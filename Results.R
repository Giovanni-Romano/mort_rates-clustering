library(ggplot2)
library(reshape2)
library(dplyr)

load("C:/Users/RomanoGi/Desktop/Bocconi/Ricerca/mort_rates-clustering/1st_sim_data.RData")

theme_set(
  theme_light() +
    theme(text = element_text(size = 14),
      strip.background = element_rect(color = 'gray', fill = 'white'),
      strip.text.x = element_text(color = 'black'),
      strip.text.y = element_text(color = 'black')
    )
)


burn_in <- 1:500





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





# Point estimate of partitions through SALSO
library(salso)
est_clust <- lapply(labels_res, 
                    function(x)
                      apply(x[ , , -burn_in], 2, 
                            function(y)
                              salso(x = t(y),
                                    loss = binder(a = NULL),
                                    maxNClusters = 4,
                                    maxZealousAttempts = 10,
                                    nRuns = 11,
                                    nCores = 11
                                    )
                            )
                    )
str(est_clust)

# Non mettendo "a = NULL" in binder() venivano tutte partizioni con un 
#   solo cluster.
lapply(est_clust, function(x) apply(x, 2, function(y) length(unique(y))))


