rm(list = ls())
library(tidyverse)
library(reshape2)
library(patchwork)
library(ggrepel)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
load("14_countries_preproc.RData")


theme_set(
  theme_light() +
    theme(strip.background = element_rect(color = 'gray', fill = 'white'),
          strip.text.x = element_text(color = 'black'),
          strip.text.y = element_text(color = 'black')
    )
)


Y_m <- lapply(rates_male, log)
Y_f <- lapply(rates_female, log)

names(Y_m)[names(Y_m) == "GBR"] <- names(Y_f)[names(Y_f) == "GBR"] <- "UK"



plot_data_m <- as_tibble(melt(Y_m$ITA, 
                              value.name = "log_rate",
                              varnames = c("Year", "Age")))
plot_data_f <- as_tibble(melt(Y_f$USA, 
                              value.name = "log_rate",
                              varnames = c("Year", "Age")))
plot_data_m %>%
  mutate(title = 'Male - Italy (1933-2020)') %>%
  bind_rows(plot_data_f %>%
              mutate(title = 'Female - United States (1933-2020)')) %>%
  ggplot(aes(x = Age, y = log_rate, group = Year, color = Year)) +
  geom_line() +
  geom_label_repel(data = . %>% 
                     filter((Year == 1943) & 
                              (title == 'Male - Italy (1933-2020)') & 
                              (Age == 40) |
                              (Year == 1933) &
                              (title == 'Female - United States (1933-2020)') & 
                              (Age == 40)),
                   aes(x = Age, y = log_rate, label = Year),
                   nudge_y = 0.9,
                   nudge_x = 1) +
  geom_label_repel(data = . %>% 
                     filter((Year == 2020) & 
                              (title == 'Male - Italy (1933-2020)') & 
                              (Age == 15) |
                              (Year == 2020) &
                              (title == 'Female - United States (1933-2020)') & 
                              (Age == 15)),
                   aes(x = Age, y = log_rate, label = Year),
                   nudge_y = -0.5,
                   nudge_x = 10) +
  facet_grid(cols = vars(title)) +
  scale_color_distiller(palette = 'BrBG') +
  scale_y_continuous(breaks = c(-1,-4,-7,-10)) +
  scale_x_continuous(n.breaks = 10) +
  guides(color = "none") +
  theme(panel.grid = element_blank(),
        text = element_text(size = 12),
        strip.background = element_rect(fill = 'gray85'),
        plot.margin = margin(r = 0),
        axis.title = element_blank()) -> plot_ages

pdf("plot_ages.pdf", width = 12, height = 5)
print(plot_ages)
dev.off()
