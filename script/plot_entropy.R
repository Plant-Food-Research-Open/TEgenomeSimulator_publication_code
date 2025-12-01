#!/software/statistical/R-4.3.3/bin/Rscript
# Read CSV with read.table
entr0g <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/original/entropy/sliding_entropy_genome.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr1g <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m2/entropy/sliding_entropy_genome.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr2g <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m0/entropy/sliding_entropy_genome.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr3g <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m0_full/entropy/sliding_entropy_genome.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr4g <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/Garlic/entropy/sliding_entropy_genome.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr0t <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/original/entropy/sliding_entropy_teonly.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr1t <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m2/entropy/sliding_entropy_teonly.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr2t <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m0/entropy/sliding_entropy_teonly.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr3t <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m0_full/entropy/sliding_entropy_teonly.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr4t <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/Garlic/entropy/sliding_entropy_teonly.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr0nt <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/original/entropy/sliding_entropy_nonte.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr1nt <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m2/entropy/sliding_entropy_nonte.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr2nt <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m0/entropy/sliding_entropy_nonte.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr3nt <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/tair10_ch1_m0_full/entropy/sliding_entropy_nonte.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
entr4nt <- read.table(file.path(WRK, 'pub/output/TE_sim_eval/Garlic/entropy/sliding_entropy_nonte.csv'), header = TRUE, sep = ",", stringsAsFactors = FALSE)
# Reshape to long format
entr0g <- entr0g |> mutate(Category = "genome", simulator = "Ori")
entr1g <- entr1g |> mutate(Category = "genome", simulator = "m2")
entr2g <- entr2g |> mutate(Category = "genome", simulator = "m2+m0")
entr3g <- entr3g |> mutate(Category = "genome", simulator = "m0full")
entr4g <- entr4g |> mutate(Category = "genome", simulator = "Garlic")
entr0t <- entr0t |> mutate(Category = "TEonly", simulator = "Ori")
entr1t <- entr1t |> mutate(Category = "TEonly", simulator = "m2")
entr2t <- entr2t |> mutate(Category = "TEonly", simulator = "m2+m0")
entr3t <- entr3t |> mutate(Category = "TEonly", simulator = "m0full")
entr4t <- entr4t |> mutate(Category = "TEonly", simulator = "Garlic")
entr0nt <- entr0nt |> mutate(Category = "nonTE", simulator = "Ori")
entr1nt <- entr1nt |> mutate(Category = "nonTE", simulator = "m2")
entr2nt <- entr2nt |> mutate(Category = "nonTE", simulator = "m2+m0")
entr3nt <- entr3nt |> mutate(Category = "nonTE", simulator = "m0full")
entr4nt <- entr4nt |> mutate(Category = "nonTE", simulator = "Garlic")
# Combine
entrdf <- rbind(entr0g, entr1g, entr2g, entr3g, entr4g)
entrdf <- rbind(entrdf, entr0t, entr1t, entr2t, entr3t, entr4t)
entrdf <- rbind(entrdf, entr0nt, entr1nt, entr2nt, entr3nt, entr4nt)
entrdf$simulator <- factor(entrdf$simulator, levels = c("Ori", "Garlic", "m2", "m2+m0", "m0full"))
# Bar plot
plot <- ggplot(entrdf |> filter(simulator != "m0full"),
               aes(x = pos/1e6, y = entropy, fill = simulator)) +   # <- fill inside aes
  geom_col(alpha = 0.6) +
  scale_fill_manual(values = c('#99CCCC', '#FFCC66', '#CC6666', '#99CC99')) +
  facet_grid(vars(simulator), vars(Category)) +
  labs(title = "Shannon entropy",
       x = "Position (Mb)",
       y = "Entropy (bits)") +
  theme(title = element_text(size = 8),
        axis.text.y = element_text(size = rel(1.3)),
        axis.text.x = element_text(size = rel(1.3)),
        axis.title.x = element_text(size = rel(1.4)),
        axis.title.y = element_text(size = rel(1.4)),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(1)),
        legend.key.size = unit(0.5, "cm"),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"))

png(filename = paste0(OUT, "entropy_spectrum.png"), width = plotwidth, height = plotheight, units = "in", res = 1200)
print(plot)
while (!is.null(dev.list()))  dev.off()