# LOAD PACKAGES
library(plyr)
library(tidyverse)
library(graph4lg)
library(cowplot)

# LOAD IN DATA
SampleInfo <- read.delim(paste0(WORKDIR, "/all_phase31.psam"), header=FALSE, comment.char="#")
Het <- read.delim(paste0(WORKDIR, "/all_phase31.het"), header=FALSE, comment.char="#")
FST <- read.delim(paste0(WORKDIR, "/all_phase31.fst.summary"), header=FALSE, comment.char="#")

# RESHAPE FST DATA INTO A SYMMATRIX SQUARE MATRIX
FST_MAT <- FST %>% 
  spread(V2, V3) %>%
  mutate(ACB = as.numeric(NA), .before = 2) %>%
  add_row(V1 = "YRI", YRI = 0, .after = 25) %>%
  column_to_rownames(var = "V1") %>%
  replace(., col(.) == row(.), 0) %>%
  as.matrix()

for(i in 1:25){
  for(j in (26-i):1){
    FST_MAT[i+j, i] <- FST_MAT[i, i+j]
  }
}

# NORMALIZE TO RANGE 0 - 1 for GRAPHIPNG PURPOSES
FST_MAT_NORM <- (FST_MAT - min(FST_MAT)) / (max(FST_MAT) - min(FST_MAT))

# CREATE A RANGE NORMALIZED HET VARIABLE
Het$V5_NORM <- (Het$V5 - min(Het$V5)) / (max(Het$V5) - min(Het$V5))

# WRANGLE AND JOIN HET AND SAMPLE TABLES
Het <- Het %>%
  rename(ID = V1, HET = V5, HET_NORM = V5_NORM) %>%
  select(all_of(c("ID", "HET", "HET_NORM")))

SampleInfo <- SampleInfo %>%
  rename(ID = V1, SUPERPOP = V5, POP = V6) %>%
  select(all_of(c("ID", "SUPERPOP", "POP"))) %>%
  right_join(Het)

# CREATE PopInfo and summarzie means of Het scores
PopInfo <- SampleInfo %>%
  group_by(POP) %>%
  summarize(HET_MEAN = mean(HET),
            HET_MED = median(HET),
            HET_NORM_MEAN = mean(HET_NORM),
            HET_NORM_MED = median(HET_NORM)) %>%
  arrange(HET_NORM_MED) %>%
  mutate(ORDER = as.numeric(row_number()),
         HETEROZYGOSITY = 1 - HET_NORM_MED) %>%
  left_join(unique(select(SampleInfo, all_of(c("POP", "SUPERPOP")))))

# REORDER FST MATRIX
FST_MAT <- reorder_mat(FST_MAT, order = PopInfo$POP)

# MELT MATRIX INTO DATAFRAME
FST_MAT[lower.tri(FST_MAT, diag = TRUE)] <- NA

FST <- as.data.frame(FST_MAT) %>%
  rownames_to_column(var = "X") %>%
  gather(key = "Y", value = "FST", MSL:CHS) %>%
  left_join(select(PopInfo, POP, ORDER), by = c("X" = "POP")) %>%
  rename(X.1 = ORDER) %>%
  left_join(select(PopInfo, POP, ORDER), by = c("Y" = "POP")) %>%
  rename(Y.1 = ORDER) %>%
  arrange(X.1, Y.1) %>%
  filter(complete.cases(.)) %>%
  mutate(FST_NORM = (FST - min(FST))/(max(FST)-min(FST)))

# PLOT THE FST DATA
W <- ggplot(PopInfo, aes(POP, y="SUP")) +
  geom_tile(aes(fill = SUPERPOP), color = "#ffffff") +
  scale_fill_manual(values = c("#fdffb6", "#caffbf", "#9EDDFF", "#bdb2ff", "#ffadad")) +
  xlim(PopInfo$POP) + 
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Superpopulation:") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust= 0.5, margin = margin(t = -20, b = 20), size = 9),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "#e9ecef"))

X <- ggplot(FST, aes(X, Y)) +
  geom_tile(aes(fill = FST_NORM), color = "#ffffff") +
  scale_fill_gradient(low="#ffffff", high="#9b2226", breaks = c(0.5)) +
  xlim(PopInfo$POP) + 
  ylim(rev(PopInfo$POP)) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Normalized FST:") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "#e9ecef"))

Y <- ggplot(PopInfo, aes(POP, y=c("HET"))) +
  geom_tile(aes(fill = HETEROZYGOSITY), color = "#ffffff") +
  scale_fill_gradient(low="#ffffff", high="#023e8a", breaks = c(0, 0.5, 1)) +
  xlim(PopInfo$POP) + 
  theme_bw() +
  labs(x = NULL, y = NULL, fill = "Heterozygosity:") +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "#e9ecef")) 

WL <- get_legend(W)
XL <- get_legend(X)
YL <- get_legend(Y)

X <- X + theme(legend.position = "none")
Y <- Y + theme(legend.position = "none")
W <- W + theme(legend.position = "none")

L <- plot_grid(XL, YL, WL, nrow = 3, align = "v", rel_heights = c(7, 3, 7))
P <- plot_grid(X, Y, W,  nrow = 3, align = "v", rel_heights = c(12, 1, 1.2))

plot_grid(P, L, ncol = 2, align = "h", axis ="t", rel_widths = c(7, 1))
