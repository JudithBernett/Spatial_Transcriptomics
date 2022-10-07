library(data.table)
library(ggplot2)

################################################
# Layout of the cells distributed over 10 runs #
################################################

mapping_df <- fread("~/PycharmProjects/SpatialTranscriptomics/dfs/mapping_df2.csv")
#mapping_df <- fread("~/PycharmProjects/SpatialTranscriptomics/dfs/mapping_df.csv")
mapping_df$run <- as.factor(mapping_df$run)
mapping_df$probability[mapping_df$probability < 0.0001] <- NA
#mapping_df <- mapping_df[complete.cases(mapping_df)]

ggplot(mapping_df, aes(x=probability, fill = run))+
  geom_histogram()+
  scale_y_log10()+
  facet_wrap(~run)+
  theme_minimal()

q <- mapping_df[, quantile(probability, na.rm=TRUE), by = run]
q$quantile <- rep(c(0,25,50,75,100), 10)
q$quantile <- as.factor(q$quantile)
ggplot(q, aes(x = quantile, y = V1))+
  geom_boxplot()+
  theme_minimal()+
  scale_y_log10()+
  labs(y = "probability")

g <- ggplot(mapping_df[cell %in% c(245,1146,8186,14441,15867)], aes(x=x, y=y))+
  geom_point(aes(color=probability, alpha = probability), size=0.8)+
  facet_grid(cell ~ run)+
  #theme_minimal()+
  scale_y_reverse()+
  scale_color_gradient(low = "navy", high = "red", na.value = NA)+
  theme(text = element_text(size = 15), line = element_blank(), axis.text = element_blank(), panel.background = element_blank(),
        rect = element_blank())+
  labs(x = "Run", y = "Cell ID")
g
ggsave("~/Downloads/mapping_slideseq.pdf", plot=g)

mapping_high <- mapping_df[probability > 0.5]
background <- mapping_df[run == 0, ]
background$probability <- 0
background$run <- as.numeric(background$run)
background$run <- -1
mapping_high <- rbind(mapping_high, background)
mapping_high[run == -1, ]$run <- "background"

cols <- c("background" = "#cfcfcf",
         "0" = "#000000",
         "1" = "#E69F00",
         "2" = "#56B4E9",
         "3" = "#009E73",
         "4" = "#F0E442",
         "5" = "#0072B2",
         "6" = "#D55E00",
         "7" = "#CC79A7",
         "8" = "#CC6666",
         "9" = "#9999CC")
ggplot(mapping_high[cell %in% c(122, 239, 753, 1375, 1457, 1914, 2433, 2999, 7594, 7933, 9116, 10622, 10932, 11061, 12923, 13506)], aes(x=x, y=y, color = run, alpha = probability))+
  geom_point()+
  facet_wrap(~ cell, ncol = 4)+
  theme_void()+
  scale_color_manual(values = cols)+
  theme(text = element_text(size = 20))
ggsave("~/Downloads/mapping_starmap.pdf")

################################################
# Correlation plots                            #
################################################

cor_sm <- fread("~/Downloads/ST/correlations_all_starmap.csv")
colnames(cor_sm) <- c("V1", "Correlation")
#cor_sm$Subset <- "Training Genes"
cor_sm$dataset <- "STARmap"
cor_sm$matrix <- "Mapping"
cor_sm <- cor_sm[, c(1,4,2,3)]

cor_slideseq <- fread("~/Downloads/ST/correlations_all_slideseq.csv")
colnames(cor_slideseq) <- c("V1", "Correlation")
#cor_slideseq$Subset <- "Training Genes"
cor_slideseq$dataset <- "SlideSeq"
cor_slideseq$matrix <- "Mapping"
cor_slideseq <- cor_slideseq[, c(1,4,2,3)]


expr_sm <- fread("~/Downloads/ST/correlations_expr_starmap.csv")
colnames(expr_sm) <- c("V1", "Expression", "Expression")
expr_sm <- melt(expr_sm[, -1], variable.name = "matrix", value.name = "Correlation")
expr_sm$dataset <- "STARmap"
#expr_sm$matrix <- "Expression"

expr_slideseq <- fread("~/Downloads/ST/correlations_expr_slideseq.csv")
colnames(expr_slideseq) <- c("V1", "Expression", "Expression")
expr_slideseq <- melt(expr_slideseq[, -1], variable.name = "matrix", value.name = "Correlation")
expr_slideseq$dataset <- "SlideSeq"
#expr_slideseq$matrix <- "Expression"

cors <- rbindlist(list(cor_sm[, -1], cor_slideseq[, -1], expr_sm, expr_slideseq))

ggplot(cors, aes(x = matrix, y = Correlation, color = dataset)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Matrix", y = "Pairwise Correlation")+
  theme(text = element_text(size = 15))
ggsave("~/Downloads/correlations.pdf", height = 3, width = 8)



probabilities <- fread("~/PycharmProjects/SpatialTranscriptomics/probabilities0_1.csv")
q0 <- quantile(probabilities$run_0)
q1 <- quantile(probabilities$run_1)
q0 <- data.table(quantile = names(q0), probabilities = q0, run = 0)
q1 <- data.table(quantile = names(q1), probabilities = q1, run = 1)
q <- rbind(q0, q1)
q$quantile <- factor(q$quantile, levels = c("0%", "25%", "50%", "75%", "100%"))
ggplot(q, aes(x = quantile, y = probabilities))+
  geom_boxplot()+
  theme_minimal()+
  scale_y_log10()+
  labs(y = "probability")


ggplot(probabilities[sample(nrow(probabilities), nrow(probabilities) * 0.3)], aes(x=run_0))+
  geom_histogram()+
  scale_y_log10()+
  theme_minimal()
