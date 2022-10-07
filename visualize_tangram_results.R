library(data.table)
library(ggplot2)

################################################
# Layout of the cells distributed over 10 runs #
################################################

mapping_SlideSeq <- fread("data/mapping_df_SlideSeq.csv")
mapping_STARmap <- fread("data/mapping_df_STARmap.csv")
mapping_SlideSeq$run <- as.factor(mapping_SlideSeq$run)
mapping_STARmap$run <- as.factor(mapping_STARmap$run)
mapping_SlideSeq$probability[mapping_SlideSeq$probability < 0.0001] <- NA

dir.create("visualizations")

g <- ggplot(mapping_SlideSeq[cell %in% c(245,1146,8186,14441,15867)], aes(x=x, y=y))+
  geom_point(aes(color=probability, alpha = probability), size=0.8)+
  facet_grid(cell ~ run)+
  #theme_minimal()+
  scale_y_reverse()+
  scale_color_gradient(low = "navy", high = "red", na.value = NA)+
  theme(text = element_text(size = 15), line = element_blank(), axis.text = element_blank(), panel.background = element_blank(),
        rect = element_blank())+
  labs(x = "Run", y = "Cell ID")
g
ggsave("visualizations/mapping_slideseq.pdf", plot=g)

mapping_high <- mapping_STARmap[probability > 0.5]
background <- mapping_STARmap[run == 0, ]
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
set.seed(1234)
random_cell <- sort(sample(unique(mapping_high$cell), 16))
ggplot(mapping_high[cell %in% random_cell], aes(x=x, y=y, color = run, alpha = probability))+
  geom_point()+
  facet_wrap(~ cell, ncol = 4)+
  theme_void()+
  scale_color_manual(values = cols)+
  theme(text = element_text(size = 20))
ggsave("visualizations/mapping_starmap.pdf", width=8, height=6)

################################################
# Correlation plots                            #
################################################

cor_sm <- fread("data/correlations_mappings_STARmap.csv")
colnames(cor_sm) <- c("V1", "Correlation")
cor_sm$dataset <- "STARmap"
cor_sm$matrix <- "Mapping"
cor_sm <- cor_sm[, c(1,4,2,3)]

cor_slideseq <- fread("data/correlations_mappings_SlideSeq.csv")
colnames(cor_slideseq) <- c("V1", "Correlation")
cor_slideseq$dataset <- "SlideSeq"
cor_slideseq$matrix <- "Mapping"
cor_slideseq <- cor_slideseq[, c(1,4,2,3)]

expr_sm <- fread("data/correlations_expr_STARmap.csv")
colnames(expr_sm) <- c("V1", "Expression", "Expression")
expr_sm <- melt(expr_sm[, -1], variable.name = "matrix", value.name = "Correlation")
expr_sm$dataset <- "STARmap"
expr_sm$matrix <- "Expression"

expr_slideseq <- fread("data/correlations_expr_SlideSeq.csv")
colnames(expr_slideseq) <- c("V1", "Expression", "Expression")
expr_slideseq <- melt(expr_slideseq[, -1], variable.name = "matrix", value.name = "Correlation")
expr_slideseq$dataset <- "SlideSeq"
expr_slideseq$matrix <- "Expression"

cors <- rbindlist(list(cor_sm[, -1], cor_slideseq[, -1], expr_sm, expr_slideseq))

ggplot(cors, aes(x = matrix, y = Correlation, color = dataset)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Matrix", y = "Pairwise Correlation")+
  theme(text = element_text(size = 15))+
  ylim(0,1)
ggsave("visualizations/correlations.pdf", height = 3, width = 8)
