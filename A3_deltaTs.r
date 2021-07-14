#^ Скрипт

library(phangorn)
library(TreeTools)

wd <- "C:/MDEV/projects/TreesModF/trees 0.01"
setwd(wd)

conSteps <- 1000

#~ Считываем деревья из файлов и получаем матрицы дистанций
tree_source0 <- read.tree("trees_complete.txt")
tree_source1 <- read.tree("trees_iqtree1.txt")
tree_source2 <- read.tree("trees_iqtree2.txt")
tree_source3 <- read.tree("trees_iqtree3.txt")
tree_source4 <- read.tree("trees_iqtree4.txt")
tree_source5 <- read.tree("trees_iqtree5.txt")

Trat <- integer(conSteps)
Tinf1 <- integer(conSteps)
Tinf2 <- integer(conSteps)
Tinf3 <- integer(conSteps)
Tinf4 <- integer(conSteps)
Tinf5 <- integer(conSteps)
Td1 <- integer(conSteps)
Td2 <- integer(conSteps)
Td3 <- integer(conSteps)
Td4 <- integer(conSteps)
Td5 <- integer(conSteps)
Ed1 <- integer(conSteps)
Ed2 <- integer(conSteps)
Ed3 <- integer(conSteps)
Ed4 <- integer(conSteps)
Ed5 <- integer(conSteps)


s <- 1
for (s in 1:conSteps) {
    tree2 <- reorder(tree_source0[[s]], order = "postorder")
    root_number <- nodepath(tree2)[[1]][1]
    tree_length <- dist.nodes(tree2)[root_number, 1]
    Trat[s] <- tree_length
}

for (s in 1:conSteps) {
    tree2 <- reorder(tree_source1[[s]], order = "postorder")
    root_number <- nodepath(tree2)[[1]][1]
    tree_length <- dist.nodes(tree2)[root_number, 1]
    Tinf1[s] <- tree_length
    Td1[s] <- Trat[s] - Tinf1[s]
}

for (s in 1:conSteps) {
    tree2 <- reorder(tree_source2[[s]], order = "postorder")
    root_number <- nodepath(tree2)[[1]][1]
    tree_length <- dist.nodes(tree2)[root_number, 1]
    Tinf2[s] <- tree_length
    Td2[s] <- Trat[s] - Tinf2[s]
}

for (s in 1:conSteps) {
    tree2 <- reorder(tree_source3[[s]], order = "postorder")
    root_number <- nodepath(tree2)[[1]][1]
    tree_length <- dist.nodes(tree2)[root_number, 1]
    Tinf3[s] <- tree_length
    Td3[s] <- Trat[s] - Tinf3[s]
}

for (s in 1:conSteps) {
    tree2 <- reorder(tree_source4[[s]], order = "postorder")
    root_number <- nodepath(tree2)[[1]][1]
    tree_length <- dist.nodes(tree2)[root_number, 1]
    Tinf4[s] <- tree_length
    Td4[s] <- Trat[s] - Tinf4[s]
}

for (s in 1:conSteps) {
    tree2 <- reorder(tree_source5[[s]], order = "postorder")
    root_number <- nodepath(tree2)[[1]][1]
    tree_length <- dist.nodes(tree2)[root_number, 1]
    Tinf5[s] <- tree_length
    Td5[s] <- Trat[s] - Tinf5[s]
}

maxTd1 <- max(max(Trat), max(Tinf1))
maxTd2 <- max(max(Trat), max(Tinf2))
maxTd3 <- max(max(Trat), max(Tinf3))
maxTd4 <- max(max(Trat), max(Tinf4))
maxTd5 <- max(max(Trat), max(Tinf5))


for (s in 1:conSteps) {
    Ed1[s] <- abs(Td1[s]) / maxTd1
    Ed2[s] <- abs(Td2[s]) / maxTd2
    Ed3[s] <- abs(Td3[s]) / maxTd3
    Ed4[s] <- abs(Td4[s]) / maxTd4
    Ed5[s] <- abs(Td5[s]) / maxTd5
}


library(ggplot2)
library(plyr)

#~ Границы Доверительного интервала
f <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#~ Границы Доверительного интервала ВСЕ ЗНАЧЕНИЯ (0 - 1)
f1 <- function(x) {
  r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#~ Границы линий
f2 <- function(x) {
  r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}




#~ ОШИБКИ
#~ Исключить 5% крайних значений со всех выборок (поставить coni <- 1)
ci1 <- 26
ci2 <- 975

coni <- 1
if (coni == 1) {
  Ed1 <- sort(Ed1)
  Ed1 <- Ed1[ci1:ci2]
  Ed2 <- sort(Ed2)
  Ed2 <- Ed2[ci1:ci2]
  Ed3 <- sort(Ed3)
  Ed3 <- Ed3[ci1:ci2]
  Ed4 <- sort(Ed4)
  Ed4 <- Ed4[ci1:ci2]
  Ed5 <- sort(Ed5)
  Ed5 <- Ed5[ci1:ci2]
}
#~ --------------------------------------------------------------------


boxEs <- data.frame("90 to 60" = Ed1, "55 to 40" = Ed2, "35 to 20" = Ed3, "15 to 5" = Ed4, "Combined" = Ed5)
boxEs_stacked <- stack(boxEs)
pe_text <- ddply(boxEs_stacked, .(ind), summarise, med = median(values), ci0d = quantile(values, probs = 0), ci1d = quantile(values, probs = 1))
pe_text$med <- signif(pe_text$med, 2)
pe_text$ci0d <- signif(pe_text$ci0d, 2)
pe_text$ci1d <- signif(pe_text$ci1d, 2)


maxRange <- max(max(boxEs), min(boxEs))
#limits1 <- maxRange
limits1 <- 0.41

ggplot(boxEs_stacked, aes(x = ind, y = values, fill = ind)) +
    labs(title = "Error T",
#subtitle = "",
    x = "Percent",
    y = "Error") +

    stat_summary(fun.data = f1, geom = "errorbar", alpha = 0.9) +
    stat_summary(fun.data = f2, geom = "linerange", alpha = 0.9) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    ylim(min = 0, max = limits1) +

    geom_text(data = pe_text, aes(x = ind, y = ci0d, label = ci0d), alpha = 0.5, size = 4, vjust = 1.5, color = "red", fontface = 15) +
    geom_text(data = pe_text, aes(x = ind, y = ci1d, label = ci1d), alpha = 0.5, size = 4, vjust = -0.5, color = "blue", fontface = 15) +
    geom_text(data = pe_text, aes(x = ind, y = med, label = med), alpha = 0.7, size = 4, vjust = -0.5, color = "black", fontface = 15) +

    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Dark2") +
    ggsave("boxplot_e2T2.png", width = 5, height = 8)




#~ Боксплот
#~ Исключить 5% крайних значений со всех выборок (поставить coni <- 1)
ci1 <- 26
ci2 <- 975

coni <- 1
if (coni == 1) {
    Td1 <- sort(Td1)
    Td1 <- Td1[ci1:ci2]
    Td2 <- sort(Td2)
    Td2 <- Td2[ci1:ci2]
    Td3 <- sort(Td3)
    Td3 <- Td3[ci1:ci2]
    Td4 <- sort(Td4)
    Td4 <- Td4[ci1:ci2]
    Td5 <- sort(Td5)
    Td5 <- Td5[ci1:ci2]
}
#~ --------------------------------------------------------------------

boxTs2 <- data.frame("90 to 60" = Td1, "55 to 40" = Td2, "35 to 20" = Td3, "15 to 5" = Td4, "Combined" = Td5)
boxTs_stacked2 <- stack(boxTs2)
p_text2 <- ddply(boxTs_stacked2, .(ind), summarise, med = median(values), ci0d = quantile(values, probs = 0), ci1d = quantile(values, probs = 1))
p_text2$med <- round(p_text2$med, 2)
p_text2$ci0d <- round(p_text2$ci0d, 2)
p_text2$ci1d <- round(p_text2$ci1d, 2)

maxRange2 <- max(max(boxTs2), min(boxTs2) * (-1))
#limits2 <- maxRange2
limits2 <- 6000

ggplot(boxTs_stacked2, aes(x = ind, y = values, fill = ind)) +
    labs(title = "Delta T",
#subtitle = "",
    x = "Percent",
    y = "dT") +

    stat_summary(fun.data = f1, geom = "errorbar", alpha = 0.9) +
    stat_summary(fun.data = f2, geom = "linerange", alpha = 0.9) +
    geom_boxplot(alpha = 0.3, outlier.shape = NA) +
    ylim(min = limits2 * (-1), max = limits2) +

    geom_text(data = p_text2, aes(x = ind, y = ci0d, label = ci0d), alpha = 0.5, size = 4, vjust = 1.5, color = "red", fontface = 15) +
    geom_text(data = p_text2, aes(x = ind, y = ci1d, label = ci1d), alpha = 0.5, size = 4, vjust = -0.5, color = "blue", fontface = 15) +
    geom_text(data = p_text2, aes(x = ind, y = med, label = med), alpha = 0.7, size = 4, vjust = -0.5, color = "black", fontface = 15) +


    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Dark2") +
    ggsave("boxplot_dT_onlyCI2.png", width = 5, height = 8)


#~======================================================================



#~ График линейный
graphTs <- data.frame("Original" = Trat, "90 to 60" = Tinf1, "55 to 40" = Tinf2, "35 to 20" = Tinf3, "15 to 5" = Tinf4, "Combined" = Tinf5)
graphTs2 <- graphTs[order(graphTs[, "Original"]),]

ggplot(graphTs2, aes(x = Original, y = X15.to.5)) +
  theme_bw() +
  xlab("Real") +
  ylab("Model") +
  ggtitle("Root ages") +
  geom_smooth(data = graphTs2, aes(x = Original, y = X90.to.60, color = "90% - 60%"), method = "lm", size = 1) +
  geom_smooth(data = graphTs2, aes(x = Original, y = X15.to.5, color = "15% - 5%"), method = "lm", size = 1) +
  geom_smooth(data = graphTs2, aes(x = Original, y = X35.to.20, color = "35% - 20%"), method = "lm", se = FALSE, size = 0.25) +
  geom_smooth(data = graphTs2, aes(x = Original, y = X55.to.40, color = "55% - 40%"), method = "lm", se = FALSE, size = 0.25) +
  geom_smooth(data = graphTs2, aes(x = Original, y = Combined, color = "Combined"), method = "lm", se = FALSE, size = 1)

ggsave("graph_dT.png", width = 8, height = 8)