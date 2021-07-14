
#* ANCHOR 1. Считываем файл с настройками ============================
wd <- "C:/MDEV/projects/TreesModF"
setwd(wd)
getwd()


#~ Настраиваем отображение чисел (1000 вместо 1e+3)
options(scipen = 10)

"Options:"
opts <- data.frame(Option = c( "DNA length",
                               "Mut.probability",
                               "Confidence Interval",
                               "Conf int step number",
                               "Number of niches",
                               "Jump probability",
                               "Jump probability increasing by",
                               "Jump probability increasing to"),
    Value = read.table("options.txt", sep = "\t", header = FALSE))

opts <- opts[, c(1, 3:ncol(opts))]
print(opts)

dnaSize <- opts[1, 2]
mutProb <- opts[2, 2]
conInt <- opts[3, 2]
conSteps <- opts[4, 2]
nichNumb <- opts[5, 2]
jumpProb <- opts[6, 2]
#* 1. ================================================================




#* ANCHOR 2. Считываем матрицы дистанций и сохраняем деревья =========
library(phangorn)
print("Matrix uploading started..")

setwd("C:/MDEV/projects/TreesModF/results")
m <- mutProb
n <- nichNumb
j <- jumpProb

#~ Открываем файлы для записи деревьев
fn1 <- file("C:/MDEV/projects/TreesModF/trees/trees_names.txt", open = "w")
fn2 <- file("C:/MDEV/projects/TreesModF/trees/trees_complete.txt", open = "w")

#~ Основной цикл
s <- 1
for (s in 1:conSteps) {
    #^ ANCHOR 2.1. Считываем матрицу из файла
    matrixPath <- paste0(s,"s_matrixResult.txt")
    matr <- matrix(NA, n, n)
    x <- scan(matrixPath)
    #file.remove(matrixPath)
    matr[upper.tri(matr, diag = TRUE)] <- x
    matr <- t(matr)

    rownames(matr) <- c(   paste0(  "sp_", c( 1 : (n) )  )   )
    colnames(matr) <- c(   paste0(  "sp_", c( 1 : (n) )  )   )

    #^ ANCHOR 2.2. Преобразуем в дерево и сохраняем на диск
    tree <- upgma(matr, method = "single")
    #plot(tree)
    #ltt.plot(tree)
    #fit <- ltt.plot.coords(tree)

    xx <- c(m, n, j, s)
    write(xx, file = fn1, sep = " ")
    write.tree(tree, file = fn2)
    #yy <- read.tree("trees_complete.txt", tree.names=TRUE)
}
close(fn1)
close(fn2)

#* 2. ================================================================




#* ANCHOR 3. Расчитывем дистанции DNA и сохраняем деревья ============
wd <- "C:/MDEV/projects/TreesModF/results"
setwd(wd)

#~ Открываем файл для записи деревьев
fn3 <- file("C:/MDEV/projects/TreesModF/trees/trees_upgma.txt", open = "w")

#^ реконструкция древа upgma
filo <- function(pas) {
  d <- dist.dna(pas, model = "JC")
  tr <- upgma(d)
  return(tr)
}

#~ Основной цикл
s <- 1
for (s in 1:conSteps) {
    DNAPath <- paste0(s,"s_DNAResult.txt")
    pas <- read.dna(DNAPath, format = "fasta")
    fit_upgma <- filo(pas)
    #edgelabels(signif(fit_upgma$edge.length, digits = 2))

    write.tree(fit_upgma, file = fn3)
}
#fit <- read.tree(fn2)
#d <- cophenetic(fit)

#~ Закрываем файл для записи деревьев
close(fn3)
#* 3. ================================================================

#* ANCHOR 4. Получаем матрицы дистанций ==============================
wd <- "C:/MDEV/projects/TreesModF/trees"
setwd(wd)

#~ Считываем деревья из файлов и получаем матрицы дистанций
library(ape)
library(nlme)
library(drc)
library(aomisc)

conSteps <- 1000
nichNumb <- 200

names_order <- c(paste0("sp_", 1:nichNumb))
tree1 <- read.tree("trees_complete.txt")

#~ Выбрать только 1 строку из 2-х
tree2 <- read.tree("trees_upgma.txt") #~ Для сравнения с деревом построенном по дистанциям DNA
#tree2 <- read.tree("trees_iqtree.txt")      #~ Для сравнения с деревом IQtree
#tree2 <- read.nexus("Onedate/1s_DNAResult.txt.timetree.nex") #~ Для сравнения с деревом IQtree
#~ ------------------------------

fn4 <- file("C:/MDEV/projects/TreesModF/regression/bics.txt", open = "w")

# bicsdf <- data.frame( f1 = double(),  mich = double(),  f2 = double(),  f3 = double(),
#                       log = double(), sqrt = double(), exp = double(), nexp = double(),
#                       asym = double() )

bicsdf <- data.frame(f1 = double(), f2 = double(), f3 = double(),
                      log = double(), sqrt = double())


write(c("f1", "f2", "f3", "log", "sqrt"), file = fn4, ncolumns = 5, sep = "\t")
#write(c("f1", "mich", "f2", "f3", "log", "sqrt", "exp", "nexp", "asym"), file = fn4, ncolumns = 9,, sep = "\t")

#bicsdf <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), c("f1", "mich", "f2", "f3", "log", "sqrt", "exp", "nexp", "asym"))

s <- 1
sss <- 1
for (s in 1:conSteps) {
  dmatrix1 <- cophenetic(tree1[[s]])
  dmatrix2 <- cophenetic(tree2[[s]])

  dmatrix1 <- dmatrix1[names_order, names_order]
  dmatrix2 <- dmatrix2[names_order, names_order]

  #~ Переводим нижнюю половину матриц без диагоналей в векторы c1 и c2
  c1 <- dmatrix1[lower.tri(dmatrix1, diag = FALSE)]
  c2 <- dmatrix2[lower.tri(dmatrix2, diag = FALSE)]

  #~ Создаём data.frame
  final_df <- data.frame(c1, c2)

  #~ Сортируем по 1-му столбцу
  final_df <- final_df[order(final_df[, "c1"]),]

  #final_df <- final_df[!duplicated(final_df[, "c1"]),]
  #write.table(final_df, "x.txt", sep = "\t")
  # final_df <- final_df[!duplicated(round(final_df[, "c2"], 8)),]

  final_df <- final_df[!(final_df[, "c1"] %in% 0),]
  # final_df <- final_df[!(final_df[, "c2"] %in% 0),]
  #* 4. ================================================================




  #* ANCHOR 5. Проводим регрессионный анализ ===========================


  # #~ Вводим функции для различных кривых
  # #^ Логарифмическая зависимость y=a*log(x)+b

  #~ Подготавливаем данные для анализа
  data <- final_df
  x <- data$c1
  y <- data$c2
  di <- data.frame(x = x, y = y)

  #m_data <- di
  #k <- 2

  # ^ Для выборки из 100 элементов
  # diold <- di
  # di <- diold
  # di2 <- di[sample(1:length(di[, 1]), 100),]
  # di2 <- di2[order(di2[, "x"]),]
  # colnames(di) <- c("S", "v")
  #^ -----------------------------

  #di <- di[!(di[, "x"] %in% 0),]

  model_1 <- lm(y ~ x, data = di)
  #model_mich <- nls(y ~ Vm * x / (K + x), data = di, start = list(K = 4200, Vm = 0.09))
  #model_mich <- drm( y ~ x, fct = MM.2(), data = di )
  model_2 <- lm(y ~ poly(x, 2), data = di)
  model_3 <- lm(y ~ poly(x, 3), data = di)
  #model_log <- lm( exp(y) ~ x, data = di )
  model_log <- lm(y ~ log(x), data = di)
  #model_exp <- nls(  y ~ a * (exp(k * x) * k), data = di, start = list( a = 0, k = ) )
  #model_exp <- drm( y ~ x, fct = DRC.expoDecay(), data = di )
  model_sqrt <- lm(y ~ I(sqrt(x)), data = di)
  #model_nexp <- drm( y ~ x, fct = DRC.negExp(), data = di )
  #model_asym <- drm(y ~ x, fct = DRC.asymReg(), data = di)


  bicsc <- c(BIC(model_1), BIC(model_2), BIC(model_3),
                BIC(model_log), BIC(model_sqrt))

  # bicsc <- c(BIC(model_1), BIC(model_mich), BIC(model_2), BIC(model_3),
  #             BIC(model_log), BIC(model_sqrt), BIC(model_exp), BIC(model_nexp),
  #             BIC(model_asym))

  #     library(rsq)
  #     rsq(model_exp)
  #     rsq(model_exp, adj = TRUE)
  #     BIC(model_sqrt)
  #     summary(model_log)

  # s <- 100
  # di2 <- di
  # plot(di2$x, di2$y)
  # lines(di2$x, di2$y, col = "blue")

  # newx <- seq(min(di2$x), max(di2$x), length.out = 100)
  # new.df <- data.frame(x = newx)

  # preds <- predict(model_1, new.df, interval = "confidence", level = 0.95)
  # lines(newx, preds[, 1], col = "red")



  # #~------------
  #     plot(di$x, di$y)
  #     lines(di$x, di$y, col = "blue")

  #     newx <- seq(min(di$x), max(di$x), length.out = 100)
  #     new.df <- data.frame(x = newx)

  #     preds <- log(predict(model_log, new.df, interval = "confidence", level = 0.95))
  #     lines(newx, preds[, 1], col = "red")




  #^ Добавляем новые строки в дата фрейм с BIC
  #bicsdf[nrow(bicsdfw) + 1,] <- list(bicsc)

  write(bicsc, file = fn4, ncolumns = 5, sep = "\t")
  # write(bicsc, file = fn4, ncolumns = 9, sep = "\t")

  print(paste("Step ", s, "completed"))
}

close(fn4)

#wd <- "C:/MDEV/projects/TreesModF/regression"
#setwd(wd)
#write.table(bicsdf, "bics.txt", sep = "\t", row.names = FALSE)

#* 5. ================================================================

library(ggplot2)
library(plyr)

wd <- "C:/MDEV/projects/TreesModF/regression 0.01"
setwd(wd)

#~ 0.01
allbics <- as.data.frame(read.table("bics.txt", sep = "\t", header = TRUE))
allbics_stacked <- stack(allbics)

#~ 0.005
ab2 <- as.data.frame(read.table("bics.txt", sep = "\t", header = TRUE))
ab2_stacked <- stack(ab2)

#~ 0.015
ab3 <- as.data.frame(read.table("bics.txt", sep = "\t", header = TRUE))
ab3_stacked <- stack(ab3)

p_meds <- ddply(allbics_stacked, .(ind), summarise, med = median(values))
p_meds2 <- ddply(ab2_stacked, .(ind), summarise, med = median(values))
p_meds3 <- ddply(ab3_stacked, .(ind), summarise, med = median(values))


ggplot(allbics_stacked, aes(x = ind, y = values, fill = ind)) +
    labs(title = "BIC's compararison",
#subtitle = "",
    x = "Functions",
    y = "BIC") +
    #~ 0.005 RED
    geom_boxplot(alpha = 0.1, data = p_meds2, aes(x = ind, y = med), color = "red", linetype = 1, size = 0.05) +
    geom_boxplot(alpha = 0.05, data = ab2_stacked, aes(x = ind, y = values), fill = "red", outlier.alpha = 0, linetype = 0) +
    #geom_text(data = p_meds2, aes(x = ind, y = med, label = med), size = 3, vjust = -0.5, color = "red") +
    #~ 0.015 BLUE
    geom_boxplot(alpha = 0.1, data = p_meds3, aes(x = ind, y = med), color = "blue", linetype = 1, size = 0.05) +
    geom_boxplot(alpha = 0.05, data = ab3_stacked, aes(x = ind, y = values), fill = "blue", outlier.alpha = 0, linetype = 0) +
    #geom_text(data = p_meds2, aes(x = ind, y = med, label = med), size = 3, vjust = -0.5, color = "red") +

    #~ 0.01 COLOR
    geom_boxplot(alpha = 0.3) +
    geom_text(data = p_meds, aes(x = ind, y = med, label = med), size = 3, vjust = -0.5) +
    theme(legend.position = "none") +
    scale_fill_brewer(palette = "Dark2")
ggsave("boxplot_iq100.png", width = 8, height = 8)

ggplot(di, aes(x = x, y = y, colour = Functions)) +
  theme_bw() +
  xlab("T rat") +
  ylab("T dna") +
  ggtitle("Regressions") +
  geom_point(alpha = 0.6, size = 1, colour = "black") +
  theme(legend.title = element_text(colour = "black", size = 10, face = "bold")) +
  geom_smooth(data = model_1, aes(x = x, y = y, colour = "Linear"), formula = "y ~ x + 0", method = "lm", se = FALSE) +
  geom_smooth(data = model_2, aes(x = x, y = y, colour = "Poly2"), formula = "y ~ poly(x, 2)", method = "lm", se = FALSE) +
  geom_smooth(data = model_3, aes(x = x, y = y, colour = "Poly3"), formula = "y ~ poly(x, 3)", method = "lm", se = FALSE) +
  geom_smooth(data = model_log, aes(x = x, y = y, colour = "Log"), formula = "y ~ log(x)", method = "lm", se = FALSE) +
  geom_smooth(data = model_sqrt, aes(x = x, y = y, colour = "Sqrt"), formula = "y ~ I(sqrt(x))", method = "lm", se = FALSE) +
  ylim(-0.02, 0.12)
#geom_line(data = model_1, aes(x = x, y = y), colour = "red")
ggsave("di.png", width = 8, height = 8)



quit()