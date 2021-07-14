#^ Скрипт подготавливает ДАТА-файлы для iqtree

library(phangorn)
library(TreeTools)

wd <- "C:/MDEV/projects/TreesModF/trees"
setwd(wd)

nichNumb <- 200
conSteps <- 2000

#~ Считываем деревья из файлов и получаем матрицы дистанций
names_order <- c(paste0("sp_", 1:nichNumb))
tree_source <- read.tree("trees_complete.txt")


edges <- c(0.9, 0.6,    #~ 1 Выборка 90% - 60%
           0.55, 0.4,   #~ 2 Выборка 55% - 40%
           0.35, 0.2,   #~ 3 Выборка 35% - 20%
           0.15, 0.05)  #~ 4 Выборка 15% - 5%

sss <- 4    #~ Номер использыемой выборки
number <- 4 #~ Сколько берём случайных узлов

for (sss in 1:(length(edges) / 2)) {
    s <- 1
    for (s in 1:conSteps) {
        tree2 <- reorder(tree_source[[s]], order = "postorder")
        #plot(tree2)
        #nodelabels()
        root_number <- nodepath(tree2)[[1]][1]
        tree_length <- dist.nodes(tree2)[root_number, 1]
        non_number <- c(  (root_number + 1) : (  root_number + length( nodepath(tree2) ) - 2  )  )
        edge_time <- dist.nodes(tree2)[root_number, (root_number + 1):(root_number + length(nodepath(tree2)) - 2)]
        edge_time <- tree_length - edge_time
        data_time <- data.frame(non_number = non_number, edge_time = edge_time)

        data_times <- data_time[(edge_time <= edges[sss * 2 - 1] * tree_length & edge_time >= edges[sss * 2] * tree_length),]

        #~ Увеличиваем окно выборки в сторону листьев, если нет ни одного узла
        mult <- 0
        while (length(data_times[,1]) == 0) {
            mult <- mult + 1
            data_times <- data_time[(edge_time <= edges[sss * 2 - 1] * tree_length & edge_time >= edges[sss * 2 + mult] * tree_length),]
        }
        #~ ----------------------------------------------------------------------

        fn5 <- file(paste0("C:/MDEV/projects/TreesModF/trees/Dates", sss, "/Date", s, ".txt"), open = "w")

        number2 <- number
        if (length(data_times[, 1]) < number) {
        number2 <- length(data_times[, 1])
        }
        fours <- data_times[sample(1:length(data_times[, 1]), number2),]

        iii <- 1
        for (iii in 1:number2) {
            tree3 <- Preorder(tree2)
            sub_tree <- Subtree(tree3, fours[iii,1])
            taxons <- c(sub_tree$tip.label)
            cat(taxons, sep = ",", file = fn5)
            write(paste0(" -",fours[iii, 2]), fn5)
        }
        close(fn5)
    }
}

#^ Для смешанной выборки узлов
number <- 1 #~ Сколько берём случайных узлов в каждой выборке
s <- 1
for (s in 1:conSteps) {
  tree2 <- reorder(tree_source[[s]], order = "postorder")
  # plot(tree2)
  # nodelabels()
  root_number <- nodepath(tree2)[[1]][1]
  tree_length <- dist.nodes(tree2)[root_number, 1]
  non_number <- c((root_number + 1):(root_number + length(nodepath(tree2)) - 2))
  edge_time <- dist.nodes(tree2)[root_number, (root_number + 1):(root_number + length(nodepath(tree2)) - 2)]
  edge_time <- tree_length - edge_time
  data_time <- data.frame(non_number = non_number, edge_time = edge_time)

  sss <- 1
  fn5 <- file(paste0("C:/MDEV/projects/TreesModF/trees/Dates5/Date", s, ".txt"), open = "w")
  for (sss in 1:(length(edges) / 2)) {
    data_times <- data_time[(edge_time <= edges[sss * 2 - 1] * tree_length & edge_time >= edges[sss * 2] * tree_length),]

    #~ Увеличиваем окно выборки в сторону листьев, если нет ни одного узла
    mult <- 0
    while (length(data_times[, 1]) == 0) {
      mult <- mult + 1
      data_times <- data_time[(edge_time <= edges[sss * 2 - 1] * tree_length & edge_time >= edges[sss * 2 + mult] * tree_length),]
    }
    #~ ----------------------------------------------------------------------

    number2 <- number
    if (length(data_times[, 1]) < number) {
      number2 <- length(data_times[, 1])
    }
    fours <- data_times[sample(1:length(data_times[, 1]), number2),]

    iii <- 1
    for (iii in 1:number2) {
      tree3 <- Preorder(tree2)
      sub_tree <- Subtree(tree3, fours[iii, 1])
      taxons <- c(sub_tree$tip.label)
      cat(taxons, sep = ",", file = fn5)
      write(paste0(" -", fours[iii, 2]), fn5)
    }
  }
  close(fn5)
}