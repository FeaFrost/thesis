#^ Скрипт объединяет все деревья полученные в iqtree в 1 файл trees/trees_iqtree.txt

library(phangorn)

conSteps <- 1000
sss <- 1
for (sss in 1:5) {
    wd <- paste0("C:/MDEV/projects/TreesModF/results/Nex", sss)
    setwd(wd)

    fn6 <- file(paste0("C:/MDEV/projects/TreesModF/trees/trees_iqtree", sss, ".txt"), open = "w")

    s <- 1
    for (s in 1:conSteps) {
        #treefile <- read.tree(paste0(s, "s_DNAResult.txt.timetree.nwk"))
        treefile <- read.nexus(paste0(s, "s_DNAResult.txt.timetree.nex"))
        write.tree(treefile, file = fn6)
        #file.remove(paste0(s, "s_DNAResult.txt.timetree.nwk"))
    }

    close(fn6)
}