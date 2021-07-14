#^ iqtree2 скрипт для получения деревьев

library(phangorn)
library(TreeTools)

conSteps <- 1000

wd <- "C:/MDEV/projects/TreesModF/results"
setwd(wd)


sss <- 1

#~ Запускаем IQTree и генерируем деревья
for (sss in 1:5) {
    s <- 1
    count <- 0

    while (s <= 1000) {
        count <- count + 1
        print(paste0("step: ", count))

        x <- system(paste0("iqtree2.exe -s C:/MDEV/projects/TreesModF/results/", count, "s_DNAResult.txt --date C:/MDEV/projects/TreesModF/trees/Dates", sss, "/Date", count, ".txt -m JC --clock-sd 0.001 --date-outlier 30 --date-tip 0 -fast -nt 4 -redo -quiet"), wait = TRUE, show.output.on.console = FALSE)

        suppressWarnings( c(file.remove(paste0(count, "s_DNAResult.txt.bionj")),
                            file.remove(paste0(count, "s_DNAResult.txt.ckp.gz")),
                            file.remove(paste0(count, "s_DNAResult.txt.iqtree")),
                            file.remove(paste0(count, "s_DNAResult.txt.log")),
                            file.remove(paste0(count, "s_DNAResult.txt.mldist")),
                            file.remove(paste0(count, "s_DNAResult.txt.timetree.lsd")),
                            file.remove(paste0(count, "s_DNAResult.txt.uniqueseq.phy")),
                            file.remove(paste0(count, "s_DNAResult.txt.timetree.nwk")),
                            file.remove(paste0(count, "s_DNAResult.txt.treefile"))))

        if (x != 0) {
            file.remove(paste0(s, "s_DNAResult.txt.timetree.nex"))
            print(paste0("NEXT! s = ", s))
            next
        }
        file.rename(paste0(count, "s_DNAResult.txt.timetree.nex"), paste0("C:/MDEV/projects/TreesModF/results/Nex", sss, "/", s, "s_DNAResult.txt.timetree.nex"))
        s <- s + 1
    }
}

