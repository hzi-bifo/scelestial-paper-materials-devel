options(warn=-1)

library(oncoNEM)
library(igraph)
set.seed(7)


calculateEdgeLengths = function(tree, seq) {
    e = as_edgelist(tree)
    seq = t(seq)
    r = c()
    for (i in 1:nrow(e)) {
        r[i] = sum(abs(seq[e[i,1],] - seq[e[i,2],]))
    }
    return(r)
}

writeClones = function(clones, f) {
    close( file( f, open="w" ) )
    for (i in seq_along(clones)) {
        cat(i, " ", clones[[i]], "\n", file=f, append=TRUE);
    }
}

writeTree = function(tree, weight, f) {
    close( file( f, open="w" ) )
    cat(V(tree), "\n", file=f, append=TRUE)
    e = as_edgelist(tree)
    for (i in 1:nrow(e)) {
        cat(e[i,2], "->", e[i,1], " ", weight[i], "\n", sep="", file=f, append=TRUE);
    }
}

writeInput = function(D, f) {
    X = apply(D, 1:2, function(x) { if (x == 0) return("A/A"); if (x == 1) return ("C/C"); return("./.");}  )
    write.table(X, sep = " ", col.names = FALSE, file=f, quote=FALSE)
}

writeInputForBitPhylogeny = function(D, f) {
    X = apply(D, 1:2, function(x) { if (x == 0) return(0); if (x == 1) return (1); return(1);}  )
    X = t(X)
    write.table(X, sep = ",", col.names = TRUE, row.names = FALSE, file=f, quote=FALSE)
}

writeInputForSCITE = function(D, f) {
    X = apply(D, 1:2, function(x) { if (x == 0) return(0); if (x == 1) return (1); return(3);}  )
    # X = t(X)
    write.table(X, sep = " ", col.names = FALSE, row.names = FALSE, file=f, quote=FALSE)
}

## Generate Simulated Data
# simData <- simulateData(N.cells = 100,
#     N.clones = 5,
#     N.unobs = 1,
#     N.sites = 30,
#     FPR = 0.2,
#     FNR = 0.1,
#     p.missing = 0.2)

options(warn=0)

args <- commandArgs(TRUE)
inputFile = args[1]
wd = args[2]
treeFile = args[3]
cloneFile = args[4]


D <- read.csv(inputFile, sep = " ", header = FALSE)
D <- D[,apply(D, 2, function(x) { sum(!is.na(x)) > 0 })]
D <- apply(D, 1:2, function(v) { if (v == 3) return(2); return(v);})

fpr.est <- 0.2
fnr.est <- 0.1

## initialize oncoNEM
oNEM <- oncoNEM$new(Data = D,
    FPR = fpr.est,
    FNR = fnr.est)
oNEM$search(delta = 200)
oNEM.expanded <- expandOncoNEM(oNEM,epsilon = 10,delta = 200,
    checkMax = 10000,app = TRUE)
oncoTree <- clusterOncoNEM(oNEM = oNEM.expanded,
    epsilon = 10)
post <- oncoNEMposteriors(tree = oncoTree$g,
    clones = oncoTree$clones,
    Data = oNEM$Data,
    FPR = oNEM$FPR,
    FNR = oNEM$FNR)

edgeLengths = colSums(post$p_theta)[-1]




## Write results
# setwd("~/Documents/HZI/sc-deep-imputation/test/test-setup-1/")
setwd(wd)

# writeInput(simData$D, "test-2/input-impute.txt")
# writeInputForBitPhylogeny(simData$D, "test-2/input-bp.txt")
# writeInputForSCITE(simData$D, "test-2/input-scite.txt")
# writeClones(simData$clones, "test-2/true-clones.txt")
# writeTree(simData$g, calculateEdgeLengths(simData$g, simData$gtyp), "test-2/true-tree.txt")

# writeClones(oncoTree$clones, "test-3/onconeme-clones.txt")
# writeTree(oncoTree$g, edgeLengths, "test-3/onconeme-tree.txt")
writeClones(oncoTree$clones, cloneFile)
writeTree(oncoTree$g, edgeLengths, treeFile)



