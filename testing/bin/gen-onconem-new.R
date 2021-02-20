options(warn=-1)

library(oncoNEM)
library(igraph)

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


options(warn=0)

args <- commandArgs(TRUE)
wd = args[1]
sciteFile = args[2]
cloneFile = args[3]
treeFile = args[4]
seqFile = args[5]

cellCnt = as.numeric(args[6])
cloneCnt = as.numeric(args[7])
siteCnt = as.numeric(args[8])
set.seed(as.numeric(args[9]))
unobserved = as.numeric(args[10])
fpr = as.numeric(args[11])
fnr = as.numeric(args[12])
mvr = as.numeric(args[13])

## Generate Simulated Data
simData <- simulateData(N.cells = cellCnt,
    N.clones = cloneCnt,
    N.unobs = unobserved,
    N.sites = siteCnt,
    FPR = fpr,
    FNR = fnr,
    p.missing = mvr,
    randomizeOrder = TRUE)

## Write results
# setwd("~/Documents/HZI/sc-deep-imputation/test/test-setup-1/")
setwd(wd)

# writeInput(simData$D, "test-2/input-impute.txt")
# writeInputForBitPhylogeny(simData$D, "test-2/input-bp.txt")
writeInputForSCITE(simData$D, sciteFile)
writeClones(simData$clones, cloneFile)
writeTree(simData$g, calculateEdgeLengths(simData$g, simData$gtyp), treeFile)

# writeClones(oncoTree$clones, "test-3/onconeme-clones.txt")
# writeTree(oncoTree$g, edgeLengths, "test-3/onconeme-tree.txt")
#writeClones(oncoTree$clones, cloneFile)
#writeTree(oncoTree$g, edgeLengths, treeFile)



