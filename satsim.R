#This function simulates data under a given condition, and runs ML on it, returning the tree and parameter estimates.

load("trees.Rdata")
source("runIQtree.R")
source("saturation.index.R")
source("get.ci.test.R")
source("get.compression.R")

satsim <- function(ntaxa = c(16, 32, 64), imbal = T, trlen = 0.01, stemminess = 0.5, seqlen = c(500, 1500, 2500), proprandsites = c(0, 0.3, 0.6, 0.9), modelsim = c("JC", "GTR+G"), modelest = "JC", gammacat = 4, gammashape = 1, qmat = c(1.3472, 4.8145, 0.9304, 1.2491, 5.5587, 1), basef = c(0.2628, 0.2605, 0.2436, 0.2331)){
       if(imbal) tr <- trs[paste0("imbal", ntaxa)][[1]] else tr <- trs[paste0("bal", ntaxa)][[1]]
       Ned <- Nedge(tr)
       
       # The following line can be removed. It makes the trlen be the mean branch length.
       trlen <- trlen * Ned
       
       Ninted <- Nedge(tr)-Ntip(tr)
       tr$edge.length[which(tr$edge[,2] %in% 1:Ntip(tr))] <- rep((trlen*(1-stemminess)) / Ntip(tr), Ntip(tr))
       tr$edge.length[which(!tr$edge[,2] %in% 1:Ntip(tr))] <- rep((trlen*stemminess) / Ninted, Ninted)       

       ##tr$edge.length <- rep(trlen / Ned, Ned)
       
       if(modelsim == "JC"){
       		al <- simSeq(tr, l = seqlen)
       } else if(modelsim == "GTR+G"){
       	        rates <- phangorn::discrete.gamma(gammashape, gammacat)
		dat <- lapply(rates, function(x) simSeq(tr, l = round(seqlen/gammacat), Q = qmat, bf = basef, rate = x))
		al <- do.call(cbind, dat)
       }
       
       al <- as.DNAbin(al)
       Nrandsites <- round(seqlen*proprandsites)
       if(proprandsites > 0) al[,1:Nrandsites] <- sample(unique(al), Nrandsites*Ntip(tr), replace = T)
       
       res <- runIQtree(sdata = al, format = "bin", aadata = F, temp_name = paste0("tempdat.", sample(1:10000, 1), ".phy"), iqtreePath = "iqtree", model = modelest, tree = NULL)
       res$trueTree <- tr
       ent <- try(calculate_index(al))
       ci <- try(ci.test(res$tree, al))
       com <- try(get.compression(res$tree, al))
       res$estEntr <- ent[[1]]
       res$expEntr <- ent[[2]]
       res$pvalEntr <- ent[[3]]
       res$tEntr <- ent[[4]]
       if(class(ci) != "try-error") res$cit <- ci$t else res$cit <- NA
       if(class(ci) != "try-error") res$cip <- ci$p else res$cip <- NA
       if(class(com) != "try-error") res$tCompress <- com$t else res$tCompress <- NA
       if(class(com) != "try-error") res$pCompress <- com$p else res$pCompress <- NA
       return(res)
}
