## Make simulations by drawing parameters from uniform distributions.

satsimcont <- function(ntaxa = c(8,16,32), imbal = c(T, F), trlen = c(0.001, 0.65), stemminess = c(0.1, 0.9), seqlen = c(250, 2500), modelsim = c("JC", "GTR+G"), modelest = "JC", gammacat = 4, gammashape = 1, qmat = c(1.3472, 4.8145, 0.9304, 1.2491, 5.5587, 1), basef = c(0.2628, 0.2605, 0.2436, 0.2331)){

       imbal <- sample(imbal, 1)
       ntaxa <- sample(ntaxa, 1)
       trlen <- runif(1, min = trlen[1], max = trlen[2])
       stemminess <- runif(1, min = stemminess[1], max = stemminess[2])
       seqlen <- runif(1, min = seqlen[1], max = seqlen[2])
       modelsim <- sample(modelsim, 1)
       
       simparams <- c(ntaxa, imbal, trlen, stemminess, seqlen, modelsim)
       names(simparams) <- c("ntaxa", "imbal", "trlen", "stemminess", "seqlen", "modelsim")

       if(imbal) tr <- trs[paste0("imbal", ntaxa)][[1]] else tr <- trs[paste0("bal", ntaxa)][[1]]
       Ned <- Nedge(tr)

       # The following line makes the trlen setting be interpreted as the mean branch length.
       trlen <- trlen * Ned

       Ninted <- Nedge(tr)-Ntip(tr)
       tr$edge.length[which(tr$edge[,2] %in% 1:Ntip(tr))] <- rep((trlen*(1-stemminess)) / Ntip(tr), Ntip(tr))
       tr$edge.length[which(!tr$edge[,2] %in% 1:Ntip(tr))] <- rep((trlen*stemminess) / Ninted, Ninted)
       
       if(modelsim == "JC"){
                al <- simSeq(tr, l = seqlen)
       } else if(modelsim == "GTR+G"){
                rates <- phangorn::discrete.gamma(gammashape, gammacat)
                dat <- lapply(rates, function(x) simSeq(tr, l = round(seqlen/gammacat), Q = qmat, bf = basef, rate = x))
                al <- do.call(cbind, dat)
       }

       al <- as.DNAbin(al)
       
       # This version of the function does not allow completely random sites in alignments.
       #Nrandsites <- round(seqlen*proprandsites)
       #if(proprandsites > 0) al[,1:Nrandsites] <- sample(unique(al), Nrandsites*Ntip(tr), replace = T)

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
       res$simparams <- simparams
       return(res)

}