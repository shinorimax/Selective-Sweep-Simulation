
library(fmatrix)

# setwd("~/Documents/Paromita/")
# setwd("~/Documents/Paromita/Frechet_Trees/Rcode/")

source("pMaxBB.R")

#simulates F matrices from BD, uniform or aldous models
#what is the BF distribution?
gen_x <- function(b, m = 100, n = 7, ...){
  x <- rEncod(m = m, n = n, b = b, ...)
  x <- lapply(x, Fmat_from_myencod)
  x <- matrix_list(x)
  x <- t(x)
  x
}

#computes the euclidean distance
dist_x_mu <- function(x,mu){
  diff <- (t(x) - mu)
  distribution <-  apply(diff, 2, function(u){sum(u**2)})
  return(distribution)
}

#Two-sample Kolmogorov-Smirnov statistic, difference between the two empirical distributions
KSstat <- function(x, y){
  ## This portion of code extracted from stat::ks.test
  n.x <- length(x)
  n.y <- length(y)
  w <- c(x, y)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))

  if (length(unique(w)) < (n.x + n.y)) {
    z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
  }

  max(abs(z))
}


cdf<-function(u,x){
  return(sum(x<u)/length(x))
}

KSstat2 <- function(x, y){
  ## This portion of code extracted from stat::ks.test
  n.x <- length(x)
  n.y <- length(y)
  w <- c(x, y)
  u<-unique(sort(w))
  for (i in 1:length(u))
  z<-sapply(1:length(u),function(i){cdf(u[i],x)})-sapply(1:length(u),function(i){cdf(u[i],y)})
  #z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
  
  #if (length(unique(w)) < (n.x + n.y)) {
   # z <- z[c(which(diff(sort(w)) != 0), n.x + n.y)]
  #}
  
 return(max(abs(z)))
}
#

pvalue <- function(stat, x, B = 5e3, verbose = TRUE){
  n <- length(x)
  cdf <- rank(x)/n

  freq <- 1/n

  pval <- 0
  if(verbose) cat("simulating brownian bridge processes\n")
  for(i in 1:B){
    if(verbose) cat("\r", i)
    gen <- e1071::rbridge(end = 1, frequency = 1/freq)
    pval <- (i-1)/i*pval + 1/i*(max(abs(gen)) >= stat)
  }
  if(verbose) cat("\ndone!\n")
  return(pval)
}

alt_pvalue<-function(stat,x){
  n <- length(x)
  cdf <- rank(x)/n
  
  freq <- 1/n
  
  pval <- 0
  if(verbose) cat("simulating brownian bridge processes\n")
  for(i in 1:B){
    if(verbose) cat("\r", i)
    gen <- e1071::rbridge(end = 1, frequency = 1/freq)
    pval <- (i-1)/i*pval + 1/i*(max(abs(gen)) >= stat)
  }
}

myHotelling <- function(x, y, ...){
  o1 <- apply(x, 2, sd) == 0
  o2 <- apply(y, 2, sd) == 0
  o3 <- apply(x, 2, mean) == apply(y, 2, mean)

  o <- which(o1 & o2 & o3)

  try({
    ret <- ICSNP::HotellingsT2(x[,-o], y[,-o], ...)$p.value < 0.05
    return(ret)
  })

  return(NA)
}


newTest <- function(x, y, flag = "test", ...){

  n.x <- nrow(x)
  n.y <- nrow(y)

  mu_x <- apply(x, 2, mean)
  mu_y <- apply(y, 2, mean)


  distxx <- dist_x_mu(x, mu_x)
  distxy <- dist_x_mu(x, mu_y)
  distyx <- dist_x_mu(y, mu_x)
  distyy <- dist_x_mu(y, mu_y)

  STATISTIC <-
    0.5*
    sqrt(n.x*n.y/(n.x+n.y))*
    (KSstat2(distxx, distyx) +
       KSstat2(distxy, distyy))

  if(flag == "stat") return(STATISTIC)

  ## Option 1
  # w <- rbind(x, y)
  # mu_w <- apply(w, 2, mean)
  #
  # distww <- dist_x_mu(w, mu_w)
  # pvalue(STATISTIC, distww, verbose = FALSE)

  ## Option 2
  # kolmim::pkolmim(STATISTIC, n.x + n.y)

  ## Option 3
  pval<-1 - pMaxBB(STATISTIC)
  return(list(pval=pval, Statistic=STATISTIC))
}

oldTest <- function(x, y, flag = "test", ...){
  n.x <- nrow(x)
  n.y <- nrow(y)

  mu_x <- apply(x, 2, mean)
  mu_y <- apply(y, 2, mean)

  w <- rbind(x, y)
  mu_w <- apply(w, 2, mean)


  distxx <- dist_x_mu(x, mu_x)
  distxy <- dist_x_mu(x, mu_y)
  distyx <- dist_x_mu(y, mu_x)
  distyy <- dist_x_mu(y, mu_y)

  distww <- dist_x_mu(w, mu_w)

  V1 <- mean(distxx)
  V1c <- mean(distxy)
  V2 <- mean(distyy)
  V2c <- mean(distyx)
  Vw <- mean(distww)

  sigma2 <- mean(distww**2) - Vw**2

  STATISTIC <- n.x*n.y*((V1-V2)**2 + (V1c - V1 + V2c - V2)**2)/sigma2

  if(flag == "stat") return(STATISTIC)


  1 - pchisq(STATISTIC, df = 1)
}

graphTest <- function(x, y, ...){

  n.x <- nrow(x)
  n.y <- nrow(y)

  dist_mat <- dist(rbind(x, y))


  stree <- ade4::mstree(dist_mat, ngmax = 5)

  gtest_out <- gTests::g.tests(stree, 1:n.x, (n.x + 1:n.y))

  c(gtest_out$original$pval.approx,
    gtest_out$generalized$pval.approx,
    gtest_out$weighted$pval.approx,
    gtest_out$maxtype$pval.approx)
}

betas <- c(seq(-1,-.4,by=.2),
           seq(-.4, -.1, by = .1)[-1],
           seq(-.1,0, by = .02)[-1],
           seq(0, 1, by = .1)[-1],
           10, 50)
# ind <- 1


# betas <- betas[1:2]


B <- 100
ret <- reth <- reto <-
  retg1 <- retg2 <- retg3 <- retg4 <-
  matrix(nrow = B, ncol = length(betas))



m_samplesize <- 1000
n_tips <- 20

  # 
  # set.seed(123)
  # for(i in 1:B){
  # 
  #   cat("\r", i)
  # 
  #   xlist <- lapply(betas, gen_x, m = m_samplesize, n = n_tips)
  #   nullx <- gen_x(0, m = m_samplesize, n = n_tips)
  # 
  #   ilist <- lapply(xlist, function(u){
  #     N <- nrow(u)
  #     dist <- matrix(NA, N, N)
  #     for(i in 1:N) {for (j in 1:i){
  #       dist[i,j] <- dist[j,i] <-
  #         mean((u[i,] - u[j])**2)}}
  #     which.min(apply(dist,1,mean))
  #   })
  # 
  #   try({
  #     ret[i, ] <- sapply(xlist, newTest, y = nullx)
  # 
  #     reth[i, ] <- sapply(xlist, myHotelling, y = nullx)
  # 
  #     reto[i, ] <- sapply(xlist, oldTest, y = nullx)
  # 
  #     retg <- sapply(xlist, graphTest, y = nullx)
  # 
  #     retg1[i, ] <- retg[1,]
  #     retg2[i, ] <- retg[2,]
  #     retg3[i, ] <- retg[3,]
  #     retg4[i, ] <- retg[4,]
  #   })
  # }
  # 
  # ALPHA <- .05
  # 
  # 
  # fret_n <- ret < ALPHA
  # fret_h <- reth
  # fret_o <- reto > ALPHA
  # fret_g1 <- retg1 < ALPHA
  # fret_g2 <- retg1 < ALPHA
  # fret_g3 <- retg1 < ALPHA
  # fret_g4 <- retg1 < ALPHA
  # 
  # 
  # 
  # library(ggplot2)
  # 
  # geom_thisfret <- function(u, label){
  #   list(geom_line(aes(y = apply(u, 2, mean, na.rm = TRUE),
  #                      colour = label)),
  #        geom_point(aes(y = apply(u, 2, mean, na.rm = TRUE),
  #                       colour = label)))
  # }
  # colours <- c("Hotelling" = "black",
  #              "New Test" = "red",
  #              # "G1" = "yellow",
  #              # "G2" = "orange",
  #              # "G3" = "purple",
  #              "Graph-based" = "brown")
  # 
  # 
  # p <- ggplot(data = NULL, aes(x = betas)) +
  #   labs(x = "beta", y = sprintf("Estimated power (%d replicates)", B),
  #        title = "Test power, null = Kingman") +
  #   theme_bw() +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         legend.position = "bottom") +
  #   geom_hline(yintercept = .05, linetype = 2) +
  #   geom_thisfret(fret_h, "Hotelling") +
  #   geom_thisfret(fret_n, "New Test") +
  #   # geom_thisfret(fret_g1, "G1") +
  #   # geom_thisfret(fret_g2, "G2") +
  #   # geom_thisfret(fret_g3, "G3") +
  #   geom_thisfret(fret_g4, "Graph-based") +
  #   scale_colour_manual(values = colours) +
  #   xlim(c(NA, 1))+
  #   NULL
  # 
  # p
  # # p <- ggplot(data = NULL, aes(x = betas, y = apply(ret2, 2, mean, na.rm = TRUE))) +
  # #   geom_line() +
  # #   geom_point() +
  # #   geom_hline(yintercept = .05, linetype = 2) +
  # #   labs(x = "beta", y = "Estimated power (5000 replicates)",
  # #        title = "BB test power, null = Kingman") +
  # #   theme_bw() +
  # #   theme(plot.title = element_text(hjust = 0.5))
  # # p
  # # ggsave(p, filename = "power_BB.png")
  # 
  # # p2 <- ggplot(data = NULL, aes(x = betas)) +
  # #   geom_line(aes(y = apply(ret2, 2, mean, na.rm = TRUE)), colour = "red") +
  # #   geom_point(aes(y = apply(ret2, 2, mean, na.rm = TRUE)), colour = "red") +
  # #   geom_line(aes(y = apply(reth, 2, mean, na.rm = TRUE)), colour = "black") +
  # #   geom_point(aes(y = apply(reth, 2, mean, na.rm = TRUE)), colour = "black") +
  # #   geom_line(aes(y = apply(reto2, 2, mean, na.rm = TRUE)), colour = "blue") +
  # #   geom_point(aes(y = apply(reto2, 2, mean, na.rm = TRUE)), colour = "blue") +
  # #   geom_hline(yintercept = .05, linetype = 2) +
  # #   labs(x = "beta", y = "Estimated power (500 replicates)",
  # #        title = "Test power, null = Kingman") +
  # #   theme_bw() +
  # #   theme(plot.title = element_text(hjust = 0.5))
  # # p2 + xlim(NA, 10)
  # 
  # 
  # # save(list = ls(), file = "dist_based_test_20210824.RData")
  # 
  # 
  # 
  # 
  # B <- 5000
  # 
  # retnull <- matrix(nrow = B, ncol = 1)
  # for(i in 1:B){
  # 
  #   cat("\r", i)
  # 
  #   nullx <- gen_x(0, m = 100, n = 7)
  #   nully <- gen_x(0, m = 100, n = 7)
  # 
  #   try({
  #     retnull[i, ] <- newTest(nullx, nully, flag = "stat")
  #   })
  # }
  # 
  # ##Analysis of SARS-Cov-2 Washington State
  # #Sample 1: All 170 sequences obtained from Gisaid, high coverage 2022-05-24
  #  #         These sequences were aligned with the following line:
  #   #        mafft --thread -1 HCOV_Washington_2023_11_12_20.fasta > HCOV_Washingtont_out.fasta
  #   #        The xml assumed mutation rate of .00074 per year: Washington_May_short.xml. 
  #   #          Generated 5000000 posterior samples and thinned every 5000 to obtain 1001 samples
  #           
  # ##Prepare F-matrices
  # ##Posterior distribution of trees
  # library("ape")
  # trees<-read.nexus(file="~/Documents/Covid/data/Washington_May_short-HCOV_Washingtont_out.trees")
  # #trees1<-read.nexus(file="~/Documents/Covid/data/trees1.trees")
  # length(trees)
  # trees_small<-trees[seq(1,5001,by=5)]
  # 
  # 
  # library("fmatrix")
  # listF<-list()
  # for (j in 1:length(trees_small)){
  #   listF[[j]]<-gen_Fmat(trees_small[[j]])
  #   print(j)
  #   print(dim(listF[[j]]))
  #   if (j==1){M<-listF[[j]]}
  #   M<-M+listF[[j]]
  # }
  # 
  # saveRDS(listF,file="listF_May")
  # M<-M/length(listF)
  # listFx<-matrix_list(listF)
  # listFx<-t(listFx)
  # 
  # ##Simulate from Kingman
  # listFKing<-gen_x(0, m =1000 , n = 170)
  # saveRDS(listFKing,file="listFKing")
  # #test
  # test_out<-newTest(listFx,listFKing)
  # 
  # 
  # 
  # tree_mcc<-read.nexus(file="~/Documents/Covid/data/MCC_Washington_May")



