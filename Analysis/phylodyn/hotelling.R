# install.packages("ICSNP")
library(ICSNP)
library(fmatrix)

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


## generates 1000 trees with 25 tips trees with beta=b
gen_x <- function(b, m = 1000, n = 25, ...){
  x <- rEncod(m = m, n = n, b = b, ...)
  x <- lapply(x, Fmat_from_myencod)
  x <- matrix_list(x)
  x <- t(x)
  x
}

# 
# 
# 
# 
# betas <- c(seq(-1,-.4,by=.2),
#            seq(-.4, -.1, by = .1)[-1],
#            seq(-.1,0, by = .02)[-1],
#            seq(0, 1, by = .1)[-1],
#            10)
# # ind <- 1
# 
# 
# 
# 
# 
# B <- 200
# ret <- matrix(nrow = B, ncol = length(betas))
# 
# 
# set.seed(123)
# for(i in 1:B){
#   cat("\r", i)
#   xlist <- lapply(betas, gen_x)
#   nullx <- gen_x(0)
#   try({
#   ret[i, ] <- sapply(xlist, myHotelling, y = nullx)
#   })
# }
# 
# 
# 
# plot(betas, apply(ret, 2, mean, na.rm = TRUE))
# 
# library(ggplot2)
# p <- ggplot(data = NULL, aes(x = betas, y = apply(ret, 2, mean, na.rm = TRUE))) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = .05, linetype = 2) +
#   labs(x = "beta", y = "Estimated power (200 replicates)",
#        title = "Hotelling test power, null = Kingman") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5))
# p
# ggsave(p, filename = "power_hotelling.png")


