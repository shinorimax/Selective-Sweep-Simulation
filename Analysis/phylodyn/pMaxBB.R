#
# B <- 1e6
#
# FREQ <- 1e3
#
# store <- c()
#
# set.seed(123)
# for(i in 1:B){
#   if(i %% 100 == 0) cat("\r", i, "      ")
#   gen <- e1071::rbridge(end = 1, frequency = FREQ)
#   store <- append(store, max(abs(gen)))
#   if(i %% 1e5 == 0) {
#     saveRDS(store, file = paste0("maxBB_", i, ".RDS"))
#     store <- c()
#   }
# }

loaded <- c()
for(i in 1:10){
  loaded <- append(loaded, readRDS(file = paste0("maxBB_", i, "00000.RDS")))
}

loaded <- sort(loaded)


pMaxBB <- function(x, data = loaded){
  mean(1*(x >= data))
}



