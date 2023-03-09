Shannon.gen_edit <- function (gInd, estimator = NULL) 
{
  MLE = function(X) {
    X = X[X > 0]
    n = sum(X)
    -sum(X/n * log(X/n))
  }
  Z = function(X) {
    X = X[X > 0]
    Y = X[X > 1]
    n = sum(X)
    -n * sum(X/n * log(X/n)) - (n - 1)/n * sum((n - X) * 
                                                 (-X/(n - 1) * log(X/(n - 1)))) - (n - 1)/n * sum(-Y * 
                                                                                                    (Y - 1)/(n - 1) * log((Y - 1)/(n - 1)))
  }
  CS = function(X) {
    x = X
    x = x[x > 0]
    n = sum(x)
    f1 = sum(x == 1)
    C_head = 1 - f1/n
    a = -sum(C_head * (x/n) * log(C_head * (x/n))/(1 - (1 - 
                                                          C_head * (x/n))^n))
    a
  }
  Ch = function(X) {
    x = X
    x = x[x > 0]
    n = sum(x)
    UE <- sum(x/n * (digamma(n) - digamma(x)))
    f1 <- sum(x == 1)
    f2 <- sum(x == 2)
    if (f1 > 0) {
      A <- 1 - ifelse(f2 > 0, (n - 1) * f1/((n - 1) * 
                                              f1 + 2 * f2), (n - 1) * f1/((n - 1) * f1 + 2))
      B = sum(x == 1)/n * (1 - A)^(-n + 1) * (-log(A) - 
                                                sum(sapply(1:(n - 1), function(k) {
                                                  1/k * (1 - A)^k
                                                })))
    }
    if (f1 == 0) {
      B = 0
    }
    if (f1 == 1 & f2 == 0) {
      B = 0
    }
    UE + B
  }
  inputs <- lapply(adegenet::seppop(gInd), 
                    # Don't split by locus
                   #function(i) sapply(adegenet::seploc(i), 
                                                              function(j) colSums(j@tab, na.rm = TRUE))#)
  out <- list()
  if (is.null(estimator)) 
    estimator <- "z"
  for (i in estimator) {
    if (i %in% c("z", "sh", "cs", "ch")) 
      print(paste(i, "Letter code OK"))
    else {
      print(paste("Wrong letter code:", i))
      stop("<estimator> option is not correct. \n          For details see help: ?Shannon.gen")
    }
  }
  if ("sh" %in% estimator) {
    out[["Shannon_1949"]] <- as.data.frame(lapply(inputs, 
                                                  function(i) sapply(i, MLE)))
  }
  if ("z" %in% estimator) {
    out[["Zahl_1977"]] <- as.data.frame(lapply(inputs, function(i) sapply(i, 
                                                                          Z)))
  }
  if ("cs" %in% estimator) {
    out[["Chao_Shen_2003"]] <- as.data.frame(lapply(inputs, 
                                                    function(i) sapply(i, CS)))
  }
  if ("ch" %in% estimator) {
    out[["Chao_et_al_2013"]] <- as.data.frame(lapply(inputs, 
                                                     function(i) sapply(i, Ch)))
  }
  return(out)
}
