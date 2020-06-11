#' @import purrr
#' @import stats
#' @import utils
#' @import parallel
#' @import future
#' @import furrr
#' @import generics
#' @import vroom
#' @import readr
#' @aliases NULL
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Construct the linear regression with bootstraps and each kind of parament might used.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#'
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param m m parts user want to split the data
#' @param B the number of replications
#' @param cl the number of cores that user want to use.
#' @export


blblm <- function(formula, data, m = 10, B = 5000,cl=1) {
  if (class(data) == "data.frame"){
    data_list <- split_data(data, m)
    n=nrow(data)
  }else{
    data_list<- map(data,read_csv,col_types="dd")
    n=length(vroom_lines(data)-length(data))
  }
  suppressWarnings(plan(multiprocess, worker=cl))
  estimates <- future_map(data,
                          ~lm_each_subsample(formula = formula, data = data, n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
#'
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param m m parts user want to split the data

split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' compute the estimates
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param n the number of random vectors to draw
#' @param B the number of replications

lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' estimate the regression estimates based on given the number of repetitions
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param n the number of random vectors to draw

lm1 <- function(formula, data, n) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  fit <- lm(formula, data, weights = freqs)
  list(coef = coef(fit), sigma = blbsigma(fit))
}



#' compute sigma from fit
#'
#' @param fit formula expression
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' Get the expression of formula
#'
#' @param x linear regression formula
#' @param ... further arguments passed from or to other methods
#' @method print blblm

print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$out))
  cat("\n")
}

#' Get the confidence interval of parament
#'
#' @param object estimated object
#'
#' @param confidence Whether return the confidence interval of sigma
#' @param level alpha value
#' @param ... further arguments passed from or to other methods
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
     alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Get the coefficient of linear regression
#'
#' @param object estimated object
#'
#' @param ... further arguments passed from or to other methods
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Get the confidence interval of parament of linear regression
#'
#' @param object The data user defined.
#'
#' @param parm The initial string
#' @param level The 95% confidence interval
#' @param ... further arguments passed from or to other methods
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$call), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' Based on the previous dataset to predict the new data
#'
#' @param object estimated object
#'
#' @param new_data the new data wanted to predict
#' @param confidence confidence interval
#' @param level alpha value
#' @param ... further arguments passed from or to other methods
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


#' Get the confidence interval of mean
#'
#' @param x The dataset
#' @param level The alpha value
mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

#' Get the mean
#'
#' @param .x A list or atomic vector.
#' @param .f A function,formula or vector
#' @param ... further arguments passed from or to other methods
map_mean <- function(.x, .f,...) {
  (map(.x, .f,...) %>% reduce(`+`)) / length(.x)
}

#' Columnbind
#'
#' @param .x the dataset
#' @param .f A function,formula or vector
#' @param ... further arguments passed from or to other methods
map_cbind <- function(.x, .f,...) {
  map(.x, .f,...) %>% reduce(cbind)
}

#' Rowbind
#'
#' @param .x the dataset
#' @param .f A function,formula or vector
#' @param ... further arguments passed from or to other methods
map_rbind <- function(.x, .f,...) {
  map(.x,.f,...) %>% reduce(rbind)
}


#' estimate the logistic estimates based on given number of repetitions
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param freqs how many multinomial probabilities exists
#'
gm<- function(formula,data,freqs){
  environment(formula)<- environment()
  fit<- glm(formula,family = gaussian,data,weights=freqs)
  list(coef = coef(fit), sigma=fit$residuals)
}

#' compute the logistic regression estimates for a blb dataset
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param n the number of random vectors to draw

gm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  gm(formula, data, freqs)
}

#'
#'
#' compute the estimates for logistic regression
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param n the number of random vectors to draw
#' @param B the number of replications
gm_each_subsample <- function(formula, data, n, B) {
  replicate(B, gm_each_boot(formula, data, n), simplify = FALSE)
}


#' Construct the logistic regression with bootstraps and each kind of parament might used.
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#'
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param m m parts user want to split the data
#' @param B the number of replications
#' @param cl the number of cores that user want to use.
#' @export

blbgm <-function(formula, data, m = 10, B = 5000,cl=1) {
  if (class(data) == "data.frame"){
    data_list <- split_data(data, m)
    n=nrow(data)
  }else{
    data_list<- map(data,read_csv,col_types="dd")
    n=length(vroom_lines(data)-length(data))
  }
  data_list <- split_data(data, m)
  suppressWarnings(plan(multiprocess, worker=cl))
  estimates <- future_map(data_list,
                          ~gm_each_subsample(formula = formula, data = data, n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbgm"
  invisible(res)
}
