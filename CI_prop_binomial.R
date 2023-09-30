#' Binomial confidence intervals for proportions
#'
#' This function uses `binom.test` to derive estimates of the confidence interval 
#' for proportions.
#'
#' @param k number of successes
#'
#' @param n number of trials
#'
#' @param conf confidence level (1 - alpha), defaults to 0.95
#'
#' @param result a `character` vector indicating which bounds of the confidence
#'   interval to return: can be `"lower"`, or `"upper"`, or `"both"`
#'
#' @param perc a `logical` indicating if results should be formatted as
#'   percentages, rounded to 2 decimal places, or not (defaults to `FALSE`)
#'
#' @param dec the number of decimal places used for rounding percentages
#'
#' @author Thibaut Jombart
#' 
#' @examples
#' 
#' ## CI for 1/10
#' > prop_ci(1, 10)
#'              [,1]
#' lower 0.002528579
#' upper 0.445016117
#' 
#' ## the function is vectorised, so 'k' and 'n' can be vectors;
#' ## shortest vector is recycled
#' > prop_ci(0:10, 10)
#'            [,1]        [,2]       [,3]       [,4]      [,5]     [,6]      [,7]
#' lower 0.0000000 0.002528579 0.02521073 0.06673951 0.1215523 0.187086 0.2623781
#' upper 0.3084971 0.445016117 0.55609546 0.65245285 0.7376219 0.812914 0.8784477
#'            [,8]      [,9]     [,10]     [,11]
#' lower 0.3475471 0.4439045 0.5549839 0.6915029
#' upper 0.9332605 0.9747893 0.9974714 1.0000000
#'

prop_ci <- function(k, n,
                    result = c("both", "lower", "upper"),
                    perc = FALSE,
                    conf = 0.95,
                    dec = 2) {
  if(n == 0){
    out <- c(0,1)
  } else{
    out <- binom.test(k, n, conf.level = conf)$conf.int
    result <- match.arg(result)
  }
  if (result == "both") {
    result <- c("lower", "upper")
  }
  names(out) <- c("lower", "upper")
  out <- out[result]
  
  if (perc) {
    out <- round(100 * out, dec)
  }
  out
}

prop_ci <- Vectorize(prop_ci)