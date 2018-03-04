##' simulation and plots explaining KS-test and Wilcox-test
##' @rdname ks.test_expl
##' @name ks.test_expl
##' @author Vitalii Kleshchevnikov
##' @description This function generates random data (rnorm, or other provided by distr argument) or uses user-provided data to explain how KS-test and Wilcox-test measure difference in distributions. It shows empirical cumulative distribution of z-scores for a given group of genes (set) and for all other genes (other). User can compare this plot to the KS-test and the Wilcox-test results.
##' @param N_other integer, size of the "other" group, ignored if x is provided
##' @param N_set integer, size of the "set" group, ignored if y is provided
##' @param alternative character, "greater" or "less". The choice of alternative hypothesis for KS-test \code{\link[stats]{ks.test}}. Greater tests is the set (y) distribution (or, in case of Wilcox test, it's median) is shifted to the right compated to other (x) distribution. This choice of "greater" or "less" is flipped for the Wilcox test due to difference in the null hypothesis formulation (\code{\link[stats]{wilcox.test}}).
##' @param z_score numeric, size of the true difference
##' @param prop_pos numeric, between 0 and 1, the fraction of "set" group that is higher than other by \code{z_score}
##' @param prop_neg numeric, between 0 and 1, the fraction of "set" group that is lower than other by \code{z_score}
##' @param x_coor_other numeric, x coordinate of "other" group label on the plot
##' @param x_coor_set numeric, x coordinate of "set" group label on the plot
##' @param cex numeric, plotting parameter, size of the group label text, axis label text and axis annotation text
##' @param title character, title of the plot
##' @param xlab character, x-axis label
##' @param lab_other character, other group label
##' @param lab_set character, set group label
##' @param seed integer, to initialise random number generator for reproducibility
##' @param x numeric vector, the distribution of values of the "other" group. These can come from the experimental data.
##' @param y numeric vector, the distribution of values of the "set" group. These can come from the experimental data.
##' @param distr function that generates a vector of random numbers of lenght n (n argument to this function). rnorm by default
##' @param ... other plotting parameters \code{\link[graphics]{plot}}
##' @return list containing the output of KS-test and Wilcox-test
##' @example # Let's have a look at a TF that regulates
##' # some of it's targets positively, some negatively and some are not affected.
##' ks.test_expl(alternative = "greater", prop_pos = 0.33, prop_neg = 0.4, seed = 1, cex.main = 1.4)
##' @export ks.test_expl
ks.test_expl = function(N_other = 1000, N_set = 100, alternative = "greater", z_score = 1.5, prop_pos = 0.3, prop_neg = 0, x_coor_other = -1.7, x_coor_set = 1.7, cex = 1.7, title = "Treatment: anti-Oct4 shRNA.\nOct4 is a pos. and a neg. regulator of its targets", xlab = "z-score, simulated", lab_other = "other genes", lab_set = "Oct4 targets", seed = 1, x = NULL, y = NULL, distr = rnorm, ...) {
  set.seed(seed)

  #data
  if(is.null(x)) random1 = distr(n = N_other) else random1 = x
  random2 = distr(n = N_set)
  random3 = sample(c(
    rep(z_score, ceiling(N_set*prop_pos)),
    rep(-z_score, ceiling(N_set*prop_neg)),
    rep(0, ceiling(N_set*(1 - prop_pos - prop_neg)))))
  if(is.null(y)) random23 = random2+random3 else random23 = y

  #plots
  par(mar = c(5,5,3,2.5))
  plot(ecdf(random1), xlab = xlab, ylab = "ECDF",
       main = title,
       cex.lab = cex, cex.axis = cex, ...)
  lines(x = c(median(random1),
              median(random1)),
        y = c(0,1), lwd = 2, col = "black")
  lines(x = c(median(random23),
              median(random23)),
        y = c(0,1), lwd = 2, col = "blue")
  lines(ecdf(random23), col = "blue")
  text(x = x_coor_other, y = 0.6, lab_other, col = "black", cex = cex)
  text(x = x_coor_set, y = 0.4, lab_set, col = "blue", cex = cex)

  # tests
  test_res = list()
  if(alternative == "greater") {
    alternative_ks = "greater"
    alternative_wilcox = "less"
  } else if(alternative == "less") {
    alternative_ks = "less"
    alternative_wilcox = "greater"
  }
  test_res$ks = ks.test(random1, random23, alternative = "greater")
  test_res$wilcox = wilcox.test(random1, random23, alternative = "less")
  test_res
}
