#' A score-based DIF test using the parametric bootstrap approach.
#'
#' \code{bootstrap_sctest} computes a score test to detect DIF in multiple 
#' item/parameters with respect to multiple person covariates (\code{DIF_covariate}).
#' A parametric bootstrap approach is applied to obtain p-values. That is, given the 
#' (item and person) parameters, new data sets are sampled to create the distribution 
#' of the test statistic under the null hypothesis. The functionality is limited to 
#' the 1-, 2-, and 3-parameter logistic models.
#' Only DIF with respect to the \code{a} and \code{b} parameters is tested for,
#' which correspond to the item discrimination and the item difficulty parameters.
#'
#' Author: Dries Debeer
#'
#' @param resp A matrix (or data frame) containing the responses, with the
#' items in the columns.
#' @param a A vector of item slopes/item discriminations.
#' @param b A vector of item locations/item difficulties.
#' @param c A vector of pseudo guessing parameters.
#' @param theta A vector with the true/estimated ability parameters or NULL
#' (the default) which leads to the ability parameters being estimated.
#' @param DIF_covariate A list with the person covariate(s) to test for as
#' element(s).
#' @param parameters A character string, either "per_item", "ab", "a", or "b",
#' to specify which parameters should be tested for.
#' @param item_selection A character vector with the column names or an integer
#' vector with the column numbers in the \code{resp}, specifying the items for
#' which the test should be computed. When set to NULL (i.t., the default),
#' all the items are tested.
#' @param nSamples An integer value with the number of permutations to be
#' sampled.
#' @param theta_method A character string, either "wle", "mle", "eap", of
#' "map" that specifies the estimator for the ability estimation. Only
#' relevant when \code{theta == NULL}.
#' @param slope_intercept A logical value indicating whether the slope-intercept
#' formulation of the 2-/3-PL model should be used.
#' @param statistic A character string, either "auto", "DM", "CvM",
#' "maxLM", "LMuo", "WDMo", or "maxLMo", specifying the test statistic to be used.
#' @param meanCenter A logical value: should the score contributions be mean
#' centered per parameter?
#' @param decorrelate A logical value: should the score contributions be
#' decorrelated?
#' @param impact_groups A vector indicating impact-group membership for
#' each person.
#' @return A list with four elements:
#' \describe{
#'   \item{\code{statistics}}{A matrix containing all the test statistics.}
#'   \item{\code{p}}{A matrix containing the obtained \emph{p}-values.}
#'   \item{\code{nSamples}}{The number of samples taken.}
#'   \item{\code{DIF_covariate}}{A list containing all the covariate(s) used to order
#'    the score contributions, as well as the used test statistics.}
#' }
#' @aliases bootstrap_sctest
#' @seealso \code{\link{permutation_sctest}}
#'
#' @examples 
#' \donttest{
#' data("toydata")
#' resp <- toydata$resp
#' group_categ <- toydata$group_categ
#' it <- toydata$it
#' discr <- it[,1]
#' diff <- it[,2]
#' 
#' bootstrap_sctest(resp = resp, DIF_covariate = group_categ, a = discr, b = diff, 
#' decorrelate = FALSE)
#' }
#'
#' @export
bootstrap_sctest <- function(resp,
                             theta = NULL,
                             a = rep(1, length(b)),
                             b,
                             c = rep(0, length(b)),
                             DIF_covariate = NULL,
                             parameters = c("per_item", "ab", "a", "b"),
                             item_selection = NULL,
                             nSamples = 1000,
                             theta_method = c("wle", "mle", "eap", "map"),
                             slope_intercept = FALSE,
                             statistic = "auto",
                             meanCenter = TRUE,
                             decorrelate = FALSE,
                             impact_groups = rep(1, dim(resp)[1])){


  # get call
  call <- match.call()

  # The responses should be in a matrix
  stopifnot(is.matrix(resp) | is.data.frame(resp))
  if(is.data.frame(resp)) resp <- as.matrix(resp)

  # number of persons
  nPerson <- nrow(resp)

  # retrieve theta (or estimate when theta == NULL)
  theta <- get_theta(resp, a, b, c, theta, theta_method, slope_intercept)

  # get list of which_col
  which_col <- get_which_col(item_selection, resp,
                             parameters = match.arg(parameters))

  # create index- matrix according to the DIF_covariates
  index_list <- get_index_list(DIF_covariate, nPerson, statistic, call)

  # get the scores, as well as the terms to compute the scores
  scores_terms <- get_scores(resp, a, b, c, theta,
                             slope_intercept, sparse = FALSE,
                             return_terms = TRUE)

  # scale generated score contributions, rather than scaling the brownian process
  scaled_scores <- scale_scores(scores_terms$scores, meanCenter, decorrelate,
                                impact_groups)

  # compute the test statistic based on the observed scores
  test_stats <- get_stats(scaled_scores, index_list, which_col)

  # get test statistic distribution based on
  bootstrapped_stats <- get_bootstrapped_stats(observed_resp = resp,
                                             terms = scores_terms$terms,
                                             meanCenter, decorrelate,
                                             impact_groups,
                                             index_list, which_col,
                                             nSamples)

  # compute the p-values
  p <- get_pvalues(test_stats, bootstrapped_stats)


  return(list(resp = resp,
              statistic = test_stats,
              p = p,
              nSamples = nSamples,
              DIF_covariate = index_list,
              theta = theta))


}


# function to compute the bootstrapped statistics
get_bootstrapped_stats <- function(observed_resp, terms, meanCenter, decorrelate,
                                  impact_groups, index_list, which_col,
                                  nSamples = 1000
                                  ){

  # cat("Generating bootstrapped samples:\n")
  # progressBar <- txtProgressBar(min = 0, max = nSamples, style = 3, char = "|")
  # env <- environment()

  bootstrapped_stats <- lapply(
    seq_len(nSamples), get_one_bootstrapped_stat, observed_resp, terms,
    meanCenter, decorrelate, impact_groups, index_list, which_col
    #, env
  )

  array(unlist(bootstrapped_stats),
        dim = c(dim(bootstrapped_stats[[1]]), nSamples))
}


# function to compute one bootstrapped statistic
get_one_bootstrapped_stat <- function(sampleNr, observed_resp, terms, meanCenter, decorrelate,
                                      impact_groups, index_list, which_col
                                      #, env
                                      ){

  # setTxtProgressBar(get("progressBar", env), sampleNr)
  # generate response matrix
  gen_resp <- generate_response_matrix(terms$P, observed_resp)

  # compute scores based no generated response matrix
  gen_scores <- get_scores_from_terms(gen_resp, terms)

  # scale generated score contributions, so that they become brownian process
  scaled_gen_scores <- scale_scores(gen_scores, meanCenter, decorrelate, impact_groups)

  # compute statistics
  stats <- get_stats(scaled_gen_scores, index_list, which_col)

  return(stats)
}


# function to generate response matrix based on response probabilities (P)
# with the same pattern of missing responses as the oberved responses
generate_response_matrix <- function(P, observed_resp){
  # generate responses
  dims <- dim(observed_resp)
  generated_resp <- (P > stats::runif(prod(dims))) * 1
  dim(generated_resp) <- dims

  # missing in real data should be missing in generated data
  generated_resp[is.na(observed_resp)] <- NA

  return(generated_resp)
}


# function to get the p-values based on the computed test statistics and
# the bootstrapped/permuted statistics
get_pvalues <- function(test_stats, sampled_stats){
  p <- array(NA, dim = dim(test_stats))

  for(orderNr in seq_len(dim(p)[1])){
    for(which_colNr in seq_len(dim(p)[2])){
      p[orderNr, which_colNr] <- mean(sampled_stats[orderNr, which_colNr, ] >=
                                        test_stats[orderNr, which_colNr])
    }
  }

  dimnames(p) <- dimnames(test_stats)
  return(p)

}

