#' A logistic regression DIF test for MSTs
#'
#' This function allows the detection of itemwise DIF for Multistage Tests. It is based on the comparison
#' of three logistic regression models for each item. The first logistic regression model (Model 1)
#' predicts the positiveness of each response solely on the estimated ability parameters. The second
#' logistic regression model (Model 2) predicts the positiveness based on the ability parameters
#' and the membership to the focal and reference group as additive predictor variables.
#' The third model (Model 3) uses the same predictors as Model 2 to predict the positiveness of the responses, but
#' also includes an interaction effect. Three model comparisons are carried out (Models 1/2, Models 1/3, Models 2/3)
#' based on two criteria: The comparison of the Nagelkerke R squared values, and the p-values of a likelihood ratio test.
#'
#' Author: Sebastian Appelbaum, with minor changes by Rudolf Debelak and Dries Debeer
#'
#' @param resp A data frame containing the response matrix. Rows correspond to respondents, columns to items.
#' @param DIF_covariate A factor indicating the membership to the reference and focal groups.
#' @param theta A vector of ability estimates for each respondent.
#'
#' @return A list with four elements. The first element is the response matrix, the second element is the name of
#' the DIF covariate, and the third element is the name of the test. The fourth element is a data frame where
#' each row corresponds to an item. The columns of this data frame correspond to the following entries:
#' \describe{
#'   \item{\code{N}}{The number of responses observed for this item.}
#'   \item{\code{overall_chi_sq}}{The chi squared statistic of the likelihood ratio test comparing Model 1 and Model 3.}
#'   \item{\code{overall_p_value}}{The p-values of the likelihood ratio test comparing Model 1 and Model 3 as
#'   an indicator for the overall DIF effect.}
#'   \item{\code{Delta_NagelkerkeR2}}{The difference of the Nagelkerke R squared values for Model 1 and Model 3.}
#'   \item{\code{UDIF_chi_sq}}{The chi squared statistic of the likelihood ratio test comparing Model 1 and Model 2.}
#'   \item{\code{UDIF_p_value}}{The p-values of the likelihood ratio test comparing Model 1 and Model 2.}
#'   \item{\code{UDIF_Delta_NagelkerkeR2}}{The difference of the Nagelkerke R squared values for Model 1 and Model 2.}
#'   \item{\code{CDIF_chi_sq}}{The chi squared statistic of the likelihood ratio test comparing Model 2 and Model 3.}
#'   \item{\code{CDIF_p_value}}{The p-values of the likelihood ratio test comparing Model 2 and Model 3.}
#'   \item{\code{CDIF_Delta_NagelkerkeR2}}{The difference of the Nagelkerke R squared values for Model 2 and Model 3.}
#' }
#'
#' @examples 
#' data("toydata")
#' resp <- toydata$resp
#' group_categ <- toydata$group_categ
#' theta_est <- toydata$theta_est
#' log_reg(resp, DIF_covariate = factor(group_categ), theta = theta_est)
#'
#'
#' @export
log_reg <- function(resp, DIF_covariate, theta = NULL){

  # get call
  call <- match.call()

  # a theta-argument is required
  if(is.null(theta)) stop("'theta'-argument is missing. Include a vector with the estimated 'theta'-values.", call. = FALSE)

  # get the DIF_covariate name
  DIF_covariate_name <- as.character(deparse(call$DIF_covariate))

  if(!is.factor(DIF_covariate)){
    DIF_covariate <- as.factor(DIF_covariate)
    message(paste0("For the DIF analysis, '", DIF_covariate_name,
                   "' was transformed into a factor."))
  }

  # Helper functions
  R2 <- function(m, n) 1 - (exp(-m$null.deviance/2 + m$deviance/2))^(2/n)
  R2max <- function(m, n) 1 - (exp(-m$null.deviance/2))^(2/n)
  R2DIF <- function(m, n) R2(m, n)/R2max(m, n)              # NagelkerkeR2
  R2DIF2 <- function(m1, m2, n) R2DIF(m2, n) - R2DIF(m1, n) # Delta_NagelkerkeR2

  # function to compute the log-reg-DIF-test for one item
  log_reg_DIF_1_item <- function(itemname, data){

    glm_2 <- stats::glm(
      stats::as.formula(paste0(itemname," ~ theta * DIF_covariate")),
      family = stats::binomial, data = d)

    # number of observations in test
    N <- summary(glm_2)$df.null + 1


    # does it make sense to test for DIF?
    if (dim(summary(glm_2)$coef)[1] == 4){
      glm_1 <- stats::glm(
        stats::as.formula(paste0(itemname," ~ theta + DIF_covariate")),
        family = stats::binomial, data = d)
      glm_0 <- stats::glm(
        stats::as.formula(paste0(itemname," ~ theta")),
        family = stats::binomial, data = d)
      test_udif <- stats::anova(glm_0, glm_1, test = "LRT")
      test_cdif <- stats::anova(glm_1, glm_2, test = "LRT")
      test_dif <- stats::anova(glm_0, glm_2, test = "LRT")

      return(data.frame(item = itemname,
                        N = N,
                        stat = test_dif$Deviance[2],
                        p_value = test_dif$`Pr(>Chi)`[2],
                        eff_size = R2DIF2(glm_0, glm_2, N),
                        stat_u = test_udif$Deviance[2],
                        p_value_u = test_udif$`Pr(>Chi)`[2],
                        eff_size_u = R2DIF2(glm_0, glm_1, N),
                        stat_nu = test_cdif$Deviance[2],
                        p_value_nu = test_cdif$`Pr(>Chi)`[2],
                        eff_size_nu = R2DIF2(glm_1, glm_2, N),
                        stringsAsFactors = FALSE))

    } else {
      return(data.frame(item = itemname,
                        N = N,
                        stat = NA,
                        p_value = NA,
                        eff_size = NA,
                        stat_u = NA,
                        p_value_u = NA,
                        eff_size_u = NA,
                        stat_nu = NA,
                        p_value_nu = NA,
                        eff_size_nu = NA,
                        stringsAsFactors = FALSE))
    }
  }

  # get number of items
  nItem <- ncol(resp)

  # get/set item names
  colnames(resp) <- itemnames <- 'if'(is.null(colnames(resp)),
                                      sprintf(paste("it%0", nchar(nItem),
                                                    "d", sep=''),
                                              seq_len(nItem)),
                                      colnames(resp))

  # combine estimated theta's with grouping variable and responses
  d <- as.data.frame(cbind(theta, DIF_covariate, resp))

  # apply the test for one item to all the items
  res <- do.call(rbind, lapply(itemnames, log_reg_DIF_1_item, data = d))
  return(list(resp = resp,
              DIF_covariate = DIF_covariate_name,
              test = "Likelihood Ratio Test",
              results = res))
}
