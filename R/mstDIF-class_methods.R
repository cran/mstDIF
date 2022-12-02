#' Methods for the mstDIF-class
#'
#' \code{print} and \code{summary} methods for objects of the
#' \code{mstDIF-class}, as returned by \code{\link{mstDIF}}. See details
#' for more information about the methods.
#'
#' The \code{print} method prints some basic information about the
#'   \code{mstDIF-class} object.
#'
#' The \code{summary} method computes a data frame with a row for each item
#' that was included in the test. The columns are:
#'   \describe{
#'      \item{item}{The name of the item}
#'      \item{statistic}{The value for the used statistic per item}
#'      \item{p_value}{The p-value per item}
#'      \item{eff_size}{An effect-size for the DIF-test, if applicable}
#'    }
#'
#'
#' @param x an object of class \code{mstDIF}
#' @param object an object of class \code{mstDIF}
#' @param DIF_type a string that should one or more of "overall", "uniform",
#'    "non-uniform", "all".
#' @param ordered logical: should the summary be ordered according to the obtained p-values (in ascending order)?
#' @param ... other arguments passed to the method.
#'
#' @examples
#'
#' ## load data
#' data("toydata")
#'
#' ## fit 2PL model using mirt
#' mirt_model <- mirt::mirt(toydata$resp, model = 1)
#'
#' ## test DIF along a contiuous covariate
#' DIFtest <- mstDIF(mirt_model, DIF_covariate = toydata$group_cont,
#' method = "analytical")
#'
#' ## print
#' DIFtest
#'
#' ## summary
#' summary(DIFtest)
#'
#'
#' @name mstDIF-Methods
NULL
#' @rdname mstDIF-Methods
#' @aliases print.mstDIF
#' @export
print.mstDIF <- function(x, ...){
  cat("    Differential Item Functioning (DIF) Detection Test \n",
      "\tMethod: \t", x$method, "\n",
      "\tTest: \t \t", x$test, "\n",
      "\tDIF covariate: \t", as.character(deparse(x$call$DIF_covariate)), "\n",
      "\tData: \t \t", as.character(deparse(x$call$resp)), "\n",
      "\tItems: \t \t", dim(x$resp)[2], "\n",
      "\tPersons: \t", dim(x$resp)[1])

}

#' @rdname mstDIF-Methods
#' @aliases summary.mstDIF
#' @export
summary.mstDIF <- function(object, DIF_type = "overall",
                           ordered = TRUE, ...){

  # match DIF_type with arguments
  DIF_type_choices <- c("overall", "uniform",
                        "non-uniform", "all")
  DIF_type <- match.arg(DIF_type, DIF_type_choices, several.ok = TRUE)

  if("all" %in% DIF_type) DIF_type <- names(object$DIF_test)

  # is DIF_type included in the object?
  in_object <- DIF_type %in% names(object$DIF_test)

  if(!all(in_object))
    warning(paste0(paste(DIF_type[!in_object], collapse = " and "),
                   " DIF tests are not included in the ",
                   object$method, "."),
            call. = FALSE)

  col_names <- c("stat", "p_value")
  if(object$method == "DIF-test using Logistic Regression") col_names <- c(col_names, "eff_size")

  summary <- do.call(cbind, lapply(DIF_type[in_object], function(name){
    out <- object$DIF_test[[name]][col_names]
    if(length(DIF_type) > 1) names(out) <- paste(name, names(out), sep = "_")
    out
  }))

  rownames(summary) <- colnames(object$resp)

  # include number of observations per test if applicable
  if(!is.null(object$DIF_test$overall$N)) summary$N <- object$DIF_test$overall$N

  # order summary according to the p_values?
  if(ordered) summary <- summary[
    do.call(order,
            as.list(summary[,grep(pattern = "p_value", names(summary)),
                            drop = FALSE])), ]
  summary
}


