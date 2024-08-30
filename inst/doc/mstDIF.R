## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(mstDIF)

## -----------------------------------------------------------------------------
data("toydata")

## -----------------------------------------------------------------------------
resp <- toydata$resp
group_categ <- toydata$group_categ
group_cont <- toydata$group_cont
it <- toydata$it
theta_est <- toydata$theta_est
see_est <- toydata$see_est

## -----------------------------------------------------------------------------
log_reg_DIF <- mstDIF(resp, DIF_covariate = factor(group_categ), method = "logreg",
                theta = theta_est)

## -----------------------------------------------------------------------------
log_reg_DIF

## -----------------------------------------------------------------------------
summary(log_reg_DIF, DIF_type = "all")

## -----------------------------------------------------------------------------
mstSIB_DIF <- mstDIF(resp, DIF_covariate = factor(group_categ), method = "mstsib",
                theta = theta_est, see = see_est)
mstSIB_DIF

## -----------------------------------------------------------------------------
summary(mstSIB_DIF)

## -----------------------------------------------------------------------------
library(mirt)
mirt_model <- mirt(as.data.frame(resp), model = 1, verbose = FALSE)

## -----------------------------------------------------------------------------
sc_DIF <- mstDIF(mirt_model, DIF_covariate = factor(group_categ), method = "analytical")
sc_DIF

## -----------------------------------------------------------------------------
summary(sc_DIF)

## -----------------------------------------------------------------------------
sc_DIF_2 <- mstDIF(mirt_model, DIF_covariate = group_cont, method = "analytical")
sc_DIF_2

## -----------------------------------------------------------------------------
summary(sc_DIF_2)

## -----------------------------------------------------------------------------
discr <- it[,1]
diff <- it[,2]

## -----------------------------------------------------------------------------
bootstrap_DIF <- mstDIF(resp = resp, DIF_covariate = group_categ, method = "bootstrap",
                a = discr, b = diff, decorrelate = F)

## -----------------------------------------------------------------------------
bootstrap_DIF

## -----------------------------------------------------------------------------
summary(bootstrap_DIF)

## -----------------------------------------------------------------------------
bootstrap_DIF_2 <- mstDIF(resp = resp, DIF_covariate = group_cont, method = "bootstrap",
                a = discr, b = diff, decorrelate = F)
bootstrap_DIF_2

## -----------------------------------------------------------------------------
summary(bootstrap_DIF_2)

## -----------------------------------------------------------------------------
permutation_DIF <- mstDIF(resp = resp, DIF_covariate = group_categ, method = "permutation",
                a = discr, b = diff, decorrelate = F)
permutation_DIF_2 <- mstDIF(resp = resp, DIF_covariate = group_cont, method = "permutation",
                a = discr, b = diff, decorrelate = F)

## -----------------------------------------------------------------------------
summary(permutation_DIF)

## -----------------------------------------------------------------------------
summary(permutation_DIF_2)

