### Testcode (Updated 22th April 2020)
library(mstDIF)
library(eRm)
library(scDIFtest)
library(mirt)

###  Preparation

### Generate beta, theta, see and response matrix
set.seed(24042020)
i <- 30
p <- 500
beta <- rnorm(i)
theta <- rnorm(p)
group <- factor(c(rep(0, times=p/2), rep(1, times=p/2)))

resp <- eRm::sim.rasch(persons = theta, items = beta)

### For checking: Calculate raw scores. Are there perfect scores?
sum_resp <- apply(resp, 1, sum)

### Estimation of person parameters and standard estimation errors
RM <- eRm::RM(resp)
ppar <- eRm::person.parameter(RM)
theta_est <- unlist(ppar$thetapar)
see <- unlist(ppar$se.theta)

mirt_obj <- mirt::multipleGroup(data = as.data.frame(resp), model=1, itemtype = "2PL", group = group,
                                 invariance = c("free_means","free_var","slopes","intercepts"),
                                 verbose=F)


### Testing
### Test 1: log_reg function
test_that("log_reg is working correctly", {
  test1 <- log_reg(resp=as.data.frame(resp), theta=theta, DIF_covariate = group)
  expect_equal(round(test1$results[1,4],3), 0.646)
  expect_equal(round(test1$results[1,3],3), 0.875)
  expect_equal(round(test1$results[2,4],3), 0.323)
})

### Test 2: mstSIB function
### (Respondents with perfect scores are excluded because eRm provides no estimates for them)
test_that("mstSIB is working correctly", {
  test2 <- mstSIB(resp = as.data.frame(resp[sum_resp > 0 & sum_resp < i,]), theta = theta_est,
                  DIF_covariate = group[sum_resp > 0 & sum_resp < i], see = see)

  expect_equal(round(test2$results$p_value[1],3), 0.956)
  expect_equal(round(test2$results$p_value[2],3), 0.71)
  expect_equal(round(test2$results$p_value[3],3), 0)
})

### Test 3: Score-based permutation tests (calculation with decorrelate = F is somewhat stabler)
# test3 <- permutation_sctest(resp = resp, DIF_covariate = group, a = rep(1,times = i), b = beta, decorrelate = F)

# test_that("Score-based permutation tests are working correctly", {
#  test3 <- permutation_sctest(resp = resp, DIF_covariate = group, a = rep(1,times = i), b = beta, decorrelate = F)

#  expect_equal(round(test3$p[1],3), 0.441)
#  expect_equal(round(test3$p[2],3), 0.643)
#  expect_equal(round(test3$p[3],3), 0.054)
#})

### Test 4: Score-based bootstrap tests (calculation with decorrelate = F is somewhat stabler)
#test4 <- bootstrap_sctest(resp = resp, DIF_covariate = group, a = rep(1,times = i), b = beta, decorrelate = F)

# ### Test 5: Analytical score-based test - needs a mirt object
# test_that("Analytical score-based tests are working correctly", {
#   test5 <- scDIFtest::scDIFtest(mirt_obj, DIF_covariate = group)
#   expect_equal(round(summary(test5)$p_value[1],3), 0.561)
#   expect_equal(round(summary(test5)$p_value[2],3), 0.599)
#   expect_equal(round(summary(test5)$p_value[3],3), 0.168)
# })

### Test 6: mstDIF function with logreg option
test_that("log_reg option of mstDIF is working correctly", {
  test6 <- mstDIF(as.data.frame(resp), DIF_covariate = group, method = "logreg", theta = theta)

  expect_equal(round(test6$method_results$results[1,4],3), 0.646)
  expect_equal(round(test6$method_results$results[1,3],3), 0.875)
  expect_equal(round(test6$method_results$results[2,4],3), 0.323)
})

### Test 7: mstDIF function with mstSIB option
test_that("mstSIB option of mstDIF is working correctly", {
  test7 <- mstDIF(as.data.frame(resp[sum_resp > 0 & sum_resp < i,]),
                  DIF_covariate = group[sum_resp > 0 & sum_resp < i],
                  method = "mstsib", theta = theta_est, see=see)

  expect_equal(round(test7$method_results$results$p_value[1],3), 0.956)
  expect_equal(round(test7$method_results$results$p_value[2],3), 0.71)
  expect_equal(round(test7$method_results$results$p_value[3],3), 0)
})

### Test 8: mstDIF function with bootstrap option
#test8 <- mstDIF(resp = resp, DIF_covariate = group, method = "bootstrap",
#                a = rep(1,times = i), b = beta, decorrelate = F)

### Test 9: mstDIF function with permutation option
#test9 <- mstDIF(resp = resp, DIF_covariate = group, method = "permutation",
#                a = rep(1,times = i), b = beta, decorrelate = F)

### Test 10: mstDIF function with analytical option
test_that("Analytical option of mstDIF is working correctly", {
  test10 <- mstDIF(mirt_obj, DIF_covariate = group, method = "analytical")

  expect_equal(round(summary(test10)$p_value[1],3), 0.007)
  expect_equal(round(summary(test10)$p_value[2],3), 0.016)
  expect_equal(round(summary(test10)$p_value[3],3), 0.052)
})
