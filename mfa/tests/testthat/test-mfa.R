context("mfa object")

test_that("mfa creates object of class mfa", {
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  mfa1 <- mfa(df, sets=list(2:7,8:13,14:19,20:23,24:30,31:35,36:39,40:45,46:50,51:54))
  expect_is(mfa1, "mfa")
})

context("check parameters")

test_that("chek_params returns errors for improper inputs", {
  expect_error(check_params(1:5, sets=list(2:7,8:13)))
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  expect_error(check_params(df, sets=FALSE))
  #expect_error(check_params(df, sets=1:(ncol(df)+1)))
})

context("check eigen table")

test_that("eigentbl returns a table with correct number of columns", {
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  mfa1 <- mfa(df, sets=list(2:7,8:13,14:19,20:23,24:30,31:35,36:39,40:45,46:50,51:54))
  expect_equal(ncol(eigentbl(mfa1)), 5)
})

context("check ob2dim")

test_that("ob2dim returns a table with correct number of columns", {
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  mfa1 <- mfa(df, sets=list(2:7,8:13,14:19,20:23,24:30,31:35,36:39,40:45,46:50,51:54))
  expect_equal(ncol(ob2dim(mfa1)), 12)
})

context("check var2dim")

test_that("var2dim returns a table with correct number of columns", {
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  mfa1 <- mfa(df, sets=list(2:7,8:13,14:19,20:23,24:30,31:35,36:39,40:45,46:50,51:54))
  expect_equal(ncol(var2dim(mfa1)), 12)
  expect_is(var2dim(mfa1),"matrix")
})

context("check tbl2dim")

test_that("tbl2dim returns a table with correct class and number of columns", {
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  mfa1 <- mfa(df, sets=list(2:7,8:13,14:19,20:23,24:30,31:35,36:39,40:45,46:50,51:54))
  expect_equal(ncol(tbl2dim(mfa1)), 12)
  expect_is(tbl2dim(mfa1),"matrix")
})

context("check RV")

test_that("RV returns correct value for two columns", {
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  expect_equal(RV(df[,2:7],df[,8:13]), 0.9701797)
})

context("check RVtbl")

test_that("RVtbl returns correct value for two columns", {
  df <- read.csv("../../data/wines.csv", stringsAsFactors = FALSE)
  expect_equal(dim(RVtbl(df, list(2:7,8:13,20:23))), c(3,3))
  expect_equal(diag(RVtbl(df, list(2:7,8:13,20:23))), c(1,1,1))
})
