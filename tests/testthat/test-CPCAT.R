test_that("hypotheses function returns a list", {
  expect_type(hypotheses(3),"list")
})

test_that("hypotheses function returns where list with each list entry is of type double", {
  expect_type(hypotheses(3)[[1]],"double")
})

test_that("hypotheses function returns a pre-specified output", {
  expect_equal(hypotheses(3)[[1]],matrix(c(1,0,0,1,1,0,1,0,1,1,1,1),byrow=TRUE,nrow=4),ignore_attr=TRUE)
})

test_that("poisson.sub.test function returns a type double value", {
  dathelp=list()
  dathelp[[1]]=testdata2[1:4,]
  dathelp[[2]]=testdata2[5:8,]
  dathelp[[3]]=testdata2[9:12,]
  dathelp[[4]]=testdata2[13:16,]
  expect_type(poisson.sub.test(dat=dathelp, contrast=matrix(c(1,0,0),nrow=1)),"double")
})

test_that("poisson.test function returns a type double value  (output part 1)", {
  expect_type(poisson.test(testdata,matrix(c(1,0,0,1,1,0,1,0,1,1,1,1),byrow=TRUE,nrow=4))[[1]],"double")
})

test_that("poisson.test function returns a type double value (output part 2)", {
  expect_type(poisson.test(testdata,matrix(c(1,0,0,1,1,0,1,0,1,1,1,1),byrow=TRUE,nrow=4))[[2]],"double")
})

test_that("CPCAT function returns a type double value", {
  dathelp=data.frame(Poissondata=c(1:9),Group=factor(rep(c("A","B","C"),each=3)))
  expect_type(CPCAT(testdata),"list")
})
