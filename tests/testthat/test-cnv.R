context("inferring CNV")


test_that("CNV must be of class cnv", {
    data(dataMLPA)
    CNV  <-  cnv(x  =  dataMLPA$Gene2,  threshold.0  =  0.01,  
               mix.method  =  "mixdist")  
    expect_equal(class(CNV), "cnv")
})
