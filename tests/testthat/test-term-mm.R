#  File tests/testthat/test-term-mm.R in package ergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

# a bipartite nw
set.seed(143)
b1 <- floor(runif(60, 1,100))
b2 <- floor(runif(60, 101, 130))
exbip.el <- cbind(b1,b2)
bipnw <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw %v% "Letter" <- letters[1:3]
bipnw %v% "Cost" <- c(3,2,1)
                          

# another bipartite nw with more ties and 2 attributes
set.seed(258)
b1 <- floor(runif(150, 1,200))
b2 <- floor(runif(150, 201, 400))
exbip.el <- cbind(b1,b2)
bipnw2 <- as.network(exbip.el, matrix.type="edgelist", bipartite=100, directed=FALSE)
bipnw2 %v% "Letter" <- letters[1:2]
color <- rbinom(400, 1, .4)
color[color ==1] <- "Purple"
color[color ==0] <- "Gold"
bipnw2 %v% "Color" <- color


# a directed nw
data(sampson)
set.seed(42)
set.edge.attribute(samplike, "YearsTrusted", rbinom(88, 4, .5))
set.seed(296)
set.vertex.attribute(samplike, "YearsServed", rbinom(18, 10, .5))
samplike %v% "Trinity" <- c("F", "S", "H")


# an undirected nw
data(faux.mesa.high)
fmh <- faux.mesa.high
set.seed(7)
set.edge.attribute(fmh, "GradeMet", rbinom(203, 6, .5))


# a small undirected nw w/ lots o' triangles
set.seed(20)
t<-trunc(runif(160, 1, 20))
set.seed(21)
h<-trunc(runif(160, 1, 20))
el <- cbind(t,h)
bad <- which(el[,2]==el[,1])
el[bad,2] = el[bad,2]+1
unnw <- network(el, directed=FALSE)
unnw %v% "Pet" <- c("dog", "cat")

test_that("Undirected mm() summary", {
  s.a <- summary(fmh ~ mm("Grade"))
  expect_equivalent(s.a, c(0, 33, 0, 2, 23, 1, 4, 7, 9, 1,
                           2, 6, 1, 17, 1, 1, 4, 5, 5, 6))
})

test_that("Directed mm() ERGM", {
  e.a <- ergm(samplike ~ mm("group", levels2=TRUE))
  expect_equivalent(coef(e.a),
                    c(0.191055236762708, -3.29583686600359,
                      -2.17475172140943, -2.5649493574587,
                      1.6094379124341, -3.29583686600359,
                      -1.49165487677766, -1.09861228866811,
                      0.916290731874155))
})

test_that("Bipartite mm() summary", {
  s.ab <- summary(bipnw ~ mm("Letter"))
  expect_equivalent(s.ab, c(8,8,7,7,5,4,6,6))
})


test_that("Bipartite mm() ERGM with level2 filter", {
  e.ab <- ergm(bipnw ~ mm("Letter", levels2=-(2:6)))
  expect_equivalent(coef(e.ab),
                    c(-3.49650756126405, -4.43081679479004,
                      -3.98898404450688, -3.98898404450688 ))
})

test_that("Undirected mm() summary with level2 filter", {
  s.ab2 <- summary(fmh ~ mm("Race", levels2=-1))
  expect_equivalent(s.ab2, c(8,53,13,41,46,0,1,0,0,5,22,10,0,4))
})

test_that("Undirected mm() summary with level2 filter by logical matrix", {
  M <- matrix(TRUE, 5, 5)
  M[1,1] <- M[3,2] <- M[2,3] <- FALSE
  s.ab2 <- summary(fmh ~ mm("Race", levels2=M))
  expect_equivalent(s.ab2, c(0,8,53,13,41,46,0,1,0,0,5,22,10,0,4)[M[upper.tri(M,TRUE)]])
})

test_that("Undirected mm() summary with level2 filter by numeric matrix", {
  M <- cbind(2,3)
  s.ab2 <- summary(fmh ~ mm("Race", levels2=M))
  expect_equivalent(s.ab2, c(0,8,53,13,41,46,0,1,0,0,5,22,10,0,4)[5])
})

test_that("Directed mm() ERGM with level2 filter", {
  e.ab2 <- ergm(samplike ~ mm("Trinity", levels2=-(3:9)))
  expect_equivalent(coef(e.ab2),
                    c(-1.01160091056776, -0.693147180549145))
})

test_that("Directed mm() ERGM with level2 filter", {
  e.ab2 <- ergm(samplike ~ mm("Trinity", levels2=-(3:9)))
  expect_equivalent(coef(e.ab2),
                    c(-1.01160091056776, -0.693147180549145))
})

test_that("Directed mm() summary with level2 matrix filter", {
  M <- matrix(FALSE, 3, 3)
  M[1,2] <- M[1,3] <- TRUE
  s.ab2 <- summary(samplike ~ mm("Trinity", levels2=M))
  expect_equivalent(s.ab2, c(8,13))
})

test_that("Undirected mm() marginal summary", {
  s.a <- summary(fmh ~ mm(.~Grade))
  expect_equivalent(s.a, c(75, 65, 36, 49, 28))
})

test_that("Undirected mm() asymmetric two-sided summary", {
  s.a.tab <- as.matrix(fmh, matrix.type="edgelist") %>%
    rbind(.,.[,2:1,drop=FALSE]) %>%
    (function(x) data.frame(Race=factor(fmh%v%"Race")[x[,1]], Grade=factor(fmh%v%"Grade")[x[,2]])) %>%
    table()
  s.a <- summary(fmh ~ mm(Race~Grade, levels2=TRUE))
  expect_equivalent(s.a, c(s.a.tab))
})

test_that("Undirected mm() asymmetric two-sided summary and pipe operator on one side", {
  Grade <- fmh%v%"Grade" %>% factor()
  Race <- fmh%v%"Race" %>% replace(., . %in% c("Black","White","Other"), "BWO") %>% factor()

  s.a.tab <- as.matrix(fmh, matrix.type="edgelist") %>%
    rbind(.,.[,2:1,drop=FALSE]) %>%
    (function(x) data.frame(Grade=Grade[x[,2]], Race=Race[x[,1]])) %>%
    table()
  s.a <- summary(fmh ~ mm(Grade~(~Race) %>% COLLAPSE_SMALLEST(3,"BWO"), levels2=TRUE))
  expect_equivalent(s.a, c(s.a.tab))
})

test_that("Undirected mm() marginal summary with fixed levels set", {
  s.a <- summary(fmh ~ mm(.~Grade, levels=I(c(7,6,9,8)), levels2=TRUE))
  expect_equivalent(s.a, c(153, 0, 65, 75))
})

test_that("Undirected mm() summary with fixed levels set", {
  s.a <- summary(fmh ~ mm("Grade", levels=I(c(7,6,9,8)), levels2=TRUE))
  expect_equivalent(s.a, c(75, 0, 0, 0, 0, 23, 0, 0, 2, 33))
})


test_that("Undirected valued mm() sum summary", {
  s.a <- summary(fmh ~ mm("Grade"), response="GradeMet")
  expect_equivalent(s.a,
                    c(0, 96, 0, 6, 75, 3, 15, 21, 24, 3, 6, 19,
                      4, 54, 2, 5, 11, 17, 18, 16))
})

test_that("Undirected valued mm() nonzero summary", {
  s.a <- summary(fmh ~ mm("Grade", form="nonzero"), response="GradeMet")
  expect_equivalent(s.a,
                    c(0, 33, 0, 2, 22, 1, 4, 7, 8, 1, 2, 6, 1, 17,
                      1, 1, 4, 5,  5, 6))
})

