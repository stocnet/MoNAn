# goodness of fit: getIndegree and getTieWeights in gofMobilityNetwork

dv <- myResDN$state$dep.var
nets <- c(list(myResDN$cache[[dv]]$valuedNetwork),
          lapply(myResDN$deps[1:20], function(d) {
            s <- myResDN$state
            s[[dv]] <- d
            createInternalCache(s)[[dv]]$valuedNetwork
          }))
ansSmall <- myResDN
ansSmall$deps <- myResDN$deps[1:20]
# level k counts the value k - 1
countIndegree <- function(m, lvls) sapply(lvls, function(k) sum(colSums(m) == k - 1))
countTieWeights <- function(m, lvls) sapply(lvls, function(k) sum(m == k - 1))

test_that("gof functions count the values of the observed and simulated networks", {
  for (lvls in list(1:100, 20:70, c(3, 1, 5))) {
    g <- gofMobilityNetwork(ansSmall, getIndegree, lvls = lvls)
    expect_equal(g$observed, countIndegree(nets[[1]], lvls))
    expect_equal(g$simulated, lapply(nets[-1], countIndegree, lvls = lvls))
  }
  for (lvls in list(1:30, 2:10)) {
    g <- gofMobilityNetwork(ansSmall, getTieWeights, lvls = lvls)
    expect_equal(g$observed, countTieWeights(nets[[1]], lvls))
    expect_equal(g$simulated, lapply(nets[-1], countTieWeights, lvls = lvls))
  }
})

test_that("without lvls, all networks get the same levels", {
  for (f in list(list(getIndegree, countIndegree, function(m) max(colSums(m))),
                 list(getTieWeights, countTieWeights, max))) {
    g <- gofMobilityNetwork(ansSmall, f[[1]])
    lvls <- 1:(max(sapply(nets, f[[3]])) + 1)
    expect_equal(g$observed, f[[2]](nets[[1]], lvls))
    expect_equal(g$simulated, lapply(nets[-1], f[[2]], lvls = lvls))
  }
})

test_that("gof results can be plotted with and without lvls", {
  g <- gofMobilityNetwork(ansSmall, getTieWeights)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(g))
  expect_no_error(plot(g, lvls = 1:5))
})
