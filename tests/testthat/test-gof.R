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
# level k counts the value k - 1; the counts are named by the values
countIndegree <- function(m, lvls) {
  setNames(sapply(lvls, function(k) sum(colSums(m) == k - 1)), lvls - 1)
}
countTieWeights <- function(m, lvls) {
  setNames(sapply(lvls, function(k) sum(m == k - 1)), lvls - 1)
}

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

test_that("gof results can be plotted with and without lvls, labelled by value", {
  g <- gofMobilityNetwork(ansSmall, getIndegree, lvls = 20:70)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  p <- plot(g, lvls = 1:5)
  expect_equal(p$names, as.character(19:23))
  expect_no_error(plot(gofMobilityNetwork(ansSmall, getTieWeights)))
  expect_no_error(plot(g, xlab = "indegree", main = "Indegree distribution"))
  # results of earlier versions (unnamed) are plotted by position
  unnamed <- g
  unnamed$observed <- unname(g$observed)
  unnamed$simulated <- lapply(g$simulated, unname)
  expect_equal(plot(unnamed, lvls = 1:3)$names, c("1", "2", "3"))
})
