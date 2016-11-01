
test_that("dbscan of invalid input throws correct error", {
  expect_error(dbscan(data.frame(mz = 1:10)),
               "dipps::dbscan: Non-existent variable selected for clustering.")
  expect_error(dbscan(data.frame(mz = 1:1000)),
               "dipps::dbscan: Non-existent variable selected for clustering.")
})

test_that("dbscan of sparse points produces correct warning", {
  expect_warning(dbscan(data.frame(m.z = 1:10)),
                 "dipps::dbscan: No core points.")
  expect_warning(dbscan(data.frame(m.z = 1:1000)),
                 "dipps::dbscan: No core points.")
})
