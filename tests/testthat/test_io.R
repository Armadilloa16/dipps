
test_that("combine_peaklists returns correct number of empty spectra.", {
  expect_equal(combine_peaklists(system.file("extdata", "test1",
                                             package = "dipps")),
               0)
  expect_equal(combine_peaklists(system.file("extdata", "test2",
                                             package = "dipps"),
                                 i.name = "."),
               0)
  expect_equal(combine_peaklists(system.file("extdata", "test3",
                                             package = "dipps")),
               2)
})

# Cleanup
unlink("*_peaklist.csv")
unlink("*_speclist.csv")

test_that("combine_peaklists of invalid input produces correct warning", {
  expect_warning(combine_peaklists("."),
                 "dipps::combine_peaklists: No peaklist files found, aborting.")
})
