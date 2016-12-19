# Copyright (C) 2016 Lyron Winderbaum
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



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

test_that("combine_peaklists of invalid input produces correct warning", {
  expect_warning(combine_peaklists("."),
                 "dipps::combine_peaklists: No peaklist files found, aborting.")
})

# Cleanup
unlink("*_peaklist.txt")
unlink("*_speclist.txt")
