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


test_that("dipps of invalid inputs produces appropriate errors", {
  expect_error(dipps("Hello", 1:10, 1:10 < 6),
               "dipps::dipps: invalid input: obs is not a numeric vector.")
  expect_error(dipps(1:10 < 6, 1:10, 1:10 < 6),
               "dipps::dipps: invalid input: obs is not a numeric vector.")
  expect_error(dipps(1:10, "Hello", 1:10 < 6),
               "dipps::dipps: invalid input: var is not a numeric vector.")
  expect_error(dipps(1:10, 1:10 < 6, 1:10 < 6),
               "dipps::dipps: invalid input: var is not a numeric vector.")
  expect_error(dipps(1:10, 1:10, "Hello"),
               "dipps::dipps: invalid input: sub is not a logical vector.")
  expect_error(dipps(1:10, 1:10, 1:10),
               "dipps::dipps: invalid input: sub is not a logical vector.")
  
  expect_error(dipps(1:11, 1:10, 1:10 < 6),
               "dipps::dipps: invalid inputs: inputs are unequal lengths.")
  expect_error(dipps(1:10, 1:11, 1:10 < 6),
               "dipps::dipps: invalid inputs: inputs are unequal lengths.")
  expect_error(dipps(1:10, 1:10, 1:11 < 6),
               "dipps::dipps: invalid inputs: inputs are unequal lengths.")

  expect_error(dipps(rep(1:5, 2), rep(1:5, 2), rep(1:5, 2) < 3),
               paste("dipps::dipps: invalid inputs: occurences not represented",
                                                   "uniquely."))
})

