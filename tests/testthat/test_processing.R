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

test_that("dbscan of sparse points produces appropriate warning", {
  expect_warning(dbscan(1:10),
                 "dipps::dbscan: No core points.")
  expect_warning(dbscan(1:1000),
                 "dipps::dbscan: No core points.")
})


test_that("dbscan of invalid inputs produces appropriate errors", {
  expect_error(dbscan("Hello"),
               "dipps::dbscan: invalid input: x is not a numeric vector.")
  expect_error(dbscan(c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)),
               "dipps::dbscan: invalid input: x is not a numeric vector.")

  expect_error(dbscan(1:10, eps = "Hello"),
               "dipps::dbscan: invalid input: eps is not a numeric vector.")
  expect_error(dbscan(1:10, eps = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)),
               "dipps::dbscan: invalid input: eps is not a numeric vector.")
  expect_error(dbscan(1:10, eps = 1:3),
               "dipps::dbscan: invalid input: multiple eps values provided.")
  expect_error(dbscan(1:10, eps = -3),
               "dipps::dbscan: invalid input: non-positive eps provided.")
  expect_error(dbscan(1:10, eps = 0),
               "dipps::dbscan: invalid input: non-positive eps provided.")

  expect_error(dbscan(1:10, mnpts = "Hello"),
               "dipps::dbscan: invalid input: mnpts is not a numeric vector.")
  expect_error(dbscan(1:10, mnpts = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE)),
               "dipps::dbscan: invalid input: mnpts is not a numeric vector.")
  expect_error(dbscan(1:10, mnpts = 1:3),
               "dipps::dbscan: invalid input: multiple mnpts values provided.")
  expect_error(dbscan(1:10, mnpts = -3),
               "dipps::dbscan: invalid input: non-positive mnpts provided.")
  expect_error(dbscan(1:10, mnpts = 0),
               "dipps::dbscan: invalid input: non-positive mnpts provided.")

  expect_error(dbscan(1:10, pp = "Hello"),
               "dipps::dbscan: invalid input: pp is not a logical vector.")
  expect_error(dbscan(1:10, pp = 1:3),
               "dipps::dbscan: invalid input: pp is not a logical vector.")
  expect_error(dbscan(1:10, pp = c(TRUE, FALSE, TRUE)),
               "dipps::dbscan: invalid input: multiple pp values provided.")
})
