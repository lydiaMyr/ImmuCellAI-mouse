## code to prepare `DATASET` dataset goes here

load("data-raw/l1_marker.Rda")
load("data-raw/l2_marker.Rda")
load("data-raw/l3_marker.Rda")
load("data-raw/marker_exp.Rda")
load("data-raw/l1_cell_correction_matrix_new.Rda")
load("data-raw/l2_cell_correction_matrix_new.Rda")
load("data-raw/l3_cell_correction_matrix_new.Rda")


usethis::use_data(l1_marker,overwrite = T)
usethis::use_data(l2_marker,overwrite = T)
usethis::use_data(l3_marker,overwrite = T)
usethis::use_data(marker_exp,overwrite = T)
usethis::use_data(l1_cell_correction_matrix_new,overwrite = T)
usethis::use_data(l2_cell_correction_matrix_new,overwrite = T)
usethis::use_data(l3_cell_correction_matrix_new,overwrite = T)
