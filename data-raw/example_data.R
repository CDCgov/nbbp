# https://en.wikipedia.org/wiki/Borealpox_virus, accessed 2024-09-25
# "As of February 2024, there are seven reported cases of illness"
borealpox <- rep(1, 7)

# Table 1 of https://doi.org/10.1016/j.jtbi.2011.10.039
pneumonic_plague <- c(8, 5, 4, 5009, 2, 13, 35, 18, 39, 16, 42, 3, 18, 12, 30, 10, 2, 2, 6)

# Text S2 of https://doi.org/10.1371/journal.ppat.1004452
measles_us_97 <- c(
  rep(1, 122),
  rep(2, 13),
  rep(3, 10),
  rep(4, 6),
  rep(5, 5),
  rep(6, 2),
  rep(8, 2),
  rep(9, 1),
  rep(11, 1),
  rep(13, 1),
  rep(15, 1),
  rep(33, 1)
)

measles_canada_98 <- c(
  rep(1, 35),
  rep(2, 5),
  rep(3, 3),
  rep(4, 1),
  rep(6, 1),
  rep(8, 1),
  rep(17, 1),
  rep(30, 1),
  rep(155, 1)
)

# apparent typo removed: duplicated row for 7
mers_pre_june <- c(
  rep(1, 11),
  rep(2, 1),
  rep(3, 3),
  rep(4, 2),
  rep(5, 1),
  rep(7, 1),
  rep(26, 1)
)

mers_post_june <- c(
  rep(1, 16),
  rep(3, 1),
  rep(4, 1),
  rep(5, 1)
)

usethis::use_data(borealpox, overwrite = TRUE)
usethis::use_data(pneumonic_plague, overwrite = TRUE)
usethis::use_data(measles_us_97, overwrite = TRUE)
usethis::use_data(measles_canada_98, overwrite = TRUE)
usethis::use_data(mers_pre_june, overwrite = TRUE)
usethis::use_data(mers_post_june, overwrite = TRUE)
