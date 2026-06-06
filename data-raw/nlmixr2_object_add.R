library(usethis)
nlmixr2Keywords <- read.csv("nlmixr2_object.csv");
nlmixr2Keywords <- nlmixr2Keywords[order(nlmixr2Keywords$Note, nlmixr2Keywords$Field), ]
rownames(nlmixr2Keywords) <-NULL
usethis::use_data(nlmixr2Keywords, overwrite=TRUE)
