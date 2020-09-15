# vector with names of packages required

packages = c(("tidyverse"),
             ("Metrics"),
             ("caret"),
             ("MuMIn"),
             ("lubridate"),
             ("mgcv"),
             ("parallel"),
             ("brms"),
             ("rlist"),
             ("e1071"),
             ("ranger"),
             ("tidymv"),
             ("cowplot"),
             ("rms"))


## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)