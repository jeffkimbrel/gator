library(tidyverse)
library(reticulate)

use_condaenv("base")

jakomics <- import("jakomics")

print(jakomics$version)

# hmmsearch
jakomics$hmm$run_hmmsearch("/Users/kimbrel1/Library/CloudStorage/OneDrive-LLNL/Documents/Biofuels_SFA/ARW/data/nrMAGs/faa/mARW1_16.faa",
                           "test.log",
                           "test_output",
                           "/Users/kimbrel1/Library/CloudStorage/Dropbox/Lab/Resources/gator/phospho.hmm",
                           cut_tc=FALSE,
                           echo=TRUE,
                           run=TRUE)

# blast
l = jakomics$blast$run_blast(type = "prot",
                         q = "/Users/kimbrel1/Library/CloudStorage/OneDrive-LLNL/Documents/Biofuels_SFA/ARW/data/nrMAGs/faa/mARW1_16.faa",
                         db = "/Users/kimbrel1/Library/CloudStorage/Dropbox/Lab/Resources/gator/gator.faa")

jakomics$blast$blast_to_df(l) |>
  as_tibble() |>
  unnest(cols = everything())

