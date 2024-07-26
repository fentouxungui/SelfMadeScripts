library(memes)
library(magrittr)
options(meme_db = system.file("extdata/flyFactorSurvey_cleaned.meme", package = "memes", mustWork = TRUE))
library(universalmotif)
view_motifs(example_tomtom$best_match_motif)
data("example_tomtom")

flyFactorDb <- MotifDb::MotifDb %>%
  MotifDb::query("FlyFactorSurvey")


MotifDb::query(MotifDb::MotifDb, andStrings = "SANGER", orStrings = c("pad_SANGER_5","tj_SANGER_5",
                                                                      "Hnf4_SANGER_5", "usp_SANGER_5",
                                                                      "sens_SANGER_10","Eip78C_SANGER_5",
                                                                      "Hr83_SANGER_5","cwo_SANGER_5",
                                                                      "Hr78_SANGER_5","phol_SANGER_5")) %>%
  view_motifs( names.pos = "right", text.size = 20)

