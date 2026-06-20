# Build the bundled example data sets from the original Le Roux et al. scripts.
# Run from the package root with:  Rscript data-raw/prepare_data.R
src <- "LeRoux Scripts/Data"

read_le_roux <- function(file) {
  read.table(file.path(src, file), header = TRUE, sep = ";", dec = ".",
             row.names = 1)
}

Target       <- read_le_roux("Target.txt")          # 10 x 2 reference cloud
Target_group <- read_le_roux("Target_group.txt")    # 4 x 2 group cloud
Student      <- read_le_roux("Student.txt")         # 10 x 1 (one-dimensional)

dir.create("data", showWarnings = FALSE)
save(Target,       file = "data/Target.rda",       version = 2)
save(Target_group, file = "data/Target_group.rda", version = 2)
save(Student,      file = "data/Student.rda",       version = 2)

message("Wrote data/Target.rda, data/Target_group.rda, data/Student.rda")
