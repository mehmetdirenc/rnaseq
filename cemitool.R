if (!("CEMiTool" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("CEMiTool", update = FALSE)
}

library("CEMiTool")