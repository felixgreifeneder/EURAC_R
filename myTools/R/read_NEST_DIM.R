#this routine reads a NEST .dim file and returns a structure holding the information

library(date)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

read_NEST_DIM <- function(path){


  #define structure
  metadata <- list(DATASET_NAME = as.character(), 
                   ACQUISITION_DATE = as.numeric(), 
                    PASS = as.character(), 
                    POLARIZATION = as.character())

  dim_file <- readLines(con=path)

  #DATASET NAME
  tmp_line <- dim_file[9]
  tmp_line <- trim(tmp_line)
  metadata$DATASET_NAME <- substr(tmp_line, 15, 73)

  #ACQUISITION DATE
  tmp_line <- dim_file[17]
  tmp_line <- trim(tmp_line)
  metadata$ACQUISITION_DATE <- as.numeric(as.date(substr(tmp_line, 34, 44), order="dmy")) - 3653 + 2440588

  #PASS - ASCENDING/DESCENDING
  tmp_line <- dim_file[123]
  tmp_line <- trim(tmp_line)
  metadata$PASS <- substr(tmp_line, 75, 77)

  #POLARIZATION
  tmp_line <- dim_file[125]
  tmp_line <- trim(tmp_line)
  metadata$POLARIZATION <- substr(tmp_line, 76, 77)

  return(metadata)
}
