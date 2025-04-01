# Data wrangling

fil <- here::here("inst","extdata","rds.xls")


descriptor <-
  here::here("inst","extdata","description.xlsx") %>%
  file.path() %>%
  readxl::read_excel(skip = 0)

# make a list to be used as labels for 'labelled::'
labs <-
  # take a 1xN matrix
  matrix( descriptor$description,
          ncol = length(descriptor$description)
  ) %>%
  as.list() %>%
  # add column names as 'names'
  `names<-`(descriptor$name_new)

sheets <- readxl::excel_sheets(fil)

for (dataset in sheets) {

  datachunk <- fil |>
    readxl::read_xlsx(sheet = dataset) |> # or read_xls() as appropriate
    # handle duplicate colnames
    janitor::clean_names() |>
    mutate( across( .cols = which( descriptor$trf == "factor"),
                    .fns = as.factor
    ),
    across( .cols = which( descriptor$trf == "numeric"),
            .fns = as.numeric # removing potential '?', 'NA', '.' etc.
    ),
    across( .cols = which( descriptor$trf == "date"),
            .fns = lubridate::as_datetime
    )) %>%
    `colnames<-`( descriptor$name_new) %>%
    # select only columnswhose name isnt NA due to duplicates or whatnot
    select( descriptor$name_new[!is.na(descriptor$name_new)])%>%
    labelled::`var_label<-`(   labs  ) %>%
    mutate(dataset = dataset)

  # binding the actual dataset to a 'master list'
  if ( dataset == sheets[1]) {
    data <- datachunk
  } else {
    data <- bind_rows(data, datachunk)
  }

}

# deselect the "not_relevant" columns
try( {
  data <- data %>%
    select( -not_relevant)
}, silent = TRUE)



