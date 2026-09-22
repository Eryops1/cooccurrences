# customfunctions script


# Length unique -----------------------------------------------------------

lunique = function(x){length(unique(x))}

theme_maps <- function() {
  theme_void()+
    theme(
      plot.title = element_text(
        size = 10),
    )
}




# sort species pairs ------------------------------------------------------

sort_pairs_vectorized <- function(x) {
  # a function that sorts species pair names string into alphabetic order
  sapply(x, function(item) {
    spl <- unlist(strsplit(item, split = "|", fixed = TRUE))
    paste0(sort(spl), collapse = "|")
  })
}






