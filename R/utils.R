#' Addition of two 1-D base::table objects
#'
#' @param t1 a base::table object
#' @param t2 a second base::table object
#' @param keep_zeros whether we should retain in the result values with total count 0 (possible if
#' negating one of the added tables).
#' @returns a base::table where the count of every label is the sum of the counts in t1 and t2
#'
#' @export
#'
#' @examples
#' table_add_1d(table(1:10), table(5:15))
table_add_1d <- function(t1, t2, keep_zeros = FALSE) {
  stopifnot(is.table(t1) && is.table(t2))
  stopifnot(length(attr(t1, "dim")) == 1 && length(attr(t2, "dim")) == 1)

  all_labels <- c(names(t1), names(t2)) |> unique()

  res <- sapply(all_labels, function(lab) {
    sum(t1[lab], t2[lab], na.rm = TRUE)
  })
  attr(res, "dim") <- length(res)
  attr(res, "dimnames") <- list(all_labels)
  if (!keep_zeros) {
    keep <- res != 0
    res <- res[keep]
    attr(res, "dimnames") <- list(all_labels[keep])
  }
  class(res) <- "table"

  return(res)
}
