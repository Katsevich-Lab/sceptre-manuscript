# Processing functions

#' Extract vector from compressed matrix
#'
#' A function to extract a vector(s) from a matrix stored in compressed column or row format
#'
#' @param vector_nos indexes of the vectors to extract
#' @param x sparse-matrix data
#' @param i sparse-matrix indexes
#' @param p sparse-matrix pointers
#' @param vector_length length of the vectors to extract
#' @param zero_based_idx (boolean; default TRUE) do i and p use zero-based indexing?
#' @param vector_subset (optional) integers indicating a subset of the vectors to return
#'
#' @return the extracted vectors
#' @export
extract_vector_from_compressed_matrix <- function(vector_nos, x, i, p, vector_length, zero_based_idx = TRUE, vector_subset = NULL) {
  if (zero_based_idx) {
    i <- i + 1
    p <- p + 1
  }
  f <- function(no) {
    curr_ptr_range <- p[no:(no + 1)] + c(0,-1)
    out <- rep(0, vector_length)
    if (curr_ptr_range[1] <= curr_ptr_range[2]) {
      x_out <- x[curr_ptr_range[1]:curr_ptr_range[2]]
      i_out <- i[curr_ptr_range[1]:curr_ptr_range[2]]
      out[i_out] <- x_out
    }
    if (!is.null(vector_subset)) out <- out[vector_subset]
    return(out)
  }
  ret_val <- if (length(vector_nos) == 1) f(vector_nos) else sapply(X = vector_nos, FUN = f)
  return(ret_val)
}


#' Extract column from compressed sparse column matrix
#'
#' @param col_nos idxs of the columns to extract
#' @param csc_mat a list representing a matrix stored in csc form containing the entries p, i, x, and dim.
#' @param zero_based_idx (boolean) do i and p use zero-based indexing
#' @param row_subset (optional) subset the extracted columns by row index
#'
#' @return a matrix containing the extracted columns
#' @export
extract_column_from_csc_matrix <- function(col_nos, csc_mat, zero_based_idx = TRUE, row_subset = NULL) {
  p <- csc_mat$p
  i <- csc_mat$i
  x <- csc_mat$x
  n_rows <- csc_mat$dim[1]
  out <- extract_vector_from_compressed_matrix(vector_nos = col_nos, x = x, i = i, p = p, vector_length = n_rows, zero_based_idx = zero_based_idx, vector_subset = row_subset)
  return(out)
}

#' Extract row from compressed sparse row matrix
#'
#' This function is similar to extract_column_from_csc_matrix; see the documentation of that function.
#'
#' @export
extract_row_from_csr_matrix <- function(row_nos, csr_mat, zero_based_idx = TRUE, col_subset = NULL) {
  x <- csr_mat$x
  i <- csr_mat$i
  p <- csr_mat$p
  n_cols <- csr_mat$dim[2]
  out <- extract_vector_from_compressed_matrix(vector_nos = row_nos, x = x, i = i, p = p, vector_length = n_cols, zero_based_idx = zero_based_idx, vector_subset = col_subset)
  if (is.matrix(out)) out <- t(out)
  return(out)
}
