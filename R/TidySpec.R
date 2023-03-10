# tidy spectra classes, methods and functions

#' Make a spec Object
#'
#' @param x A named list with (at least) the elements freq, spec, and dof
#'
#' @return Object of class spec
#' @export
#' @family TidySpec
as.spec <- function(x){
  class(x) <- unique(append(c("spec", "list"), class(x)) )
  x
}


#' Make a spec_df Object
#'
#' @param x A dataframe or tibble with (at least) columns freq, spec, and dof
#'
#' @return Object of class spec_df
#' @export
#' @family TidySpec
as_spec_df <- function(x) {
  class(x) <- unique(append(c("spec_df", "list"), class(x)))
  x
}



#' Transform a spec Object Into a Dataframe
#'
#' @param x A spec object
#'
#' @return A dataframe or tibble (if package tibble is installed)
#' @export
#' @family TidySpec
#' @examples
#' library(PaleoSpec)
#' ts1 <- ts(rnorm(100))
#' sp1 <- SpecMTM(ts1)
#' sp1_df <- as.data.frame(sp1)
as.data.frame.spec <- function(x){

  df <- data.frame(
    freq = x$freq,
    spec = x$spec
  )

  if (exists("dof", x)){
    df$dof <- x$dof
  }

  if (exists("lim.1", x)){

    df$lim.1 <- x$lim.1
    df$lim.2 <- x$lim.2

  }

  if (require("tibble", character.only = TRUE)){
    df <- tibble::as_tibble(df)
  }

  df <- as_spec_df(df)

  return(df)

}


#' Transform Spec Object(s) Into a Dataframe
#'
#' @param x A spec object or (optionally named) list of spec objects.
#'
#' @return A spec_df object
#' @export
#' @family TidySpec
#' @importFrom tibble as_tibble
#' @importFrom data.table rbindlist
#'
#' @examples
#' library(PaleoSpec)
#' ts1 <- ts(SimPLS(1000, beta = 1))
#' sp1 <- SpecMTM(ts1)
#' sp1 <- AddConfInterval(sp1)
#' ts2 <- ts(rnorm(1000))
#' sp2 <- SpecMTM(ts2)
#' sp_lst <- list(sp1 = sp1, sp2 = sp2)
#' sp_df <- Spec2DF(sp_lst)
#' sp_df
Spec2DF <- function(x){

  if ("spec" %in% class(x)){
    x <- list(x)
  }

  df.lst <- lapply(x, as.data.frame.spec)

  df <- data.table::rbindlist(df.lst, fill = TRUE, idcol = "spec_id")

  if (require("tibble", character.only = TRUE)){
  df <- dplyr::as_tibble(df)
  }

  class(df) <- c("spec_df", class(df))

  return(df)
}


#' Transform a spec_df Object into a spec Object or List of spec Objects
#'
#' @param spec_df
#'
#' @return A spec object or a list of spec objects
#' @export
#' @family TidySpec
#' @examples
#' library(PaleoSpec)
#' ts1 <- ts(SimPLS(1000, beta = 1))
#' sp1 <- SpecMTM(ts1)
#' ts2 <- ts(rnorm(1000))
#' sp2 <- SpecMTM(ts2)
#' sp_lst <- list(sp1 = sp1, sp2 = sp2)
#' sp_df <- Spec2DF(sp_lst)
#' sp_df
#' DF2Spec(sp_df)
DF2Spec <- function(spec_df){

  stopifnot(c("freq", "spec") %in% names(spec_df))

  if ("spec_id" %in% names(spec_df)) {
    lst <- split(spec_df, spec_df$spec_id)
  } else {
    lst <- list(spec_df)
  }

  lst <- lapply(lst, function(x) as.spec(as.list(x)))


  if (length(lst) == 1) {
    lst <- lst[[1]]
  }

  class(lst) <- unique(class(lst))


  return(lst)
}
