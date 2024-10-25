#' Plot One or More Spectra with ggplot2
#'
#' @param x An object of class "spec" or "spec_df"
#' @param gg An existing ggplot object on which to add a new spec layer
#' @param spec_id Name for this plot layer
#' @param conf Plot shaded confidence interval if it exists in the spec object
#' @param alpha.line Alpha level for the spectra line(s)
#' @param alpha.ribbon Alpha level for the confidence region(s)
#' @param colour Name of variable map to colour, unquoted
#' @param linetype Name of variable map to linetype, unquoted
#' @param group Name of variable to map to group, unquoted
#' @param min.colours Minimum number of spectra before starting to colour them
#' separately
#' @param removeFirst,removeLast Remove first or last "n" values on the low or
#'   high frequency side respectively. Operates on a per group basis.
#' @param force.lims,force.CI Force the plotting of confidence regions when the total number
#' of frequencies exceeds 10000. force.CI is deprecated, use force.lims. Defaults to FALSE
#' @param quantiles Plot uncertainty quantiles if they exist in the spec object
#' @param time_unit Optional string giving the time unit for axes labels, e.g. "years"
#'
#' @return a ggplot object
#' @family functions to plot power spectra
#' @export
#'
#' @examples
#' library(PaleoSpec)
#' N <- 1e03
#' beta <- 1
#' alpha <- 0.1
#'
#' ts1 <- SimPLS(N = N, b = beta, a = alpha)
#' ts2 <- SimPLS(N = N, b = beta, a = alpha)
#' sp1 <- SpecMTM(ts1, deltat = 1)
#' sp1 <- AddConfInterval(sp1)
#' sp2 <- SpecMTM(ts2, deltat = 1)
#'
#' # plot single spectrum
#' p <- gg_spec(sp1, spec_id = "df1")
#' p
#'
#' # Add additional second spectra
#' p <- gg_spec(sp2, p, spec_id = "df2", removeFirst = 2)
#' p
#'
#' p <- gg_spec(sp1, p, spec_id = "df3", removeLast = 200)
#' p
#'
#' sp2 <- LogSmooth(sp1)
#' p <- gg_spec(sp1, spec_id = "df1")
#' p <- gg_spec(sp2, p, spec_id = "df2")
#' p <- p + ggplot2::geom_abline(intercept = log10(alpha), slope = -beta, colour = "red")
#' p
#'
#' # Or directly plot named or unnamed list
#'
#' gg_spec(list(sp1, sp2))
#' gg_spec(list(raw = sp1, smoothed = sp2))
#'
#' # without setting any names all spectra will be black
#' p <- gg_spec(sp1)
#' sp2 <- LogSmooth(sp1)
#' p <- gg_spec(sp2, p)
#' p
gg_spec <- function(x, gg = NULL,
                    conf = TRUE,
                    spec_id = NULL,
                    colour = spec_id,
                    group = spec_id,
                    linetype = NULL,
                    alpha.line = 1,
                    alpha.ribbon = c(0.166, 0.333),
                    removeFirst = 0, removeLast = 0,
                    min.colours = 2,
                    force.lims = FALSE,
                    force.CI = NULL,
                    quantiles = FALSE,
                    time_unit = NULL) {
  if (!missing("force.CI")) {

    warning("the force.CI argument is deprecated please use force.lims")
    force.lims <- force.CI
  }

  gg_installed <- requireNamespace("ggplot2", quietly = TRUE)

  if (gg_installed == FALSE) {
    stop("package ggplot2 is required to use gg_spec(). To install ggplot2, run install.packages(\"ggplot2\") from the console")
  }

  if (class(x)[1] != "spec_df" & "list" %in% class(x) == FALSE) {
    x <- list(x)
    names(x) <- spec_id
  }

  if (class(x)[1] != "spec_df") {
    df <- Spec2DF(x)
  } else {
    df <- x
  }


  if (is.numeric(df$spec_id)){
    df$spec_id <- as.character(df$spec_id)
  }


  if (removeFirst > 0) {
    # Convert 'df' to a data.table
    data.table::setDT(df)

    # Group by 'spec_id' and filter rows where 'rank(freq)' is greater than 'removeFirst'
    df <- df[, .SD[rank(freq) > removeFirst], by = spec_id]

    # return to data.frame
    df <- as.data.frame(df)
  }

  if (removeLast > 0) {
    # Convert 'df' to a data.table
    data.table::setDT(df)

    # Group by 'spec_id' and filter rows where 'rank(-freq)' is greater than 'removeLast'
    df <- df[, .SD[rank(-freq) > removeLast], by = spec_id]

    # return to data.frame
    df <- as.data.frame(df)

  }

  if (is.null(time_unit) == FALSE) {
    lab_timescale <- paste0("Timescale [", time_unit, "]")
    lab_freq <- paste0("Frequency [1/", time_unit, "]")
  } else {
    lab_timescale <- "Timescale"
    lab_freq <- "Frequency"
  }


  if (is.null(gg)) {
    p <- ggplot2::ggplot(data = df, ggplot2::aes(group = {{ group }}))
    } else {
    p <- gg
  }

  # rename to PSD so that y axis label can be overwritten later
  df$PSD <- df$spec
  df$Frequency <- df$freq

  df <- as.data.frame(df)

  if (conf == TRUE & exists("lim.1", df)) {
    if (nrow(df) > 1e04 & force.lims == FALSE) {
      warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
    } else {
    p <- p +
      ggplot2::geom_ribbon(
          data = df, ggplot2::aes(
            x = Frequency, ymin = lim.2,
            ymax = lim.1, fill = {{ colour }}
          ),
          alpha = alpha.ribbon[1], colour = NA
        )
    }
    }


  if (quantiles == TRUE & exists("X2.5.", df)) {
    if (nrow(df) > 1e04 & force.lims == FALSE) {
      warning("geom_ribbon is very slow when the number of points > 1e04, skipping the confidence region")
    } else {
      p <- p +
        ggplot2::geom_ribbon(
          data = df, ggplot2::aes(
            x = Frequency, ymin = `X2.5.`,
            ymax = `X97.5.`, fill = {{ colour }}),
          alpha = alpha.ribbon[1], colour = NA
        ) +
        ggplot2::geom_ribbon(
          data = df, ggplot2::aes(
            x = Frequency, ymin = `X15.9.`,
            ymax = `X84.1.`, fill = {{ colour }}
          ),
          alpha = alpha.ribbon[2], colour = NA
        )
    }
  }

  p <- p + ggplot2::geom_line(
    data = df, ggplot2::aes(
      x = Frequency, y = PSD,
      linetype = {{ linetype }},
      colour = {{ colour }}
    ),
    alpha = alpha.line
  ) +
    ggplot2::scale_x_continuous(
      name = lab_freq, trans = "log10",
      sec.axis = ggplot2::sec_axis(~ 1 / ., lab_timescale, labels = scales::comma)
    ) +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::annotation_logticks(sides = "tlb") +
    ggplot2::theme_bw() +
    ggplot2::scale_colour_brewer("",
      type = "qual",
      palette = "Dark2",
      aesthetics = c("colour", "fill")
    ) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::scale_alpha() +
    ggplot2::labs(linetype = "")

  g <- ggplot2::ggplot_build(p)
  colrs <- unlist(unique(sapply(g$data, function(x) unique(x["colour"])$colour)))
  colrs <- colrs[is.na(colrs) == FALSE]
  ncolrs <- length(colrs)


  if (ncolrs <= min.colours) {
    p <- p + ggplot2::scale_colour_manual("",
      values = "black",
      aesthetics = c("colour", "fill")
    )

    if (is.null({{ colour }}) & is.null({{ linetype }})) {
      p <- p + ggplot2::theme(legend.position = "none")
    }
  }

  p
}
