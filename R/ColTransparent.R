#' @title Modify a color to get brighter and tranparent for the confidence intervals
#' @param color color value, e.g. "red"
#' @param alpha  (0..1) transparency value
#' @param beta (0..255) to make it brighter, this value gets added on the RGB values
#' @return  modified color
#' @author Thomas Laepple
#' @export
ColTransparent<-function (color, alpha = 0.8,beta=150)
{
    x <- col2rgb(color)[, 1]
    rgb(min(x[1]+beta,255),min(x[2]+beta,255),min(x[3]+beta,255), 255 * alpha, maxColorValue = 255)
}

