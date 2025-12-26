#' LFPCA: Localized Functional Principal Component Analysis
#'
#' @description
#' The LFPCA package provides tools for performing Localized Functional Principal
#' Component Analysis, a method that combines block-diagonal structure detection
#' with FPCA to produce interpretable, localized eigenfunctions.
#'
#' @details
#' The main function is \code{\link{LFPCA}}, which:
#' \enumerate{
#'   \item Detects block-diagonal structure in functional data using the bdsvd algorithm
#'   \item Computes FPCA within each detected block
#'   \item Combines and ranks eigenfunctions across blocks
#'   \item Applies smooth window tapering for seamless transitions
#' }
#'
#' The package also exports the window functions \code{\link{window_cos}} and
#' \code{\link{window_epan}} for creating smooth tapers, and the utility function
#' \code{\link{convert_segments_to_numeric}}.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{LFPCA}}}{Main function for Localized Functional Principal Component Analysis}
#'   \item{\code{\link{window_cos}}}{Cosine window function}
#'   \item{\code{\link{window_epan}}}{Epanechnikov window function}
#'   \item{\code{\link{convert_segments_to_numeric}}}{Convert segment list to numeric}
#' }
#'
#' @docType package
#' @name LFPCA-package
#' @aliases LFPCA-package
"_PACKAGE"
