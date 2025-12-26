#' Convert Segment List to Numeric
#'
#' Converts a list of segment indices (potentially stored as character or factor)

#' to numeric values.
#'
#' @param segment_list A list where each element contains segment indices.
#'
#' @return A list of numeric vectors containing the segment indices.
#'
#' @export
#'
#' @examples
#' segments <- list(c("1", "2", "3"), c("10", "11", "12"))
#' convert_segments_to_numeric(segments)
convert_segments_to_numeric <- function(segment_list) {
  numeric_segments <- lapply(segment_list, as.numeric)
  return(numeric_segments)
}


#' Cosine Window Function
#'
#' Computes a smooth window function using cosine tapering for transitioning
#' between blocks. The function is 1 within the block and smoothly transitions
#' to 0 outside using a cosine taper.
#'
#' @param block_indices Integer vector of time indices belonging to the block.
#' @param num_timepoints Total number of time points in the data.
#' @param delta Numeric value controlling the width of the transition region.
#'   Default is 1.
#'
#' @return A numeric vector of length \code{num_timepoints} containing the
#'   window function values.
#'
#' @details
#' For each contiguous segment within \code{block_indices}, the function creates
#' a window that:
#' \itemize{
#'   \item Transitions from 0 to 1 over \code{[a - delta, a]} using a cosine taper
#'   \item Equals 1 over the interior \code{(a, b)}
#'
#'   \item Transitions from 1 to 0 over \code{[b, b + delta]} using a cosine taper
#' }
#' where \code{a} and \code{b} are the start and end of each contiguous segment.
#'
#' @export
#'
#' @examples
#' # Create a window for indices 15-20 out of 50 timepoints
#' w <- window_cos(15:20, num_timepoints = 50, delta = 0.2)
#' plot(w, type = "l", main = "Cosine Window Function")
window_cos <- function(block_indices, num_timepoints, delta = 1) {
  t <- 1:num_timepoints
  # Ensure block_indices is sorted
  block_indices <- sort(block_indices)

  # Identify contiguous segments within block_indices
  # A new segment starts where the gap between consecutive indices is > 1
  diff_indices <- c(TRUE, diff(block_indices) > 1)
  # TRUE for the first element and wherever a gap > 1 occurs

  # End of each segment is right before the next segment start,
  # or the end of block_indices if no more segments
  # We can get segments by splitting at these break points
  segments <- split(block_indices, cumsum(diff_indices))

  # Initialize the overall window w
  w <- rep(0, num_timepoints)

  # Loop through each segment and compute its window contribution
  for (s in seq_along(segments)) {
    seg <- segments[[s]]
    a <- min(seg)
    b <- max(seg)

    w_seg <- rep(0, num_timepoints)

    # Transition from 0 to 1 over [a - delta, a]
    idx1 <- (t >= (a - delta)) & (t <= a)
    w_seg[idx1] <- 0.5 * (1 + cos(pi * (t[idx1] - a + delta) / delta))

    # w = 1 over [a, b]
    idx2 <- (t > a) & (t < b)
    w_seg[idx2] <- 1

    # Transition from 1 to 0 over [b, b + delta]
    idx3 <- (t >= b) & (t <= (b + delta))
    w_seg[idx3] <- 0.5 * (1 + cos(pi * (t[idx3] - b) / delta))

    # Combine this segment's window with the overall window
    # Using pmax ensures w doesn't exceed 1 and handles overlaps gracefully
    w <- pmax(w, w_seg)
  }

  return(w)
}


#' Epanechnikov Window Function
#'
#' Computes a smooth window function using Epanechnikov kernel tapering for
#' transitioning between blocks. The function is 1 within the block and
#' smoothly transitions to 0 outside using an Epanechnikov kernel.
#'
#' @param block_indices Integer vector of time indices belonging to the block.
#' @param num_timepoints Total number of time points in the data.
#' @param delta Numeric value controlling the width of the transition region.
#'   Default is 1.
#'
#' @return A numeric vector of length \code{num_timepoints} containing the
#'   window function values.
#'
#' @details
#' For each contiguous segment within \code{block_indices}, the function creates
#' a window that:
#' \itemize{
#'   \item Transitions from 0 to 1 over \code{[a - delta, a]} using an
#'     Epanechnikov kernel (1 - u^2)
#'   \item Equals 1 over the interior \code{(a, b)}
#'   \item Transitions from 1 to 0 over \code{[b, b + delta]} using an
#'     Epanechnikov kernel
#' }
#' where \code{a} and \code{b} are the start and end of each contiguous segment.
#'
#' @export
#'
#' @examples
#' # Create a window for indices 15-20 out of 50 timepoints
#' w <- window_epan(15:20, num_timepoints = 50, delta = 0.2)
#' plot(w, type = "l", main = "Epanechnikov Window Function")
window_epan <- function(block_indices, num_timepoints, delta = 1) {
  t <- 1:num_timepoints
  # Ensure block_indices is sorted
  block_indices <- sort(block_indices)

  # Identify contiguous segments within block_indices
  segments <- split(block_indices, cumsum(c(TRUE, diff(block_indices) > 1)))

  # Initialize the overall window w
  w <- rep(0, num_timepoints)

  # Loop through each segment and compute its window contribution using the Epanechnikov kernel
  for (seg in segments) {
    a <- min(seg)
    b <- max(seg)

    w_seg <- rep(0, num_timepoints)

    # Left transition: from 0 to 1 over [a - delta, a]
    idx1 <- (t >= (a - delta)) & (t <= a)
    u1 <- (t[idx1] - a) / delta  # u1 ranges from -1 to 0
    w_seg[idx1] <- 1 - u1^2       # Epanechnikov kernel

    # Central region: constant value 1 over (a, b)
    idx2 <- (t > a) & (t < b)
    w_seg[idx2] <- 1

    # Right transition: from 1 to 0 over [b, b + delta]
    idx3 <- (t >= b) & (t <= (b + delta))
    u3 <- (t[idx3] - b) / delta  # u3 ranges from 0 to 1
    w_seg[idx3] <- 1 - u3^2       # Epanechnikov kernel

    # Combine this segment's window with the overall window,
    # using pmax to handle overlaps gracefully
    w <- pmax(w, w_seg)
  }

  return(w)
}


#' Localized Functional Principal Component Analysis
#'
#' Performs Localized Functional Principal Component Analysis (LFPCA) by first
#' detecting block-diagonal structure in functional data, then computing
#' principal components within each localized block. 
#' 
#'
#' @param data A numeric matrix where rows are observations (curves) and
#'   columns are time points.
#' @param Tmax Integer specifying the number of time points.
#' @param range Numeric vector of length 2 specifying the range of the time
#'   domain, e.g., \code{c(0, 1)}.
#' @param Nsel Integer specifying the number of principal components to select.
#' @param smooth.data Logical indicating whether to smooth the data using
#'   \code{\link[refund]{fpca.face}} before block detection. Default is
#'   \code{TRUE}.
#' @param pve Numeric value between 0 and 1 specifying the proportion of
#'   variance explained for the initial smoothing step. Default is 0.99.
#' @param delta Numeric value controlling the width of the window function
#'   transition region. Default is 1.
#' @param spar Smoothing parameter passed to \code{\link[stats]{smooth.spline}}
#'   for smoothing the tapered eigenfunctions. Default is 0.5.
#' @param epsilon Numeric threshold for minimum block size as a proportion of
#'   total time points. Blocks smaller than \code{epsilon * Tmax} are ignored.
#'   Default is 0.05.
#' @param window_type Character string specifying the window function type.
#'   Either \code{"cos"} for cosine tapering or \code{"epan"} for Epanechnikov
#'   kernel. Default is \code{"cos"}.
#' @param anp Character string passed to \code{\link[bdsvd]{bdsvd}} controlling
#'   the regularization function used for the HBIC. Default is \code{"2"}.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{phi.df}{A data frame containing the smoothed eigenfunctions with columns:
#'     \code{nr} (PC number as factor), \code{phi} (eigenfunction values),
#'     \code{t} (time points), \code{Subprocess} (block number as factor),
#'     and \code{eigenvalue}.}
#'   \item{blocks}{A list of numeric vectors containing the time indices for
#'     each detected block.}
#'   \item{eigenvalues}{Numeric vector of the top \code{Nsel} eigenvalues in
#'     decreasing order.}
#'   \item{subprocess_numbers}{Integer vector indicating which block each of
#'     the top eigenvalues came from.}
#' }
#'
#' @details
#' The LFPCA algorithm proceeds as follows:
#' \enumerate{
#'   \item Optionally smooth the input data using FPCA FACE.
#'   \item Detect block-diagonal structure using the bdsvd algorithm.
#'   \item For each block larger than \code{epsilon * Tmax}, compute local
#'     FPCA eigenfunctions.
#'   \item Select the top \code{Nsel} eigenfunctions across all blocks based
#'     on eigenvalue magnitude.
#'   \item Apply window tapering and smoothing to create smooth, localized
#'     eigenfunctions.
#' }
#'
#' @references
#' For the bdsvd block detection method, see the \pkg{bdsvd} package.
#' For FPCA FACE, see the \pkg{refund} package.
#'
#' @importFrom bdsvd bdsvd
#' @importFrom refund fpca.face
#' @importFrom stats smooth.spline predict
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load the pinch force data
#' library(fda)
#' data(pinch)
#' X <- t(pinch)  # 301 time points x 24 curves
#'
#' # Run LFPCA
#' result <- LFPCA(
#'   data = X,
#'   Tmax = ncol(X),
#'   range = c(1, ncol(X)),
#'   Nsel = 6,
#'   smooth.data = TRUE,
#'   window_type = "cos"
#' )
#'
#' # Access the eigenfunctions
#' phi_df <- result[[1]]
#'
#' # Access the eigenvalues
#' eigenvalues <- result[[3]]
#' }
LFPCA <- function(data, Tmax, range, Nsel, smooth.data = TRUE, pve = 0.99,
                  delta = 1, spar = 0.5, epsilon = 0.05, window_type = "cos",
                  anp = "2") {

  t <- seq(range[1], range[2], length.out = Tmax)
  num_timepoints <- Tmax

  # Smooth data if requested
  if (smooth.data == TRUE) {
    data <- refund::fpca.face(data, pve = pve)$Yhat
  }

  # Detect the blocks using bdSVD on the smoothed data
  bl <- bdsvd::bdsvd(data, anp = anp)

  # Convert the blocks to indices
  bl.num <- convert_segments_to_numeric(bl)

  eigenvalues_list <- list()
  eigenvectors_list <- list()
  block_info_list <- list()
  eigen_df_list <- list()

  # Loop over all blocks to compute eigenvalues and eigenvectors
  for (i in seq_along(bl.num)) {

    block_indices <- bl.num[[i]]

    if (length(block_indices) >= round(num_timepoints * epsilon, digits = 0)) {
      data_block <- data[, block_indices]

      # Obtain eigenvalues and eigenfunctions in the block of data
      if (length(block_indices) < 39) {
        fpca <- refund::fpca.face(data_block, pve = 0.99,
                                  knots = (length(block_indices) - 4))
      } else {
        fpca <- refund::fpca.face(data_block, pve = 0.99)
      }

      # Store the results for this block
      eigenvalues_list[[i]] <- fpca$evalues
      eigenvectors_list[[i]] <- fpca$efunctions
      block_info_list[[i]] <- block_indices

      # Create a data frame for this block: note that eigenvalue here is a proxy for lambda.
      n_eigenvalues <- length(fpca$evalues)
      eigen_df_block <- data.frame(
        eigenvalue = fpca$evalues,
        block_number = i,
        eigenvector_index = 1:n_eigenvalues
      )

      eigen_df_list[[i]] <- eigen_df_block
    }
  }

  # Combine all eigenvalues into a single data frame
  eigen_df <- do.call(rbind, eigen_df_list)

  # Sort the eigenvalues in decreasing order and select the top N
  eigen_df_sorted <- eigen_df[order(eigen_df$eigenvalue, decreasing = TRUE), ]

  N <- Nsel
  if (N > nrow(eigen_df_sorted)) {
    message("N exceeds maximum number of PCs, using ", nrow(eigen_df_sorted),
            " components instead.")
    N <- nrow(eigen_df_sorted)
  }

  top_N_eigen_df <- eigen_df_sorted[1:N, ]

  # Instead of explained variance, we now want the eigenvalues (lambdas)
  # and the corresponding subprocess numbers (block numbers) in decreasing order.

  top_eigenvalues <- top_N_eigen_df$eigenvalue
  subprocess_numbers <- top_N_eigen_df$block_number

  # Process eigenfunctions for the top N eigenvalues
  long.phi <- c()
  long.block <- c()
  nr.phi <- c()
  long.eigenvalues <- c()  # New vector for eigenvalues

  for (j in 1:N) {
    # Get eigenvalue and associated indices from the top_N list
    eigenvalue <- top_N_eigen_df$eigenvalue[j]
    block_number <- top_N_eigen_df$block_number[j]
    eigenvector_index <- top_N_eigen_df$eigenvector_index[j]

    # Retrieve the eigenvector from its block
    if (N > 1 & !is.null(dim(eigenvectors_list[[block_number]]))) {
      phi <- eigenvectors_list[[block_number]][, eigenvector_index]
    } else {
      phi <- eigenvectors_list[[block_number]]
    }

    # Initialize a full-length eigenvector and place the block eigenvector at the correct indices
    phi_full <- rep(0, num_timepoints)
    block_indices <- block_info_list[[block_number]]
    phi_full[block_indices] <- phi

    # Define the window function w(t) for the block
    if (window_type == "cos") {
      w <- window_cos(block_indices, num_timepoints = num_timepoints, delta = delta)
    }
    if (window_type == "epan") {
      w <- window_epan(block_indices, num_timepoints = num_timepoints, delta = delta)
    }

    # Multiply by window function
    phi_w <- phi_full * w

    # Smooth the tapered eigenfunction
    spline_fit <- stats::smooth.spline(t, phi_w, spar = spar)
    phi_smoothed <- stats::predict(spline_fit, t)$y

    # Save the smoothed eigenfunction, the block and the PC number (1,...,N) in long format
    long.phi <- c(long.phi, phi_smoothed)
    long.block <- c(long.block, rep(block_number, num_timepoints))
    nr.phi <- c(nr.phi, rep(j, num_timepoints))
    long.eigenvalues <- c(long.eigenvalues, rep(eigenvalue, num_timepoints))
  }

  # Create a data frame for the processed eigenfunctions with the extra eigenvalue column
  phi.df <- data.frame(
    nr = factor(nr.phi),
    phi = long.phi,
    t = rep(t, N),
    Subprocess = factor(long.block),
    eigenvalue = long.eigenvalues
  )

  res <- list(
    phi.df = phi.df,
    blocks = bl.num,
    eigenvalues = top_eigenvalues,
    subprocess_numbers = subprocess_numbers
  )

  return(res)
}
