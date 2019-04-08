#' @title Format MS2 spectra for export in GNPS-MGF format
#'
#' @description
#'
#' Re-format MS2 spectrum information for export of the data in Mgf format
#' supported by GNPS. In detail, the function replaces the acquisition number
#' of each spectrum with the feature ID (expected to be present in the
#' `"feature_id"` column of `mcols(x)`) converted to an integer by removing
#' the ID's leading `"FT"`.
#'
#' @param x `Spectra`.
#'
#' @return `Spectra` with the acquisition number replaced.
#'
#' @author Johannes Rainer
#' 
#' @md
formatSpectraForGNPS <- function(x) {
    fids <- mcols(x)$feature_id
    if (!length(fids))
        stop("No column named 'feature_id' present in 'mcols(x)'")
    fids <- as.integer(sub("^FT", "", fids))
    mendoapply(x, fids, FUN = function(z, id) {
        z@acquisitionNum <- id
        z
    })
}

#' @title Plot multiple spectra into the same plot
#'
#' @description
#'
#' Plot multiple spectra into the same plot.
#'
#' @param x `Spectra` that should be plotted.
#'
#' @param col color to be used for the individual peaks.
#'
#' @param type `character(1)` defining the plot type. Defaults to `"h"` to plot
#'     vertical lines. For more details see documentation of `plot`.
#'
#' @param main `character(1)` defining the title.
#' 
#' @param ... additional arguments to be passed to `points`.
#'
#' @author Johannes Rainer
#' 
#' @md
plotSpectra <- function(x, col = "#00000040", type = "h", main, ...) {
    plot(3, 3, pch = NA, xlab = "m/z", ylab = "intensity", xlim = range(mz(x)),
         ylim = range(intensity(x)), main = main)
    tmp <- lapply(x, function(z) points(mz(z), intensity(z), type = type, col = col,
                                        ...))
}

#' @title Select spectrum with maximal intensity from a list of spectra
#'
#' @description
#'
#' `maxTic` can be used with the `combineSpectra` method to select the spectrum
#' with the largest overall signal from a list of spectra.
#'
#' @param z `Spectra` object.
#'
#' @return `Spectrum`
#'
#' @author Johannes Rainer
#'
#' @md
maxTic <- function(z) {
    z[[which.max(lapply(intensity(z), sum))]]
}
