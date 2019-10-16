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

#' @title Convert CAMERA output to an edge list for GNPS
#'
#' @description
#'
#' `getEdgelist` takes the output from the `getPeaklist` function from `CAMERA`
#' and converts it to a `data.frame` with edge definitions for *GNPS*. Each
#' row in that `data.frame` contains in columns `"ID1"` and `"ID2"` the
#' identifiers (i.e. `rownames` of the input `data.frame`) of the features.
#' Column `"EdgeType"` is always `"MS1 annotation"` and column `"Score"` `NA`
#' (since no score can be exported from `CAMERA`). Columns `"Annotation"`
#' contains the adduct names and their difference in m/z if **both** edges
#' (features) were predicted to be an adduct of the **same** compound.
#' Column `"Isotope"` provides the same information as for adducts, only for
#' isotope definitions. Column `"pcgroup"` provides the information which
#' features were grouped by `CAMERA` into the same group.
#' 
#' @param peaklist `data.frame` as returned by the [getPeaklist()] function
#'   from `CAMERA` package or an `xsAnnotate` object.
#'
#' @return `data.frame` with edge definitions (see description for more
#'     details).
#'
#' @author Mar Garcia-Aloy, Johannes Rainer
#' 
#' @examples
#' 
#' res <- getEdgelist(getPeaklist(xsaFA))
getEdgelist <- function(peaklist) {
    if (is(peaklist, "xsAnnotate")) {
        peaklist <- getPeaklist(peaklist)
        if (!nrow(peaklist))
            stop("Got an empty peak list.")
    }
    pl <- split(peaklist, factor(peaklist$pcgroup,
                                 levels = unique(peaklist$pcgroup)))
    res <- do.call(rbind, lapply(pl, .process_pcgroup))
    rownames(res) <- NULL
    res
}

#' @title Extract feature annotations from CAMERA results
#'
#' @description
#'
#' Similar to the `getEdgelist` function, this function extracts information
#' from a `CAMERA` result for use in GNPS.
#'
#' @param x `xsAnnotate` object after calling `findAdducts`.
#'
#' @return
#'
#' `data.frame` with columns:
#' - `"annotation network number"`: ion identify network (IIN) number. All
#'   features predicted by `CAMERA` to be an adduct of a (co-eluting) compound
#'   with the same mass are part of this IIN. If a feature was predicted to be
#'   an adduct of two different compounds (with different masses) the ID of the
#'   larger network is reported. All features for which no adduct annotation is
#'   available will have an `NA` in this column.
#' - `"best ion"`: the adduct definition of the feature.
#' - `"correlation group ID"`: this corresponds to the `"pcgroup"` column in
#'   `getPeaklist(x)`.
#' - `"auto MS2 verify"`: always `NA`.
#' - `"identified by n="`: the size of the IIN.
#' - `"partners"`: all other features (rows in the feature table) of this IIN.
#' - `"neutral M mass"`: the mass of the compound.
#'
#' @author Johannes Rainer
#'
#' @noRd
getFeatureAnnotations <- function(x) {
    if (!length(x@annoID))
        stop("No adduct information present. Please call 'findAdducts' on ",
             "the object.")
    corr_group <- rep(seq_along(x@pspectra), lengths(x@pspectra))
    corr_group <- corr_group[order(unlist(x@pspectra, use.names = FALSE))]

    ## Get the all ids (feature rows) for which an adduct was defined
    ids <- unique(x@annoID[, "id"])

    ## Note: @pspectra contains the "correlation groups", @annoID the
    ## adduct groups, but it can happen that two ids are in the same adduct
    ## group without being in the same correlation group!

    ## loop through the adduct definitions and build the output data.frame
    adduct_def <- lapply(ids, function(id) {
        ## IDs of the same correlation group
        ids_pcgroup <- x@pspectra[[corr_group[id]]]
        ## ID of the adduct annotation groups this id is part of 
        anno_grp <- x@annoID[x@annoID[, "id"] == id, "grpID"]
        ## Subset the adduct annotation to rows matching the annotation group
        ## of the present ID and to ids present in the same correlation group.
        adduct_ann <- x@annoID[x@annoID[, "id"] %in% ids_pcgroup &
                               x@annoID[, "grpID"] %in% anno_grp, , drop = FALSE]
        ## if we have more than one annotation group, select the bigger one
        if (length(anno_grp) > 1) {
            cnts <- table(adduct_ann[, "grpID"])
            adduct_ann <- adduct_ann[
                adduct_ann[, "grpID"] ==
                names(cnts)[order(cnts, decreasing = TRUE)][1], ,
                drop = FALSE]
        }
        grp_id <- adduct_ann[1, "grpID"]
        ## different adduct rules can match the same m/z - we're just taking
        ## the first one (with lower ruleID)
        rule_id <- adduct_ann[adduct_ann[, "id"] == id, "ruleID"][1]
        ids_grp <- adduct_ann[, "id"]
        df <- data.frame(
            `row ID` = id,
            `annotation network number` = unname(grp_id),
            `best ion` = as.character(x@ruleset[rule_id, "name"]),
            `identified by n=` = nrow(adduct_ann),
            `partners` = paste0(ids_grp[ids_grp != id], collapse = ";"),
            `neutral M mass` = unname(
                x@annoGrp[x@annoGrp[, "id"] == grp_id, "mass"]),
            stringsAsFactors = FALSE, check.names = FALSE)
    })
    adduct_def <- do.call(rbind, adduct_def)
    res <- adduct_def[rep("other", length(corr_group)), ]
    res[, "row ID"] <- seq_along(corr_group)
    rownames(res) <- as.character(res[, "row ID"])
    res[ids, ] <- adduct_def
    res$`correlation group ID` <- corr_group
    res$`auto MS2 verify` <- NA
    res
}

#' Helper function to extract the adduct annotation from a pair of adducts
#' from the same *pcgroup*.
#'
#' @author Mar Garcia-Aloy
#'
#' @noRd
.define_annot <- function(y) {
    if (any(y$adduct == "")) return(NA)
    mass_1 <- .extract_mass_adduct(y$adduct[1])
    mass_2 <- .extract_mass_adduct(y$adduct[2])
    mass <- intersect(mass_1, mass_2)
    if (length(mass)) {
        def1 <- unlist(strsplit(y$adduct[1], " "))
        def2 <- unlist(strsplit(y$adduct[2], " "))
        paste0(def1[grep(mass[1], def1) - 1], " ",
               def2[grep(mass[1], def2) - 1], " dm/z=",
               round(abs(y$mz[1] - y$mz[2]), 4))
    } else NA
}

#' Helper function to extract the isotope annotation from a pair of adducts
#' from the same *pcgroup*.
#'
#' @author Mar Garcia-Aloy
#'
#' @noRd
.define_isotop <- function(w) {
  if (any(w$isotopes == "")) return(NA)
  if (unlist(strsplit(w$isotopes[1], '\\]\\[') )[1] == 
     unlist(strsplit(w$isotopes[2], '\\]\\[') )[1]) {
    a = paste0("[", do.call(rbind, strsplit(w$isotopes, "\\]\\["))[, 2], 
               collapse = " ")
    b = paste0(unlist(strsplit(w$isotopes[2], '\\]') )[1], "]")
    paste0(b, a, " dm/z=", round(abs(w$mz[1] - w$mz[2]), 4))
  } else NA
}

#' Simple helper to extract the mass(es) from strings such as
#' [M+NH4]+ 70.9681 [M+H]+ 87.9886
#'
#' @author Johannes Rainer
#' 
#' @noRd
#' 
#' @examples
#'
#' .extract_mass_adduct("[M+NH4]+ 70.9681 [M+H]+ 87.9886")
#' .extract_mass_adduct("some 4")
.extract_mass_adduct <- function(x) {
    if (!length(x) || x == "") return(NA)
    spl <- unlist(strsplit(x, " ", fixed = TRUE))
    spl[seq(2, by = 2, length.out = length(spl)/2)]
}

#' Helper function to process features from the same *pcgroup*
#'
#' @author Mar Garcia-Aloy
#'
#' @noRd
.process_pcgroup <- function(x) {
    if (nrow(x) > 1) {
        res <- combn(seq_len(nrow(x)), 2, FUN = function(z) {
            data.frame(ID1 = rownames(x)[z[1]],
                       ID2 =  rownames(x)[z[2]],
                       EdgeType = "MS1 annotation",
                       Score = NA,
                       Annotation = .define_annot(x[z,]),
                       Isotope = .define_isotop(x[z,]),
                       pcgroup = x$pcgroup[1],
                   stringsAsFactors = FALSE)
        }, simplify = FALSE)
        do.call(rbind, res)
    } else NULL
}
