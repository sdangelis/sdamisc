fractionAltered <- \(segs) {
  #' Calculates the fraction of the genome altrered in a given segment table
  #'
  #' @param segs data frame with CNA length,
  #' chromosome and direction ("gain","loss","normal") columns
  #' @return A named list with the fraction gained and lost
  #' @note This function estimates the length of the genome from the length
  #' of all segments in the table
  #' @export

  tot <- sum(segs$length)
  gain <- sum(dplyr::filter(segs, direction == "gain")$length)
  loss <- sum(dplyr::filter(segs, direction == "loss")$length)

  return(list("gain" = gain / tot, "loss" = loss / tot))
}

callWcna <- \(segs, chr, threshold = 0.95) {
  #' Call Whole chromosome aneuploidies from segments based on a threshold
  #'
  #' @note This function estimates the length of the genome
  #' from the length of all segments in the table
  #' @note if the given chromosome is not found returns NA.
  #' @note thresholds below < 0.5 always return gains as long as it passes the threshold,
  #' regardless of whether there are more losses.
  #' @param chr the chromosome of interest. must match the format
  #' of the segment table
  #' @param segs data frame with CNA length,
  #' chromosome and direction ("gain","loss","normal") columns
  #' @param thershold float, the fraction of the chromosome needed to be
  #' in a direction to call. Default = 0.95
  #' @return "gain", "loss", "normal" or NA, a string with the status
  #' of the chromosome
  #' @seealso \code{\link{callArmCna}} \code{\link{fractionAltered}}
  #' @export

  filtered_segs <- dplyr::filter(segs, chromosome == chr)
  total_length <- sum(filtered_segs$length)
  if (total_length < 1) {
    warning(paste("chromosome: ", chr, "not found"))
    return(NA_character_)
  }
  if (threshold < 0.5) {
    warning("thresholds below 0.5 can bias the results towards gains")
  }
  if (sum(dplyr::filter(filtered_segs, direction == "gain")$length)
  / total_length >= threshold) {
    return("gain")
  } else if (sum(dplyr::filter(filtered_segs, direction == "loss")$length)
  / total_length >= threshold) {
    return("loss")
  } else {
    return("normal")
  }
}



callArmCna <- \(segs, chr, threshold = 0.95, arms = bands) {
  #' Wrapper around CallWcna() that calls chromosome arm aneuploidies
  #'
  #' @note This assumes that there is a gap around the centromere and
  #' that no segment can ever overlaps it.
  #' @note for full documentation see callWcna.
  #' @note This function estimates the length of the genome
  #' from the length of all segments in the table
  #' @param chromosome the chromosome of interest.
  #' must match the format of the segment table
  #' @param segments data frame with CNA length,
  #' chromosome and direction ("gain","loss","normal") columns
  #' @param thershold float, the fraction of the chromosome
  #' needed to be in a direction to call. Default = 0.95
  #' @return a named vector with the arm and "gain", "loss", "normal" or NA.
  #' @seealso \code{\link{callWcna}}
  #' @export
  arms <- dplyr::filter(arms, chrom == chr)
  lapply(arms$arm, \(x) {
    arm_seg <- dplyr::filter(
      segs,
      start >= as.numeric(arms[arms$arm == x, "chromStart"]) &
        end <= as.numeric(arms[arms$arm == x, "chromEnd"])
    )
    if (length(arm_seg) == 0) {
      return(NA)
    }
    callWcna(segs = arm_seg, chr = chr, threshold = threshold)
  })
}
