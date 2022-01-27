# Julia Interop -----------------------------------------------------------
# This section introduces some helper functions to interact with
# GenomePermutations.jl

bed2Julia <- \(segs = NULL, file = NULL, name = "mycollection") {
  #' send a bed file or df to Julia as a GenomicFeatures.IntervalCollection
  #' with the given name
  #'
  #' @param segs a dataframe with the segments.
  #' The first columns must contain chr, start, end
  #' @param file must be a tab separated file with
  #' first columns must contain chr, start, end
  #' @param name the name to assign the object to
  #' @note Make sure GenomicFeatures is loaded in Julia
  #'
  #' @returns (in Julia) a GenomicFeatures.IntervalCollection without metadata
  #' @export

  # open file
  if (xor(is.null(segs), is.null(file))) {
    if (is.null(segs)) {
      segs <- as.data.frame(data.table::fread(file, sep = "\t"))
    }
  } else {
    stop("Please specify either a seg object or a file")
  }

  # Sort the file
  segs <- dplyr::arrange(segs, get(colnames(segs)[1]), get(colnames(segs)[2]))

  # assign tmp vectors
  JuliaCall::julia_assign("tmp1", segs[, 1])
  JuliaCall::julia_assign("tmp2", as.integer(segs[, 2]))
  JuliaCall::julia_assign("tmp3", as.integer(segs[, 3]))

  # generate interval collection
  JuliaCall::julia_command(paste0(
    name,
    "=GenomicFeatures.IntervalCollection([GenomicFeatures.Interval(tmp1[i],
    tmp2[i], tmp3[i]) for i in 1:length(tmp1)])"
  ))
  return(invisible(NULL))
}


# Unexported as GenomePermutations.jl is not publicly released yet
# and breaking changes can occur
plotPermResults <- \(results) {
  #' Plots the permutation test in GenomePermutations.jl in ggplot2
  #'
  #' @param results the Results from the permutation test in Julia
  #' @returns the ggplot object
  ggplot2::ggplot() +
    aes(field(results, "rand")) +
    geom_histogram(binwidth = 1, fill = "grey41") +
    ggplot2::geom_vline(
      xintercept = field(results, "obs"),
      linetype = "solid", color = "red", size = 1.5
    ) +
    ggplot2::geom_vline(
      xintercept = mean(field(results, "rand")),
      linetype = "solid", color = "Black", size = 1.5
    ) +
    ggplot2::annotate("text",
      x = field(results, "obs") + 1,
      y = max(table(field(results, "rand"))), label = "Ev"
    ) +
    ggplot2::annotate("text",
      x = mean(field(results, "rand")) + 1,
      y = max(table(field(results, "rand"))), label = "Ov"
    ) +
    ggplot2::ggtitle(
      label = paste0(
        "Feature results of ", field(results, "tested_regions"), " with ",
        field(results, "randomised_regions")
      ),
      subtitle = paste0(
        "p-value: ", field(results, "p_val"), "\t minimum p-value ",
        field(results, "min_p_val"),
        "\nNumber of Permutations: ", field(results, "iterations"),
        "\nAlternative: ", field(results, "alternative")
      )
    )
}
