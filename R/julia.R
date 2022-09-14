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


.parseC <- \(results){
  #' private function to parse common stuff out of any GenomnePermutaiton.Jl permtest
  #'
  #' @param results the julia object from GenomePermutation's permtest
  #' @seealso \code{\link{.parseZ}}, \code{\link{.parseS}}

  return(list(
    "iter" = JuliaCall::field(results, "iterations"),
    "tested" = JuliaCall::field(results, "tested_regions"),
    "randomised" = JuliaCall::field(results, "randomised_regions"),
    "obs" = JuliaCall::field(results, "obs"),
    "ran" = JuliaCall::field(results, "ran"),
    "eval_function" = JuliaCall::field(results, "eval_function"),
    "random_strategy" = JuliaCall::field(results, "random_strategy")
  ))
}

.parseS <- \(results) {
  #' internal functions to parse the results of a GenomnePermutaiton.Jl permtest
  #' This version deals with the simple_P test
  #'
  #' @param results the julia object from GenomePermutation's permtest
  #' @return a list with the fields:
  #' iterm, tested, randomised, obs, ran, eval_function, random_strategy, P, test, notes
  #' @note the notes field contains a line of information for plotting specific to each test
  #' @seealso  \code{\link{.parseZ}}
  return(purrr::flatten(list(
    .parseC(results),
    "p" = JuliaCall::field(results, "p_val"),
    "test" = "Simple empirical p-test",
    "notes" = paste(
      "min p-value:", JuliaCall::field(JuliaCall::field(results, "test"), "min_p_val"),
      "alternative:", JuliaCall::field(JuliaCall::field(results, "test"), "alternative")
    )
  )))
}


.parseZ <- \(results) {
  #' internal functions to parse the results of a GenomnePermutaiton.Jl permtest
  #' This version deals with the HypothesisTets.Z-test
  #'
  #' @param results the julia object from GenomePermutation's permtest
  #' @return a list with the fields:
  #' iterm, tested, randomised, obs, ran, eval_function, random_strategy, P, test, notes
  #' @seealso  \code{\link{.parseS}}

  JuliaCall::julia_install_package_if_needed("HypothesisTests")
  JuliaCall::julia_library("HypothesisTests")
  JuliaCall::julia_assign("tmp", results)
  out <- (purrr::flatten(list(
    .parseC(results),
    "p" = JuliaCall::julia_eval("pvalue(tmp.test)"),
    "test" = "Two sample z-test",
    "notes" = paste(
      "mean difference (ran-obs):",
      JuliaCall::field(JuliaCall::field(results, "test"), "xbar"),
      "CI:", julia_eval("confint(tmp.test)")[[1]], "-", julia_eval("confint(tmp.test)")[[2]]
    )
  )))
  JuliaCall::julia_command("tmp = nothing") # can't delete tmp but can make tiny
  return(out)
}



# Now exported as GenomePermutations.jl is not publicly released yet
# parsing is abstracted away in .parseC, .parseS and .parseZ
plotPermResults <- function(results, test) {
  #' Plots the permutation test in GenomePermutations.jl in ggplot2
  #'
  #' @param results the julia object from GenomePermutation's permtest
  #' @param test the type of test, can (partially) match either "simple" or z-test"
  #' @returns the ggplot object
  #' @export

  # get fields
  tests <- c("simple p test", "z-test")
  t <- match.arg(tolower(test), tests)
  if (t == tests[[1]]) {
    res <- .parseS(results)
  } else if (t == tests[[2]]) {
    res <- .parseZ(results)
  }

  ggplot2::ggplot() +
    aes(res[["ran"]]) +
    geom_histogram(fill = "grey41") +
    ggplot2::ggtitle(
      label = paste0(
        "Feature results of ", res[["tested"]], " with ",
        res[["randomised"]]
      ),
      subtitle = paste0(
        "p-value: ", res[["p_val"]], "iterations:", res[["iter"]], res[["notes"]]
      )
    ) +
    ggplot2::geom_vline(
      xintercept = mean(res[["obs"]]),
      linetype = "solid", color = "red", size = 1.5
    ) +
    ggplot2::geom_vline(
      xintercept = mean(res[["ran"]]),
      linetype = "solid", color = "Black", size = 1.5
    ) +
    ggplot2::annotate("text",
      x = mean(res[["obs"]]) + 10,
      y = res[["iter"]] / 20, label = "Ov"
    ) +
    ggplot2::annotate("text",
      x = mean(res[["ran"]]) + 10,
      y = res[["iter"]] / 20, label = "Ev"
    ) +
    ggplot2::ggtitle(
      label = paste0(
        "permutation result for ", res[["tested"]], " with ",
        res[["randomised"]]
      ),
      subtitle = paste0(
        "p-value: ", res[["p"]], " iterations: ", res[["iter"]], "\n", res[["notes"]]
      )
    )
}
