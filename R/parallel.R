makeBPparam <- function(BPparam=NULL,RNGseed = NULL){
    if (is.null(BPparam)) {
            if (.Platform$OS.type == "windows") {
                # windows doesn't support multicore, using snow instead
                result <- SerialParam(progressbar = TRUE,RNGseed = RNGseed)
            } else {
                chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")

                if (nzchar(chk) && chk == "TRUE") {
                    # use 2 cores in CRAN/Travis/AppVeyor
                    num_workers <- 2L
                } else {
                    # use all cores in devtools::test()
                    num_workers <- multicoreWorkers()
                }
                result <- MulticoreParam(workers=num_workers,progress = TRUE,
                                         RNGseed = RNGseed)
            }
        return(result)
    }
    else {
        return(BPparam)
    }
}
