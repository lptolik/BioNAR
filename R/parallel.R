makeBPparam <- function(BPparam=NULL,RNGseed = NULL){
    if (is.null(BPparam)) {
            if (.Platform$OS.type == "windows") {
                # windows doesn't support multicore, using snow instead
                result <- SerialParam(progressbar = TRUE,RNGseed = RNGseed)
            } else {
                result <- MulticoreParam(progress = TRUE,RNGseed = RNGseed)
            }
        return(result)
    }
    else {
        return(BPparam)
    }
}
