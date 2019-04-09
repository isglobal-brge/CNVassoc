getProbsRegions <- function(probs, regions, 
        intensities, nclass = 3) {
        blocks <- regions@featureData@data
        annotation <- intensities[, seq_len(4)]
        probsBlock <- lapply(seq_len(nrow(blocks)), 
            getProbsRegions.i, blocks = blocks, 
            probs = probs, annotation = annotation, 
            nclass = nclass)
        probsBlock
}
