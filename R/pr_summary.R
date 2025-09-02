pr_summary <- function(epreds, group) {
    names <- unique(epreds$group)
    cmbs <- t(utils::combn(names, 2, simplify = TRUE))
    res <- dplyr::tibble(
        Comparison = character(),
        Difference = double(),
        CrI.low = double(),
        CrI.high = double(),
        Pr = double()
    )
    for (i in 1:nrow(cmbs)) {
        comp <- paste(cmbs[i, 1], cmbs[i, 2], sep = " - ")
        draws1 <- epreds %>%
            dplyr::filter(group %in% cmbs[i,1])
        draws2 <- epreds %>%
            dplyr::filter(group %in% cmbs[i,2])
        diff <- draws1$.epred - draws2$.epred
        cri <- quantile(diff, probs = c(0.05, 0.95))
        res <- res %>%
            dplyr::add_row(
                Comparison = comp,
                Difference = mean(diff),
                CrI.low = cri[1],
                CrI.high = cri[2], 
                Pr = bayestestR::p_direction(diff)$pd
            )
    }
    return(res)
}
