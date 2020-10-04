#' Organize estimates and likelihood from probability model fits
#'
#' @param x list of probability distribution fits. Probably a list of outputs from fitdistrplus::fit_dist()
#'
#' @return
#' @export
#'
#' @examples
#'
get_fits_params_df <- function(x){

  nm <- deparse(substitute(x))

  x %>%
    map(., `[`, c("estimate", "sd", "loglik", "n")) %>%
    map(as_tibble) %>%
    bind_rows()
}
