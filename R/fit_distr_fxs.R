#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
build.fits.df <- function(x){

  nm <- deparse(substitute(x))

  x %>%
    map(., `[`, c("estimate", "sd", "loglik", "n")) %>%
    map(as_tibble) %>%
    bind_rows()
}
