#' Combines sample and population data
#'
#' This function combines different data sets.
#'
#' @param pop_data a data frame with population data.
#' @param pop_domains population domains.
#' @param smp_data sample data.
#' @param smp_domains sample domains.
#' @param vardir direct variance.
#' @return combined data set.
#' @export


combine_data <- function(pop_data, pop_domains, smp_data, smp_domains, vardir) {

  smp_domains_vec <- smp_data[, smp_domains]
  pop_domains_vec <- pop_data[, pop_domains]

  obs_dom <- pop_domains_vec %in% smp_domains_vec

  smp_data$idD <- pop_domains_vec[obs_dom == TRUE]
  pop_data$idD <- pop_domains_vec
  data <- merge(smp_data, pop_data, by = "idD", all = TRUE)

  if (!is.null(vardir)) {
    data$vardir <- NA
    data$vardir[obs_dom == TRUE] <- vardir
  }

  return(data)
}
