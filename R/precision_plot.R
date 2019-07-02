precision_plot <- function(direct, model, indicators, MSE = TRUE, CV = FALSE,
                           label, color, shape, line_type, gg_theme) {

  plotList <- vector(mode = "list", length = 2)
  names(plotList) <- c("ordered", "boxplot")


  if (inherits(model, "ebp")) {
    direct <- as.data.frame(estimators(direct, indicators = indicators,
                                       MSE = MSE, CV = CV))
    model <- as.data.frame(estimators(model, indicators = indicators,
                                      MSE = MSE, CV = CV))
  } else if (inherits(model, "fh")) {
    estims <- as.data.frame(estimators(model, indicators = indicators,
                                      MSE = MSE, CV = CV))

  }






}
