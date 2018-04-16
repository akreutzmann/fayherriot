#' CV plots
#'
#' This function returns two different plots for the CV.
#'
#' @param fit_FH fitted FH model.
#' @param label_direct label for the direct estimates.
#' @param label_FH label for FH estimates.
#' @param line benchmark CV level.
#' @param colline color of benchmark line.
#' @param color color of the two estimates.
#' @return two CV plots.
#' @export

CV_plots <- function(fit_FH, label_direct = "Direct",
                     label_FH = "FH-MI", line = 20, colline = "red",
                     color = c("dodgerblue4", "#99CC00")) {

  # Define a data frame with the different variables
  x <- data.frame(Direct = fit_FH$ind$Direct,
                  Var = fit_FH$MSE$Var,
                  FH = fit_FH$ind$EBLUP,
                  MSE = fit_FH$MSE$PR_MSE)

  x <- x[order(abs(sqrt(x$Var) / x$Direct)), ]
  x <- x[!is.na(x$Direct),]
  class(x)

  ggDat <- data.frame(
    predictions = c(x$Direct, x$FH),
    mse = c(x$Var, x$MSE),
    method = c(rep(label_direct, NROW(x)), rep(label_FH, each = NROW(x))),
    area = 1:NROW(x)
  )

  ggDat$cv <- 100 * abs(sqrt(ggDat$mse) / ggDat$predictions)
  class(ggDat)
  is.data.frame(ggDat)

  # Save plots
  if(!is.null(line)) {
    cv_ordered <- ggplot(ggDat, aes(x = area, y = cv, colour = method)) +
      geom_point() + labs(x = "Domain (sorted by increasing CV of Direct)",
                          y = "CV",
                          colour = "Method") +
      scale_color_manual(values = color) +
      theme_minimal(base_size = 25) + geom_hline(yintercept = line, size = 1,
                                                 linetype = 2, col = colline)

    cv_boxplot <- ggplot(ggDat, aes(x = method, y = cv, fill = method)) +
      geom_boxplot() +
      coord_flip() +
      labs(x = NULL, y = "CV") +
      scale_fill_manual(name = "Method",
                        values = color) +
      theme_minimal(base_size = 25) + geom_hline(yintercept = line, size = 1,
                                                 linetype = 2, col = colline)
  } else if (is.null(line)) {
    cv_ordered <- ggplot(ggDat, aes(x = area, y = cv, colour = method)) +
      geom_point(size = 2.5) +
      labs(x = "Domain (sorted by increasing CV of Direct)",
           y = "CV",
           colour = "Method") +
      scale_color_manual(values = color) +
      theme_minimal(base_size = 25)

    cv_boxplot <- ggplot(ggDat, aes(x = method, y = cv, fill = method)) +
      geom_boxplot() +
      coord_flip() +
      labs(x = NULL, y = "CV") +
      scale_fill_manual(name = "Method",
                        values = color) +
      theme_minimal(base_size = 25)

  }



  return(list(cv_ordered = cv_ordered, cv_boxplot = cv_boxplot))
  invisible()
}
