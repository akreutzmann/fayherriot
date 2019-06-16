#' Fitted emdiObject
#'
#' An object of class emdi that represents point predictions of regional
#' disaggregated indicators. Optionally it also contains corresponding MSE
#' estimates. Depending on the estimation the object is also of class direct
#' or model. For each provided model-based approach an additional class is assigned:
#' the standard Fay-Herriot approach ("fh"), and the empirical best prediction
#' ("ebp"). Objects of these classes have methods for the generic functions
#' \code{\link{estimators}}, \code{\link{print}}, \code{\link{plot}} (only for
#' class model), and \code{\link{summary}}.
#'
#' @return
#' The following components are always included in an emdi object:
#' \item{\code{call}}{a list containing an image of the function call that
#'                    produced the object.}
#' \item{\code{fixed}}{a formula of fixed effects used in the nested error linear
#'              regression (see also \code{fixed} in \code{\link{fh}} and
#'              \code{\link{ebp}}). Not filled for class direct.}
#'  \item{\code{framework}}{a list with following components:}
#'
#'  \tabular{rll}{
#'      \tab \code{N_dom_smp} \tab number of domains in the sample. \cr
#'      \tab \code{N_dom_unobs} \tab number of out-of-sample domains. Not
#'              filled for class direct. \cr
#'      \tab \code{N_pop} \tab total number of units in population. Not
#'              filled for class direct. \cr
#'      \tab \code{N_smp} \tab total number of units in sample. \cr
#'      \tab \code{pop_domains_vec} \tab an arranged vector
#'                     of the domain indicator variable. Not
#'              filled for class direct. \cr
#'      \tab \code{smp_data} \tab an arranged data set of sample data. Not
#'              filled for class direct. \cr
#'      \tab \code{smp_domains} \tab a character naming the
#'                     domain indicator variable.\cr
#'      \tab \code{smp_domains_vec} \tab an arranged vector of
#'                     the domain indicator variable. \cr
#' }
#' \item{\code{ind}}{data frame containing estimates for indicators per domain.}
#' \item{\code{method}}{character returning the method for the estimation of
#'    the optimal lambda (for class ebp), here "reml", or a list returning method for
#'    the estimation of the variance of the random effect and the applied MSE
#'    estimation (for class fh). Not filled for class direct.}
#' \item{\code{model}}{for class ebp, an object returned by the lme function of type "lme"
#'    and representing a fitted linear mixed-effects model (for
#'    further explanations see \code{\link[nlme]{lme}} and
#'    \code{\link[nlme]{lmeObject}}). For class fh, a selection of model components.
#'    Not filled for class direct.}
#' \item{\code{MSE}}{data frame containing MSE estimates corresponding to the
#' point predictions in \code{ind} per indicator per domain if MSE is selected
#' in function call. If \code{FALSE}, \code{MSE} is \code{NULL}.}
#' \item{\code{transformation}}{character returning the selected transformation
#' type (see also \code{transformation} in \code{\link{fh}}, and \code{\link{ebp}}).
#' Not filled for class direct.}
#' \item{\code{transform_param}}{a list with two elements, \code{optimal_lambda}
#'    and \code{shift_par}, where the first contains the optimal parameter for a
#'    Box-Cox transformation or NULL for no and log transformation and the
#'    second the potential shift parameter in the log or Box-Cox transformation
#'    and NULL for no transformation. Not filled for class fh and direct.}
#' \item{\code{successful_bootstraps}}{a matrix with domains as rows and indicators
#'      as columns. The cells contain the number of successful bootstraps for each
#'      combination. Not filled for class model.}
#' @references
#' Alfons, A. and Templ, M. (2013). Estimation of Social Exclusion Indicators
#' from Complex Surveys: The \R Package \pkg{laeken}. Journal of
#' Statistical Software, 54(15), 1-25.  \cr \cr
#' Fay R.E., Herriot R.A. (1979) Estimates of income for small places: An
#' application of James–Stein procedures to census data. Journal of the American
#' Statistical Association, Vol. 74, 269–277.  \cr \cr
#' Molina, I. and Rao, J.N.K. (2010). Small area estimation of poverty
#' indicators. The Canadian Journal of Statistics, Vol. 38, No.3, 369-385.
#' @seealso \code{\link{direct}}, \code{\link{ebp}}, \code{ \link[nlme]{lme}},
#' \code{ \link[nlme]{lmeObject}}
#'
#' @name emdiObject
NULL
