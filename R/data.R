#' Ammonia Oxidation Data (Brownlee, 1965)
#'
#' Data from an industrial ammonia oxidation process recorded over 21 days
#' (Brownlee, 1965, p. 45).  The response is the proportion of ammonia not
#' converted to nitric acid; three process variables are included.
#'
#' @format A data frame with 21 rows and 4 columns:
#' \describe{
#'   \item{\code{perda}}{Proportion of ammonia loss (response), in (0, 1).}
#'   \item{\code{corr_ar}}{Air flow rate (arbitrary units).}
#'   \item{\code{temp_agua}}{Cooling water inlet temperature (degrees Celsius).}
#'   \item{\code{conc_acido}}{Nitric acid concentration (percent).}
#' }
#'
#' @details
#' The model selected in Espinheira et al. (2026) is:
#' \deqn{\mathrm{logit}(\mu_t) = \beta_1 + \beta_2 x_{t2} + \beta_3 x_{t3}
#'       + \beta_4 (x_{t2} \times x_{t3})}
#' \deqn{\log(\sigma^2_t) = \gamma_1 + \gamma_2 x_{t3} +
#'       \gamma_3 (x_{t2} \times x_{t3})}
#' where \eqn{x_{t2}} = \code{corr_ar}, \eqn{x_{t3}} = \code{temp_agua}.
#' The variable \code{conc_acido} is included in the dataset for completeness
#' but was not used in the selected model.
#'
#' @source Brownlee, K.A. (1965). \emph{Statistical Theory and Methodology in
#'   Science and Engineering}, 2nd ed. Wiley, New York, p. 45.
#'
#' @references
#' Espinheira P.L., Silva F.C., Barros M., Ospina R. (2026).
#' A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
#' Simplex Regression Models.
#'
#' @examples
#' data(ammonia)
#' str(ammonia)
#' with(ammonia, plot(corr_ar, perda, xlab = "Air flow", ylab = "Ammonia loss"))
"ammonia"


#' Peripheral Blood Stem Cell (PBSC) Transplant Data
#'
#' Data from a study of 242 patients at the Edmonton Hematopoietic Institute -
#' Alberta Health Services (Espinheira et al., 2026).  The response is the
#' CD34+ cell recovery rate (ratio of viable CD34+ cells post-thaw to
#' pre-freeze), a continuous proportion in (0, 1).
#'
#' @format A data frame with 242 rows and 4 columns:
#' \describe{
#'   \item{\code{recovery}}{CD34+ recovery rate (response), in (0, 1).}
#'   \item{\code{age}}{Patient age (years).}
#'   \item{\code{adj_age}}{Adjusted patient-age variable.}
#'   \item{\code{chemo}}{Chemotherapy protocol indicator:
#'     0 = one-day protocol, 1 = three-day protocol.}
#' }
#'
#' @details
#' The model selected in Espinheira et al. (2026) is:
#' \deqn{\mathrm{logit}(\mu_t) = \beta_1 + \beta_2 x_{t2} + \beta_3 x_{t3}}
#' \deqn{\log(\sigma^2_t) = \gamma_1 + \gamma_2 x_{t2} + \gamma_3 x_{t1}}
#' where \eqn{x_{t1}} = \code{chemo}, \eqn{x_{t2}} = \code{age},
#' \eqn{x_{t3}} = \code{adj_age}.
#'
#' @source Edmonton Hematopoietic Institute - Alberta Health Services.
#'   Used under the data-sharing policy of that institution.
#'
#' @references
#' Espinheira P.L., Silva F.C., Barros M., Ospina R. (2026).
#' A Bootstrap-Calibrated Local Influence Goodness-of-Fit Procedure for
#' Simplex Regression Models.
#'
#' @examples
#' data(pbsc)
#' str(pbsc)
#' hist(pbsc$recovery, breaks = 30, main = "CD34+ Recovery Rate")
"pbsc"
