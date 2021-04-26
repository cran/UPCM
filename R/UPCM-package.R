

#' Uncertainty in Partial Credit Models
#' 
#' Performs UPCM, a method to model uncertainty in (Generalized) Partial Credit Models
#' 
#' 
#' @name UPCM-package
#' @docType package
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}\cr
#' \url{https://www.sg.tum.de/epidemiologie/team/schauberger/}
#' @seealso \code{\link{UPCM}}
#' @keywords package Partial Credit Model Uncertainty UPCM
#' @references Tutz, Gerhard and Schauberger, Gunther (2020): Uncertainty in Latent Trait Models, 
#' \emph{Applied Psychological Measurement}, \url{https://journals.sagepub.com/doi/abs/10.1177/0146621620920932?journalCode=apma}
#' @examples
#' \donttest{
#' data(tenseness)
#' 
#' Y <- data.matrix(tenseness[,1:4])
#' X <- model.matrix(~ Gender + Age, data = tenseness)[,-1]
#' 
#' m_upcm <- UPCM(Y = Y, X = X, cores = 2, GPCM = FALSE)
#' m_upcm
#' plot(m_upcm)
#' }
#' \dontshow{
#' set.seed(1860)
#' n <- 50
#' I <- 2
#' Y <- matrix(sample(1:3, I*n, replace = TRUE), ncol = I)
#' m_upcm <- UPCM(Y = Y, cores = 1, GPCM = FALSE, se = FALSE, ctrl.nlminb = list(rel.tol = 1e-06))
#' m_upcm
#' }
NULL

#' Tenseness data from the Freiburg Complaint Checklist
#' 
#' Data from the Freiburg Complaint Checklist. 
#' The data contain all 8 items corresponding to the scale \emph{Tenseness} for 2042 participants of the 
#' standardization sample of the Freiburg Complaint Checklist. 
#' 
#' @name tenseness
#' @docType data
#' @format A data frame containing data from the Freiburg Complaint Checklist with 1847 observations. 
#' All items refer to the scale \emph{Tenseness} and are measured on a 5-point Likert scale where low numbers 
#' correspond to low frequencies or low intensitites of the respective complaint and vice versa. 
#' \describe{ 
#' \item{Clammy_hands}{Do you have clammy hands?}
#' \item{Sweat_attacks}{Do you have sudden attacks of sweating?}
#' \item{Clumsiness}{Do you notice that you behave clumsy?}
#' \item{Wavering_hands}{Are your hands wavering frequently, e.g. when lightning a cigarette or when holding a cup?}
#' \item{Restless_hands}{Do you notice that your hands are restless?}
#' \item{Restless_feet}{Do you notice that your feet are restless?}
#' \item{Twitching_eyes}{Do you notice unvoluntary twitching of your eyes?}
#' \item{Twitching_mouth}{Do you notice unvoluntary twitching of your mouth?}
#' \item{Gender}{Gender of the person}
#' \item{Household}{Does the person live alone in a household or together with somebody?}
#' \item{Income}{Income, categorized to levels from 1 (low income) to 11(high income). For simplicity,
#' due to the high number of categories income can be treated as a metric variable.}
#' \item{WestEast}{Is the person from East Germany (former GDR)?}
#' \item{Abitur}{Does the person have Abitur (A-levels)?}
#' \item{Age}{Age of the person}
#'  }
#' @source 
#' ZPID (2013). PsychData of the Leibniz Institute for Psychology Information ZPID. Trier: Center for Research Data in Psychology.
#' 
#' Fahrenberg, J. (2010). Freiburg Complaint Checklist [Freiburger Beschwerdenliste (FBL)]. Goettingen, Hogrefe.
#' @keywords datasets
#' @examples
#' 
#' data(tenseness)
#' 
NULL
