\name{fluA.inc.per}
\alias{fluA.inc.per}
\docType{data}
\title{
Coarse incubation period data for influenza A
}
\description{
These observations on the incubation period of influenza A come from a variety of sources, and were gathered for a literature review.  They report doubly interval-censored, single interval-censored or exact observations for the incubation period. 
}
\usage{data(fluA.inc.per)}
\format{
  A data frame with 151 observations on the following 7 variables.
  \describe{
    \item{\code{author}}{the name of the primary author for the source of the observation}
    \item{\code{year}}{the year of the study which is the source of the observation}
    \item{\code{EL}}{the earliest possible time of infection}
    \item{\code{ER}}{the latest possible time of infection}
    \item{\code{SL}}{the earliest possible time of symptom onset}
    \item{\code{SR}}{the latest possible time of symptom onset}
    \item{\code{type}}{an indicator of the type of observation: 0 for doubly interval-censored, 1 for single-interval censored, 2 for exact}
  }
}
\source{
Lessler J, Reich NG, Brookmeyer R, Perl TM, Nelson KE, Cummings DAT. (2009) A systematic review of the incubation periods of acute respiratory viral infections. Lancet Infectious Diseases. 9(5):291-300.
}
\examples{
data(fluA.inc.per)
head(fluA.inc.per)
}
\keyword{datasets}
