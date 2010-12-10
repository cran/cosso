\name{ozone}
\docType{data}
\alias{ozone}
\title{Ozone pollution data in Los Angels, 1976}
\description{
This is the ozone data used in Breiman and Friedman (1985). This dataset contains 330 observations, 
and each observation is a daily measurement.
}

\usage{
data(ozone)
}

\format{
\tabular{ll}{
ozone \tab Ozone reading\cr
temp \tab Temperature (degree F)\cr
invHt \tab Inversion base height (feet) \cr
press \tab Pressure gradient (mm Hg)\cr
vis  \tab Visibility (miles)\cr
milPress \tab 500 millibar pressure height (m)\cr
hum \tab Humidity (percent) \cr
invTemp \tab Inversion base temperature (degrees F)\cr
wind \tab Wind speed (mph) \cr
}
}

\details{
All predictors have been scaled to range from 0 to 1.
}

\source{
Breiman, L and Friedman, J. (1985), "Estimating Optimal Transformations for Multiple Regression
and Correlation", JASA, \bold{80}, 580--598.
}
