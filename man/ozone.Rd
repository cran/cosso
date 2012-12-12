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
ozone   \tab Ozone reading\cr
temp    \tab Temperature (degree C). minimum 25 and maximum 93 in original scale. \cr
invHt   \tab Inversion base height (feet). minimum 111 and maximum 5000 in original scale.  \cr
press   \tab Pressure gradient (mm Hg). minimum -69 and maximum 107 in original scale. \cr
vis     \tab Visibility (miles). minimum 0 and maximum 350 in original scale. \cr
milPress\tab 500 millibar pressure height (m). minimum 5320 and maximum 5950 in original scale. \cr
hum     \tab Humidity (percent). minimum 19 and maximum 93.  \cr
invTemp \tab Inversion base temperature (degrees F). minimum -25 and maximum 332 in original scale. \cr
wind    \tab Wind speed (mph). minimum 0 and maximum 21 in original scale.  \cr
}
}

\details{
All the variables, except for the response, have been sacled to [0,1] interval. To transform
back to the original scale, use the formula:
\deqn{x=min+(max-min)*z.}
}

\source{
Breiman, L. and Friedman, J. (1985), "Estimating Optimal Transformations for Multiple Regression
and Correlation", Journal of the American Statistical Association, \bold{80}, 580--598.
}
