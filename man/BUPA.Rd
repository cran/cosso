\name{BUPA}
\docType{data}
\alias{BUPA}
\title{BUPA Liver Disorder Data}
\description{
345 male patients' blood test result and liver disorder status.
}

\usage{
data(BUPA)
}

\format{
\tabular{ll}{
CLASS   \tab  0: no liver disorder 1: liver disorder\cr
MCV     \tab  mean corpuscular volume. minimum 65 and maximum 103 in original scale.    \cr
ALKPHOS \tab  alkaline phosphotase. minimum 23 and maximum 138 in original scale.       \cr
SGPT    \tab 	alamine aminotransferase. minimum 4 and maximum 155 in original scale.   \cr
SGOT    \tab  aspartate aminotransferase. minimum 5 and maximum 82 in original scale. \cr
GAMMAGT \tab  gamma-glutamyl transpeptidase. minimum 5 and maximum 297 in original scale. \cr
DRINKS  \tab  number of alcoholic beverages drunk per day. minimum 0 and maximum 20 in original scale. \cr
}
}

\details{
All the variables, except for the response, have been sacled to [0,1] interval. To transform
back to the original scale, use the formula:
\deqn{x=min+(max-min)*z.}
}


\source{
Richard S. Forsyth at BUPA Medical Research Ltd.
}
