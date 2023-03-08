\name{veteran}
\docType{data}
\alias{veteran}
\title{Veterans' Administration Lung Cancer study}
\description{
Randomized trial of two treatment regimens for lung cancer.
}

\usage{
data(veteran)
}

\format{
\tabular{ll}{
time    \tab survival time      \cr
status  \tab censoring status   \cr
trt     \tab 0=standard 1=test  \cr
celltype\tab 1=squamous, 2=smallcell, 3=adeno, 4=large. \cr
karno   \tab Karnofsky performance score. minimum 10 and maximum 99 in original scale. \cr
diagtime\tab months from diagnosis to randomization. minimum 1 and maximum 87 in original scale. \cr
age     \tab in years. minimum 34 and maximum 81 in original scale. \cr
prior   \tab prior therapy 0=no, 1=yes.
}
}

\details{
All the variables, except for the response, have been scaled to [0,1] interval. To transform
back to the original scale, use the formula:
\deqn{x=min+(max-min)*z.}
}


\source{
Kalbfleisch, J. and Prentice, R.L. (2002), The Statistical Analysis of Failure Time Data (Second Edition) Wiley: New Jersey.
}
