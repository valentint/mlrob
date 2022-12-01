######
##  VT::22.11.2022
##
##
##  roxygen2::roxygenise("C:/projects/statproj/R/mlrob", load_code=roxygen2:::load_installed)
##
#'
#'
#'
#' Colon gene expression data (Alon et al., 1999).
#'
#' Gene expression data (2000 genes for 62 samples) from the microarray experiments of
#'  Colon tissue samples of Alon et al. (1999).
#'
#' @name colon
#' @docType data
#' @usage data(colon)
#' @format A data frame with 62 rows and 2001 variables: 2000 gene expressions and one
#'  (the last, \code{Class}) grouping variable: 40 tumor tissues,
#'  coded \code{'t'} and 22 normal tissues, coded \code{'n'}.
#'
#' @source
#'  The data set was taken from the Bioconductor R package \pkg{colonCA}.
#'
#'  Almost identical version of the data set (up to rounding) is available
#'      in the package \pkg{plsgenomics}.
#'
#'  Alon, U. and Barkai, N. and Notterman, D.A. and Gish, K. and Ybarra,
#'      S. and Mack, D. and Levine, A.J. (1999). Broad patterns of gene
#'      expression revealed by clustering analysis of tumor and normal colon
#'      tissues probed by oligonucleotide arrays, Proc. Natl. Acad. Sci. USA, 96(12), 6745--6750.
#'
#' @keywords datasets
NULL
#' Johns Hopkins University Ionosphere database.
#'
#' ''This radar data was collected by a system in Goose Bay, Labrador.  This
#'   system consists of a phased array of 16 high-frequency antennas with a
#'   total transmitted power on the order of 6.4 kilowatts.  The targets
#'   were free electrons in the ionosphere.
#'   "good" radar returns are those showing evidence of some type of structure
#'   in the ionosphere.  "bad" returns are those that do not; their signals pass
#'   through the ionosphere.
#'   Received signals were processed using an autocorrelation function whose
#'   arguments are the time of a pulse and the pulse number.  There were 17
#"   pulse numbers for the Goose Bay system.  Instances in this databse are
#'   described by 2 attributes per pulse number, corresponding to the complex
#'   values returned by the function resulting from the complex electromagnetic
#'   signal.'' [UCI archive]
#'
#' @name ionosphere
#' @docType data
#' @usage data(ionosphere)
#' @format A data frame with 351 rows and 33 variables: 32 measurements and one
#'  (the last, \code{Class}) grouping variable: 225 \code{'good'} and 126 \code{'bad'}.
#'
#'  The original dataset at UCI contains 351 rows and 35 columns. The first 34
#'  columns are features, the last column contains the classification label of
#'  'g' and 'b'. The first feature is binary and the second one is only 0s,
#;  therefore these two features were removed. We remain with 32 featres and
#'  one grouping variable - factor with labels 'good' and 'bad'.
#'
#' @source
#'  Source: Space Physics Group; Applied Physics Laboratory; Johns Hopkins University; Johns Hopkins Road; Laurel; MD 20723
#'
#'  Donor: Vince Sigillito (vgs@aplcen.apl.jhu.edu)
#'
#'  The data have been taken from the UCI Repository Of Machine Learning Databases at
#'  \url{https://archive.ics.uci.edu/ml/datasets/ionosphere}
#'
#'  This data set, with the original 34 features is available in the package \pkg{mlbench}
#'  and a different data set (refering to the same UCI repository) is available in
#'  the package \code{dprep} (archived on CRAN).
#' @references
#'  Sigillito, V. G., Wing, S. P., Hutton, L. V., \& Baker, K. B. (1989).
#'      Classification of radar returns from the ionosphere using neural
#'      networks. Johns Hopkins APL Technical Digest, 10, 262-266.
#' @examples
#'  data(ionosphere)
#'  ionosphere[, 1:6] |> pairs()
NULL
#' BUPA liver disorders.
#'
#' The BUPA Liver Disorders data set was created by BUPA Medical Research and
#'  Development Ltd. (hereafter ''BMRDL'') during the 1980s as part of a larger
#'  health-screening database. At the time the second author was developing
#'  machine learning software, including what may be the first tree-structured
#'  genetic programming (GP) system, and collaborating with the BMRDL
#'  researchers who collected the data. He went on to use the data set
#'  as a GP benchmark. In 1990 the data set was donated on his behalf
#'  to the UCI machine learning repository.
#'  It is hosted at \url{https://archive.ics.uci.edu/ml/datasets/Liver+Disorders}.
#'
#' @name bupa
#' @docType data
#' @usage data(bupa)
#' @format The original data set is a data frame with 345 rows and 7 variables: The first 5
#'  variables are all blood tests which are thought to be sensitive
#'  to liver disorders that might arise from excessive alcohol consumption.
#'  It appears that \code{drinks > 5} is some sort of a selector on this database, i.e.
#'  this could be considered the dependent variable. There is a seventh variable in the
#'  original data set at UCI which was used to split the data set into train and test
#'  for a specific interest. Many studies used actually this variable as a dependent one.
#'  The original data from UCI is stored in the data frame \code{bupa_original} with
#'  345 rows and 7 variables.
#'  The second data frame is \code{bupa} in which the last two variables are
#'  removed and replaced by a variable \code{Class} which is a factor with
#'  levels 1 (for \code{drinks < 5}) and 2 (for \code{drinks >= 5}). There are
#'  four records which are dupkicated and were removed. Thus, the data frame
#'  \code{bupa} is with 341 rows and 6 variables from which the last one is the grouping one.
#'
#' The variables are as follows:
#' \itemize{
#'   \item mcv: 	mean corpuscular volume
#'   \item alkphos: alkaline phosphotase
#'   \item sgpt: alamine aminotransferase
#'   \item sgot: aspartate aminotransferase
#'   \item gammagt: gamma-glutamyl transpeptidase
#'   \item drinks: number of half-pint equivalents of alcoholic beverages drunk per day
#'   \item selector: field used to split data into two sets
#'   \item Class: a factor with \code{1=drinks < 5} and \code{2=drinks >= 5}
#' }
#'
#' @source
#'  -- Creators: BUPA Medical Research Ltd.
#'
#'  -- Donor: Richard S. Forsyth, 8 Grosvenor Avenue, Mapperley Park, Nottingham NG3 5DX, 0602-621676
#'
#'  -- Date: 5/15/1990
#'
#' Available at UCI repository: \url{https://archive.ics.uci.edu/ml/datasets/liver+disorders}
#'
#' @references
#'  James McDermott and Richard S. Forsyth (2016). Diagnosing a disorder in a classification benchmark,
#'  \emph{Pattern Recognition Letters}, \bold{73}, pp 41--43.
#' @examples
#' data(bupa)
#' summary(bupa)
#' table(bupa$Class)
NULL
