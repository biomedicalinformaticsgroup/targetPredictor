\name{getPredictedTargets}
\alias{getPredictedTargets}
\title{Get aggregated ordered list of predicted targets for miRNA}
\usage{
getPredictedTargets(mirna, sources = c("pictar", "diana", "targetscan",
  "miranda"), species = "mmu", min_src = 2, method = "min",
  promote = TRUE, ...)
}
\arguments{
  \item{mirna}{mirna in a standard format}

  \item{sources}{a list of sources to use for aggregation,
  default c('pictar','diana','targetscan','miranda')}

  \item{species}{species in a standard three-letter
  acronym, default 'mmu'}

  \item{min_src}{minimum number of sources required for a
  target to be considered, default 2}

  \item{method}{method of aggregation - choose from
  c('min', 'max', 'geom'), default 'min' is a minimum of
  ranks, 'max' is a maximum of ranks, and 'geom' is based
  on geometric mean of the ranks, it proves to be the most
  accurate.}

  \item{promote}{add weights to improve accuracy of the
  method, default TRUE}

  \item{...}{any optional arguments}
}
\value{
getPredictedTargets a data.frame object where row names are
entrez IDs of target genes
}
\description{
This method performs aggregation of target lists from
multiple sources. Aggregated list is more accurate than any
list from a single source. Multiple aggregation methods are
available.
}
\examples{
##Not run:
targets = getPredictedTargets('let-7a',species='mmu')
head(targets) #top of the list
}
\author{
Maciej Pajak \email{m.pajak@sms.ed.ac.uk}
}
\references{
Friedman, R. C., Farh, K. K.-H., Burge, C. B., and Bartel,
D. P. (2009). Most mammalian mRNAs are conserved targets of
microRNAs. Genome research, 19(1):92-105.

Griffiths-Jones, S., Saini, H. K., van Dongen, S., and
Enright, A. J. (2008). miRBase: tools for microRNA
genomics. Nucleic acids research, 36(Database
issue):D154-8.

Lall, S., Grun, D., Krek, A., Chen, K., Wang, Y.-L., Dewey,
C. N., ... Rajewsky, N. (2006). A genome-wide map of
conserved microRNA targets in C. elegans. Current biology :
CB, 16(5):460-71.

Maragkakis, M., Vergoulis, T., Alexiou, P., Reczko, M.,
Plomaritou, K., Gousis, M., ... Hatzigeorgiou, A. G.
(2011). DIANA-microT Web server upgrade supports Fly and
Worm miRNA target prediction and bibliographic miRNA to
disease association. Nucleic Acids Research, 39(Web Server
issue), W145-8.
}

