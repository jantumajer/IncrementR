\name{apical}

\title{Function to calculate series of apical(=height) growth of a tree stem or shrub ramet elongation}

\description{
apical function creates datasets storing total number of tree-rings (measured + missing tree-rings nearby the pith) for individual levels along the stem. If discrepancies in estimated number of tree-rings exist between individual cores coming from the same level and none of the cores includes the pith, median of estimated number of tree-rings is assigned for entire level.
}
  
\usage{
  apical(trw.series, meta, mr.estimate)  
}
  
\arguments{
  \item{trw.series    }{  a dataframe with tree-ring width series as columns and years as rows such as that produced by read.rwl from dplR package}
  \item{meta    } { a dataframe with two columns. The first column (ID) identifies individual trees and coring levels. The other column (Level_cm) contains the distances between stem base and individual coring positions.}
  \item{mr.estimate   }{  a dataframe containing estimated number of missing tree-rings of individual cores such as that produced by EMR}
}
  
\value{
  \item{N.ring_Core   }{  dataframe storing number of measured and estimated missing tree-rings for each  core}
  \item{N.ring_Level    }{  dataframe storing estimated number of total rings for each sampling level}
}
  
\author{ Jan Tumajer, Charles University, Faculty of Science, Department of Physical Geography and Geoecology, Prague, Czech Republic}
  
\examples{
data(IncrementR_data)
  missing.rings <- EMR(trw.series = trw, p.off = po, nyrs = 5, method = "Both")
    # estimate of the number of tree-rings nearby the pith, which were missed by the core
  
  apical.growth <- apical(trw.series = trw, meta = meta, mr.estimate = missing.rings)
    # builds datframe with the total number of tree-rings (measured+missing) in individual coring levels 
  
  drawApicalData(trw.series = trw, apicalData = apical.series, plot = 8, tree = 2)
    # plots the graph of apical growth
}
