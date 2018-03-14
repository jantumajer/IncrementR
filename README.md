# Function syntax and descriptions

```Eccentricity(trw.series, complete=FALSE)```

> Function `Eccentricity` calculates eccentricity indices according to Braam et al. (1987), Alestalo (1971) and Schweingruber (1996) for each tree-ring. This function requires the input of four cross-dated series per tree and/or sampling level (argument `trw.series`). Argument complete defines whether four samples were obtained for all sampling levels. If some of required samples are missing, then `complete=F`. The output of this function is a list containing three data frames with values of eccentricity indices for each method and each orientation around the stem circumference.

```EMR(trw.series, nyrs, p.off, method="TRW")```

Function `EMR` estimates the number of missing tree-rings between the last measured tree-ring and the pith of each core. Argument `trw.series` defines data frame with series for which the estimation shall be performed. Argument `nyrs` defines the number of innermost measured tree-rings, whose average will be used to estimate mean tree-ring width of missing rings. Argument `p.off` contains data frame with an information about estimated distance in mm of the last known tree-ring to the pith (stored in the column with name “P.OFFSET”). Last parameter is `method`, which defines the method used for estimation. User can select from `“TRW”`, `“BAI”` and `“Both”`, which is and average value of `“BAI”` and `“TRW”` (for detailed description of the approaches see Altman et al. 2016). The output of this function is a data frame with the number of missing tree-rings that should be added to each core.

```RMR(trw.series, meta, mr.estimate, nyrs=5, nsph=4)```

Function `RMR` may be applied following the `EMR` processing, and serves to add the number of missing rings to pith into data frame with tree-ring width series. Arguments of `RMR` function include `trw.series` (contains tree-ring series to be modified), `mr.etimate` (an output of `EMR` function for trw.series), `nyrs` (the number of innermost measured tree-rings, whose average will be used to estimate mean tree-ring width of missing rings) and `nsph`, which defines how many cores were sampled at each level (options are 1, 2 and 4). If `nsph`>1, the estimation of the number of missing rings from `EMR` function is corrected so as to have the same number of tree-rings in all series coming from the same level (for the description of approach see `apical` function). The output of this function is a new data frame, which contains tree ring-series with added “virtual” tree-rings between pith and the first measured tree-ring.

```apical(trw.series, meta, mr.estimate)```

The function `apical` first creates datasets containing number of TotalRings (Measured + Missing rings nearby pith as derived by EMR) for each core and each level along the stem. Inputs include `trw.series`, `meta` file with sampling heights of each sample in cm (column “Level_cm”) and `mr.estimate` (the output of `EMR` function). The output is the list of two data frames; the first one (`N.ring_Core`) contains series codes with the number of measured, missing and total rings. The other one (`N.ring_Level`) contains codes of each sampling level with the number of measured, missing and total rings, together with the level height. Moreover, the mean height growth velocity between two subsequent levels (in cm.yr-1) is calculated and stored in “Speed.cmyr” column and error of its estimation (see Figure 1 for description) is stored in “MeanSpeedError.cmyr” column.
If the number of estimated tree-rings differs from the other cores of the same height level, two rules are applied to estimate number of tree-rings in a given level. First, cores with pith presence are prioritized; and if the number of tree-rings of such cores differs, median of their TotalRings is used to estimate the number of tree-rings in height level. Secondly, if none of cores contains the pith, median of TotalRings of all series from the same level is used.

```taperCalcul(trw.series, meta)```

Function `taperCalcul` calculates taper and taper angle between subsequent sampling heights for all trees. The inputs are `trw.series` (containing tree-ring series) and `meta` (metadata containing the sampling heights in cm in column “Level_cm”). 

```BAIcalculation(trw.series)```

Function `BAIcalculation` creates a series of basal area increments for each sampling level. It requires ideally four cores per one sampling point, however, calculation will run also with 3 or 2 cores per level. The main difference of this function comparing functions `bai.in` and `bai.out` (dplR package) is consideration of eccentric tree-rings. The only input argument of this function is `trw.series` (containing cross dated tree-ring series and, optionally, replaced missing tree-rings nearby the pith using `EMR` and `RMR`). The output of this function is a data frame containing series of basal area increments for each sampling level.

```drawEccentricityGraph (trw, ecc, meta, plot, tree, withEccentricity=T, method="Schweingruber")```

First graphical function of “IncrementR” package is `drawEccentricityGraph`. Argument `trw` represents a data frame with tree-ring series (possibly) previously corrected for missing rings nearby the pith using `RMR`. Argument `ecc` contains output of eccentricity function.  Arguments `plot` and `tree` specify a tree which will be drawn, argument `meta` refers to a metadata file with sampling heights of each level in cm (column “Level_cm”). Logical argument `withEccentricity` specifies whether the graph will show only stem profile or stem profile with eccentricity indices; and, finally, `method` defines index to be drawn (available options include `“Alestalo”`, `“Braam”` and `“Schweingruber”`). The output consists of two graphs of stem profiles (E-W and S-N) with respective eccentricity indices of each sampling level.

```drawCrossSectionProfile (trw.series, plot, tree, level, show.legend=T)```

Cross section profile based on approximation of cross section as an ellipse can be drawn for each sampling level. The function `drawCrossSectionProfile` has arguments specifying the cross section (`plot`, `tree` and `level`) to be plot; `trw.series` (data frame containing tree-ring series, possibly previously corrected using `RMR`); and logical argument `show.legend` specifying if the legend (colours of different calendar years) will be drawn.

```drawBAI(baiFile, plot=1,tree=1,logscale=F, show.legend=T)```

The output of function `BAIcaclulation` can be drawn using `drawBAI` function as a line-chart. Arguments of this function include `baiFile` (output of `BAIcalculation` function); arguments `plot` and `tree` specify a tree to be plot. Logical argument `logscale` specifies whether the result should be presented on a log scaled y-axis; and logical argument `show.legend` specifies whether the legend will be drawn.

```drawTaper(taperFile, plot=1,tree=1, variant="Taper")```

Taper between consecutive sampling levels can be plot by `taperCalculation` function as a line-chart. Argument `taperFile` is the output of `taperCalculation` function; arguments `plot` and `tree` specify tree to have the taper drawn; and logical argument `variant` specifies whether taper or taper angle will be drawn (available options are `“Taper”` or `“Angle”`).

```drawApicalData(apicalData, plot=1, tree=1)```

Arguments of `drawApicalData` function include `apicalData` (output of `apical` function); and `plot` and `tree` (specifying the tree to be visualized). The resulting plot consist of two parts – the left panel shows the number of “Total” tree-rings (i.e., measured tree-rings + estimated missing tree-rings nearby pith) in each sampling level along the stem; the right one shows mean height growth velocity between two consecutive sampling levels (points) together with mean errors of its estimates (error bars). 
