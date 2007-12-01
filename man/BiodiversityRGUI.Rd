\name{BiodiversityRGUI}
\alias{BiodiversityRGUI}

\title{GUI for Biodiversity Analysis and Ordination }

\description{
This function provides a GUI (Graphical User Interface) for some of the functions 
of \pkg{vegan}, some other packages and some new functions to run biodiversity 
analysis, including species accumulation curves, diversity indices, Renyi 
profiles, rank-abundance curves, GLMs for analysis of species abundance and
presence-absence, distance matrices, Mantel tests, cluster and ordination 
analysis (including constrained ordination methods such as RDA, CCA, db-RDA and
CAP). The function depends and builds on \pkg{Rcmdr}, performing all analyses
on the community and environmental datasets that the user selects. A thorough
description of the package and the biodiversity and ecological methods that it
accomodates (including examples) is provided in the freely available Tree Diversity Analysis manual (Kindt and Coe, 2005).
}

\usage{
BiodiversityRGUI()
}

\details{
The function launches the R-Commander GUI with an extra menu list for
common statistical methods for biodiversity and community ecology analysis.
 
The R-Commander is launched by changing the location of the Rcmdr "etc" folder to the "etc" folder of BiodiversityR. As the files of the "etc" folder of BiodiversityR are copied from Rcmdr 1.3-9, it is possible that newer versions of the R-Commander will not be launched properly. In such situations, it is possible that copying all files from the Rcmdr "etc" folder again and adding the BiodiversityR menu options to the Rcmdr-menus.txt is all that is needed to launch the R-Commander again.
 
BiodiversityR uses two data sets for analysis: the community dataset (or community matrix or species matrix) and the environmental dataset (or environmental matrix). The environmental dataset is the same dataset that is used as the "active dataset" of The R-Commander. (Note that you could sometimes use the same dataset as both the community and environmental dataset. For example, you could use the community dataset as environmental dataset as well to add information about specific species to ordination diagrams. As another example, you could use the environmental dataset as community dataset if you first calculated species richness of each site, saved this information in the environmental dataset, and then use species richness as response variable in a regression analysis.) Some options of analysis of ecological distance allow the community matrix to be a distance matrix (the community data set will be interpreted as distance matrix via \code{\link[stats:dist]{as.dist}} prior to further analysis). 
 
BiodiversityR provides the following menu options (each described below in greater detail):

\itemize{
  \item{\bold{Select community dataset} (Community matrix menu)}{
Selects a dataset to be the community dataset.
} 
 
  \item{\bold{Import datasets from Excel} (Community matrix menu)}{
Imports a community and environmental dataset from an Excel workbook.
} 

  \item{\bold{Import datasets from Access} (Community matrix menu)}{
Imports a community and environmental dataset from an Access database.
} 

  \item{\bold{View community data set} (Community matrix menu)}{
Invoke the R text editor to view the data of the community data set.
} 

  \item{\bold{Edit community data set} (Community matrix menu)}{
Invoke the R text editor to edit the data of the community data set.
} 

  \item{\bold{Check data sets} (Community matrix menu)}{
Check whether the community and environmental data sets have compatible dimensions.
} 

  \item{\bold{Same sites for community and environmental} (Community matrix menu)}{
Creates a new community dataset with the same sites sequence as the environmental matrix.
} 

  \item{\bold{Make community dataset} (Community matrix menu)}{
Creates a community dataset from the environmental dataset.
} 

  \item{\bold{Remove NA} (Community matrix menu)}{
Removes the same sites with NA from the environmental and community datasets.
} 

  \item{\bold{Transform community matrix} (Community matrix menu)}{
Transforms the community matrix.
} 

  \item{\bold{Select environmental data set} (Environmental matrix menu)}{
Selects a dataset to be the environmental dataset.
} 

  \item{\bold{View environmental data set} (Environmental matrix menu)}{
Invoke the R text editor to view the data of the environmental dataset.
} 

  \item{\bold{Edit environmental data set} (Environmental matrix menu)}{
Invoke the R text editor to edit the data of the environmental dataset.
} 

  \item{\bold{Summary} (Environmental matrix menu)}{
Explores variables of the environmental dataset.
} 

  \item{\bold{Box Cox transformation} (Environmental matrix menu)}{
Creates a transformed variable from one of the variables of the environmental dataset.
} 

  \item{\bold{Species accumulation curves} (Analysis of diversity menu)}{
Estimates and plots species accumulation curves.
} 

  \item{\bold{Diversity indices} (Analysis of diversity menu)}{
Calculates and plots diversity indices.
} 

  \item{\bold{Rank abundance} (Analysis of diversity menu)}{
Calculates and plots rank-abundance curves.
} 

  \item{\bold{Renyi profile} (Analysis of diversity menu)}{
Calculates and plots Renyi diversity profiles.
} 

  \item{\bold{Species abundance as response} (Analysis of species as response menu)}{
Fits and plots regression models assuming that the response variable is count data.
} 

  \item{\bold{Species presence-absence as response} (Analysis of species as response menu)}{
Fits and plots regression models transforming and analysing the response variable as presence-absence.
} 

  \item{\bold{Calculate distance matrix} (Analysis of ecological distance menu)}{
Calculates a distance matrix.} 

  \item{\bold{Unconstrained ordination} (Analysis of ecological distance menu)}{
Fits and plots unconstrained ordination models.
} 

  \item{\bold{Constrained ordination} (Analysis of ecological distance menu)}{
Fits and plots constrained ordination models.
} 

  \item{\bold{Clustering} (Analysis of ecological distance menu)}{
Calculates and plots results from clustering algorithms.
} 

  \item{\bold{Compare distance matrices} (Analysis of ecological distance menu)}{
Conducts some analysis such as Mantel, MRPP and ANOSIM tests on distance matrices.
} 

  \item{\bold{Help about BiodiversityR} (Help menu)}{
Opens the help file available for the BiodiversityR package (including this html file).
} 

  \item{\bold{Citations for loaded packages} (Help menu)}{
Provides a list of all the loaded packages and gives citation information.
} 

  \item{\bold{Go to website for BiodiversityR} (Help menu)}{
Links to the website for the BiodiversityR package and Tree Diversity Analysis manual.
} 

  \item{\bold{Tree diversity analysis manual} (Help menu)}{
Links to the PDF version of the Tree Diversity Analysis manual. Separate chapters can be downloaded from the website of BiodiversityR (see directly above).
} 
}
}


\section{Select Community Dataset}{

This window selects the community dataset to be used in the biodiversity analyses and provides the following options:

\itemize{
  \item{\bold{Data Sets (pick one)}}{
A drop-down list is provided with all the datasets that are available. The current community data set is indicated, or the first data set of the list is shown.
New datasets can be loaded through the Data menu of the Rcmdr or through the "import from Excel" option of BiodiversityR.
} 

  \item{\bold{OK}}{
Make the selected data set the community data set.
} 

  \item{\bold{Cancel}}{
Close the window and do not select a new data set.
}
}
}


\section{Import datasets from Excel}{

This window enables to import community and environmental data sets from an Excel workbook with a specific format (sheets that are named "community" and "environmental" or "stacked" and "environmental"; first row containing the names of the variables) (an example is provided in the etc folder of the BiodiversityR package for the dune meadow dataset: dune.xls). The menu provides the following options:

\itemize{
  \item{\bold{Enter name for community data set}}{
The name for the new community dataset.
} 

  \item{\bold{Enter name for environmental data set}}{
The name for the new environmental dataset.
} 

  \item{\bold{Enter name for variable for sites}}{
The name for the variable that indicates site names in the new dataset. The same variable should be available from the various sheets.
Passed as argument for "sitenames" for function \code{\link{import.from.Excel}}
} 

  \item{\bold{Import community dataset from stacked format}}{
Import the community data set from the stacked format or not. The stacked format is the only possibility for community data with more than 255 species (related to the maximum number of columns in Excel).
  
Option "Yes" will result in argument of "stacked" for "sheet" for function \code{\link{import.from.Excel}} to import the community data set.
 
Option "No" will result in argument of "community" for "sheet" for function \code{\link{import.from.Excel}} to import the community data set.
 
Either option will result in argument of "environmental" for "sheet" for function \code{\link{import.from.Excel}} to import the environmental data set.

} 

  \item{\bold{Enter variable for species}}{
This option is only available for stacked data. The list shows the variables that can be used for the names of species (shown as names for the columns). 
  
Passed as argument for "column" of function \code{\link{import.from.Excel}}.
} 

  \item{\bold{Enter variable for abundance}}{
This option is only available for stacked data. The list shows the variables that can be used for the abundance values (shown as totals for cells). 
Passed as argument for "value" of function \code{\link{import.from.Excel}}.
} 

  \item{\bold{Enter factor for subset}}{
This option is only available for stacked data. The list shows the variables that can be used for the abundance values (shown as totals for cells). 
Passed as argument for "factor" of function \code{\link{import.from.Excel}}.
} 

  \item{\bold{Enter level for subset}}{
Chooses the value for the subset variable to create the subset. 
Passed as argument for "level" of function \code{\link{import.from.Excel}}.
} 

  \item{\bold{OK}}{
Import the community and environmental datasets and make these the active datasets.
} 

  \item{\bold{Cancel}}{
Close the window and do not import new datasets.
} 
}
}


\section{Import datasets from Access}{

This window enables to import community and environmental data sets from an Access database with a specific format (sheets that are named "community" and "environmental" or "stacked" and "environmental"; first row containing the names of the variables) (an example is provided in the etc folder of the BiodiversityR package for the dune meadow dataset: dune.xls). The menu provides the following options:

\itemize{
  \item{\bold{Enter name for community data set}}{
The name for the new community dataset.
} 

  \item{\bold{Enter name for environmental data set}}{
The name for the new environmental dataset.
} 

  \item{\bold{Enter name for variable for sites}}{
The name for the variable that indicates site names in the new dataset. The same variable should be available from the various sheets.
Passed as argument for "sitenames" for function \code{\link{import.from.Access}}
} 

  \item{\bold{Import community dataset from stacked format}}{
Import the community data set from the stacked format or not.
  
Option "Yes" will result in argument of "stacked" for "sheet" for function \code{\link{import.from.Access}} to import the community data set.
 
Option "No" will result in argument of "community" for "sheet" for function \code{\link{import.from.Access}} to import the community data set.
 
Either option will result in argument of "environmental" for "sheet" for function \code{\link{import.from.Access}} to import the environmental data set.

} 

  \item{\bold{Enter variable for species}}{
This option is only available for stacked data. The list shows the variables that can be used for the names of species (shown as names for the columns). 
  
Passed as argument for "column" of function \code{\link{import.from.Access}}.
} 

  \item{\bold{Enter variable for abundance}}{
This option is only available for stacked data. The list shows the variables that can be used for the abundance values (shown as totals for cells). 
Passed as argument for "value" of function \code{\link{import.from.Access}}.
} 

  \item{\bold{Enter factor for subset}}{
This option is only available for stacked data. The list shows the variables that can be used for the abundance values (shown as totals for cells). 
Passed as argument for "factor" of function \code{\link{import.from.Access}}.
} 

  \item{\bold{Enter level for subset}}{
Chooses the value for the subset variable to create the subset. 
Passed as argument for "level" of function \code{\link{import.from.Access}}.
} 

  \item{\bold{OK}}{
Import the community and environmental datasets and make these the active datasets.
} 

  \item{\bold{Cancel}}{
Close the window and do not import new datasets.
} 
}
}


\section{Same sites for community and environmental datasets}{

This window maps the community dataset onto the rownames of the environmental dataset by function \code{\link{same.sites}}. Having the same sequence of sites is an assumption for analysis with BiodiversityR. It may be useful to use this function after making a community dataset from a stacked environmental dataset (especially as sites are ordered in an alphabetic way from the stacked dataset, which may create problems with X1, X10, X100 site names versus the X001, X010 and X100 formats; the function is also useful where some sites do not contain any species). The menu provides the following options:

\itemize{
  \item{\bold{save original community matrix}}{
If this option is selected, the original data set is saved under the name of the community dataset followed by ".orig".
} 

  \item{\bold{OK}}{
Order the sites of the community dataset in exactly the same way as the sites of the environmental data set, leaving out sites that do not have matching names in the environmental data set.
} 

  \item{\bold{Cancel}}{
Close the window and do not re-order and select the sites.
} 
}
}


\section{Make Community Dataset}{

This window selects the variables that indicates sites, species and abundance to create a new community dataset. This dataset becomes the active community dataset. The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new community dataset.
} 

  \item{\bold{Site variable (rows)}}{
The list shows the variables that can be used for the names of sites (shown as names for the rows).  
Passed as argument for "row" of function \code{\link{makecommunitydataset}}.
} 

  \item{\bold{Species variable (columns)}}{
The list shows the variables that can be used for the names of species (shown as names for the columns).  
Passed as argument for "column" of function \code{\link{makecommunitydataset}}.
} 

  \item{\bold{Abundance variable}}{
The list shows the variables that can be used for the abundance values (shown as totals for cells). 
Passed as argument for "value" of function \code{\link{makecommunitydataset}}.
} 

  \item{\bold{Subset options}}{
The list shows the variables that can be used for the abundance values (shown as totals for cells). 
Passed as argument for "factor" of function \code{\link{makecommunitydataset}}.
} 

  \item{\bold{Subset}}{
Chooses the value for the subset variable to create the subset. 
Passed as argument for "level" of function \code{\link{makecommunitydataset}}.
} 

  \item{\bold{OK}}{
Create the community data set and make it the active community dataset.
} 

  \item{\bold{Cancel}}{
Close the window and do not create a new community dataset.
} 
}
}

\section{Remove NA}{

This window removes the sites that have \code{NA} (missing values) for a selected varialbe of the environmental dataset. When environmental variables have missing values, this often creates problems with biodiversity analysis. The menu provides the following options:

\itemize{
  \item{\bold{Select variable}}{
The list shows the variables that can be used to remove sites with NA.  
Passed as argument for var for functions \code{\link{removeNAcomm}} and \code{\link{removeNAenv}}.
} 

  \item{\bold{OK}}{
Remove the sites with \code{NA}.
} 

  \item{\bold{Cancel}}{
Close the window and do not remove the sites with \code{NA}.
} 
}
}

\section{Transform community matrix}{

This window transforms the community matrix. The menu provides the following options:

\itemize{
  \item{\bold{Method}}{
Method of transforming the community dataset. 
Passed as argument for "method" for function \code{\link{disttransform}}. 
The transformed community matrix is saved under the same name of the original dataset, and the current community dataset therefore becomes the transformed community dataset.
} 

  \item{\bold{Save original community matrix}}{
This option saves the untransformed community dataset by adding .orig to the name of the community dataset, as the function replaces the original dataset with the transformed community dataset.
} 

  \item{\bold{OK}}{
Calculate the new community matrix.
} 

  \item{\bold{Cancel}}{
Close the window and do not calculate a new community matrix.
} 
}
}

\section{Select Environmental Dataset}{

This window selects the environmental dataset to be used in the biodiversity analyses. The environmental dataset is always the active dataset for non-Biodiversity Rcmdr options. By selecting the community dataset as the environmental dataset as well, you can also manipulate the community dataset with the other Rcmdr options. The menu provides the following options:

\itemize{
  \item{\bold{Data Sets (pick one)}}{
A drop-down list is provided with all the datasets that are available. The current community data set is indicated, or the first data set of the list is shown.
New datasets can be loaded through the Data menu of the Rcmdr or through the "import from Excel" option of BiodiversityR.
} 

  \item{\bold{OK}}{
Make the selected data set the environmental data set.
} 

  \item{\bold{Cancel}}{
Close the window and do not select a new data set.
} 
}
}

\section{Summary}{

This window makes a summary of all or a selection of the variables of the environmental dataset, or plots the variables. In case that you want to make a summary of the community dataset, then you need to make the community dataset the environmental dataset at the same time. The menu provides the following options:

\itemize{
  \item{\bold{Select variable}}{
A drop-down list is provided with all the variables of the environmental dataset. 
The first item of the list (all) is reserved to make a summary of all variables. datasets that are available. 
} 

  \item{\bold{OK}}{
Make a summary of all variables or the selected variable by function \code{\link[base]{summary}}.
} 

  \item{\bold{Plot}}{
Plots all variables against each other with function \code{\link[graphics]{pairs}}, plots a selected continuous variable with function \code{\link[graphics]{plot}} or plots a categorical with function \code{\link[graphics]{boxplot}}.
} 

  \item{\bold{Cancel}}{
Close the window and do provide any summary or plot.
} 
}
}

\section{Box Cox transformation}{

This window makes a Box-Cox transformation of a selected variable from the environmental dataset. The menu provides the following options:

\itemize{
  \item{\bold{Select variable}}{
A drop-down list is provided with all the variables of the environmental dataset.  Click on the variable to transform.
} 

  \item{\bold{OK}}{
Calculates a Box-Cox transformation of the selected variable with function \code{\link[car]{box.cox.powers}}. Makes a QQ-plot (function \code{\link[car]{qq.plot}}), and performs a Shapiro test (function \code{\link[stats]{shapiro.test}}) and Kolmogorov-Smirnov test (function \code{\link[stats]{ks.test}}) of the original and transformed variable.
} 

  \item{\bold{Cancel}}{
Close the window.
} 
}
}

\section{Species accumulation curves}{

This window fits and plots species accumulation curves. The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the estimated species accumulation curve after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{Accumulation method}}{
Select the method of species accumulation. 
Passed as argument for "method" of functions \code{\link{accumresult}} or \code{\link{accumcomp}}.
} 

  \item{\bold{permutations}}{
Number of permutations for random species accumulation.  
Passed as argument for "permutation" of functions \code{\link{accumresult}} or \code{\link{accumcomp}}.
} 

  \item{\bold{scale of x axis}}{
Method of scaling the horizontal axis. 
Passed as argument for "scale" of functions \code{\link{accumresult}} or \code{\link{accumcomp}}.
} 

  \item{\bold{subset options}}{
The list shows the variables that can be used for selecting subsets. Option "all" indicates that no subset will be calculated. 
In case a variable is selected, it will be passed as argument for "factor" of functions \code{\link{accumresult}} or \code{\link{accumcomp}}.
} 

  \item{\bold{Subset}}{
Subset chooses which subsets are calculated. 
In case that the value of "." (a period) is selected then function \code{\link{accumcomp}} will used to calculate the species accumulation curve and to plot the curve (you may need to click in the graph to show where the legend needs to be placed). 
In case another value is chosen, then this will be the argument for "level" of function \code{\link{accumresult}}.
} 

  \item{\bold{Plot options}}{
Options for plotting passed to function \code{\link{accumplot}}.
 
Option "addplot" sets "addit=T" meaning that the species accumulation curve will be added to an existing graph.
 
Option "x limits"sets "xlim". Providing "1,10" will plot between 1 and 10.
 
Option "y limits"sets "ylim". Providing "2,20" will plot between 2 and 20.
 
Option "ci"sets "ci".

Option "symbol"sets "pch".

Option "cex"sets "cex".

Option "colour" sets "col".
} 

  \item{\bold{OK}}{
Calculate the species accumulation curve with functions functions \code{\link{accumresult}} or \code{\link{accumcomp}}.
} 

  \item{\bold{Plot}}{
Plot the species accumulation curve with the name listed on top with function  \code{\link{accumplot}}. You may need to click in the graph to indicate where the legend needs to be placed.
} 

  \item{\bold{Cancel}}{
Close the window and do not calculate a new species accumulation curve.
} 
}
}

\section{Diversity indices}{

The window calculates and fits diversity indices from the community dataset. The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the estimated diversity indices after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you plot the result by typing in the name of previous result first in this box. To obtain a meaningful graph, you need to provide similar selections as for the original result (and it may thus be easier to recalculate first and then plot immediately).
} 

  \item{\bold{Diversity index}}{
Select the diversity index. 
Passed as argument for "index" of functions \code{\link{diversityresult}} or \code{\link{diversitycomp}}.
} 

  \item{\bold{Calculation method}}{
Select the method of calculation. 
Passed as argument for "method" of functions \code{\link{diversityresult}} or \code{\link{diversitycomp}}.
} 

  \item{\bold{subset options}}{
The list shows the variables that can be used for selecting subsets. Option "all" indicates that no subset will be calculated. 
In case a variable is selected, it will be passed as argument for "factor" of functions \code{\link{diversityresult}} or \code{\link{diversitycomp}}.
} 

  \item{\bold{Subset}}{
Subset chooses which subsets are calculated. 
In case that the value of "." (a period) is selected then function \code{\link{diversitycomp}} will used to calculate the species accumulation curve and to plot the curve (you may need to click in the graph to show where the legend needs to be placed). 
In case another value is chosen, then this will be the argument for "level" of function \code{\link{diversityresult}}.
} 

  \item{\bold{Output options}}{
Options for obtaining results with functions \code{\link{diversityresult}}, \code{\link{diversitycomp}} or for plotting results.
 
Option "save results" results in adding a new variable with the diversity indices to the environmental dataset. This method only works for calculation method "separate per site" and function \code{\link{diversityresult}}.
 
Option "sort results" results in setting option "sortit=T" for functions \code{\link{diversityresult}} or \code{\link{diversitycomp}}.
 
Option "label results" results in labeling points in the resulting graph.
 
Option "add plot" results in adding points to an existing graph.
 
Option "y limits" results in setting limits for the y axis. Providing "0,10" results in limits of 0 and 10 for the vertical axis.
 
Option "symbol" sets "pch" to choose symbols as in function \code{\link[graphics]{points}}.
} 

  \item{\bold{OK}}{
Calculate the diversity indices with  \code{\link{diversityresult}} or \code{\link{diversitycomp}}.
} 

  \item{\bold{Plot}}{
Plot the diversity results with the name listed on top (should have been calculated first). This will only provide meaningful results if similar options are provided as when calculating the results.
} 

  \item{\bold{Cancel}}{
Close the window and do not calculate new diversity indices.
} 
}
}

\section{Rank Abundance}{

The window fits and plots rank abundance curves for the community dataset. The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the estimated rank abundance curve after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{subset options}}{
The list shows the variables that can be used for selecting subsets. Option "all" indicates that no subset will be calculated. 
In case a variable is selected, it will be passed as argument for "factor" of functions \code{\link{rankabundance}} or \code{\link{rankabuncomp}}.
} 

  \item{\bold{Subset}}{
Subset chooses which subsets are calculated. 
In case that the value of "." (a period) is selected then function \code{\link{rankabuncomp}} will used to calculate and plot the rank abundance curves (you may need to click in the graph to show where the legend needs to be placed). 
In case another value is chosen, then this will be the argument for "level" of function \code{\link{rankabundance}}.
} 

  \item{\bold{Plot options}}{
The list provides options for scaling the vertical axis. The selection is passed as argument for "scale" of function \code{\link{rankabunplot}}.
 
Option "fit RAD" fits distribution models to the observed rank-abundance distribution with function \code{\link{radfitresult}} and plots the results.
 
Option "add plot" sets addit=T for function \code{\link{rankabunplot}} meaning that the rank abundance curve will be added to an existing graph.
 
Option "x limits"sets xlim for function \code{\link{rankabunplot}}. Providing "1,10" will plot between 1 and 10.
 
Option "y limits"sets ylim for function \code{\link{rankabunplot}}. Providing "2,20" will plot between 2 and 20.
} 

  \item{\bold{OK}}{
Calculate the rank abundance curve with functions \code{\link{rankabundance}} or \code{\link{rankabuncomp}}.
} 

  \item{\bold{Plot}}{
Plot the rank abundance curve with the name listed on top (should have been calculated first) with function \code{\link{rankabunplot}}, or fit models to rank abundance distribution.
} 

  \item{\bold{Cancel}}{
Close the window and do not calculate a new rank abundance curve.
} 
}
}

\section{Renyi diversity profiles}{

The window fits and plots Renyi diversity profiles from the community dataset. The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the diversity profiles after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{Calculation method}}{
The list allows to select the method of calculating the diversity profile. 
 
Options "all" and "separate per site" are passed as argument for "method" of function \code{\link{renyiresult}}. 
 
Option "accumulation" results in using function \code{\link{renyiaccumresult}}.  
  
These options are not valid when \code{\link{renyicomp}} is invoked (see Subset options). 
} 

  \item{\bold{Scale parameters}}{
The "scale parameters" are passed as argument for "scale" for functions \code{\link{renyiresult}}, \code{\link{renyiaccumresult}} or \code{\link{renyicomp}}.
} 

  \item{\bold{Permutations}}{
The "permutations" are passed as argument for "permutations" for functions \code{\link{renyiaccumresult}} or \code{\link{renyicomp}}.
} 

  \item{\bold{subset options}}{
The list shows the variables that can be used for selecting subsets. Option "all" indicates that no subset will be calculated.
   
In case a variable is selected, it will be passed as argument for "factor" of functions \code{\link{renyiresult}} or \code{\link{renyicomp}}.
} 

  \item{\bold{Subset}}{
Subset chooses which subsets are calculated. 
In case that the value of "." (a period) is selected then function \code{\link{renyicomp}} will used to calculate the diversity profile and to plot the curve (you may need to click in the graph to show where the legend needs to be placed). In case another value is chosen, then this will be the argument for "level" of function \code{\link{renyiresult}}.
} 

  \item{\bold{Plot options}}{
Options for plotting passed to function \code{\link{renyiplot}}.
 
Option "evenness profile" sets "evenness=T".
 
Option "evenness profile" sets addit=T meaning that the diversity profiles will be added to an existing graph.
  
Option "y limits"sets ylim. Providing "2,20" will plot between 2 and 20.
 
Option "symbol"sets pch.
 
Option "cex"sets cex.
 
Option "colour" sets col.
} 

  \item{\bold{OK}}{
Calculate the diversity profile with functions \code{\link{renyiresult}}, \code{\link{renyiaccumresult}} or \code{\link{renyicomp}}.
} 

  \item{\bold{Plot}}{
Plot the species accumulation curve with the name listed on top with functions \code{\link{renyiplot}} or \code{\link[vegan:renyi]{persp.renyiaccum}}. The calculation method will determine which plot function is used.
} 

  \item{\bold{Cancel}}{
Close the window and do not calculate a new diversity profile.
} 
}
}

\section{Species abundance as response}{

The window fits and plots regression models for abundance data with a response variable selected from the community dataset and explanatory variables selected from the environmental dataset. (Hint: to analysis species richness patterns, save site-specific species richness (from diversity indices menu) into the environmental data set, and then make the environmental data set to be the community dataset as well). The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the fitted regression model after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you can plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{Model options}}{
Select the method of regression analysis.
 
Option "linear model" fits a simple linear regression model with function \code{\link[stats]{lm}}.
 
Option "Poisson model" fits GLMs with Poisson variance functions and log link functions through function \code{\link[stats]{glm}}.
 
Option "quasi-Poisson model" fits GLMs with quasi-Poisson variance functions and log link functions through function \code{\link[stats]{glm}}.
 
Option "negative binomial model" fits GLMs with negative binomial variance functions and log link functions through function \code{\link[MASS]{glm.nb}}.
 
Option "gam model" fits GAMs with Poisson variance functions and log link functions through function \code{\link[mgcv]{gam}}..
 
Option "gam negbinom model" fits GAMs with negative binomial variance functions and log link functions through function  \code{\link[mgcv]{gam}}.
 
Option "glmmPQL" fits GLMMs with negative binomial variance functions and log link functions through function \code{\link[MASS]{glmmPQL}}.

Option "rpart" fits a regression tree through function \code{\link[rpart]{rpart}}.
 
} 

  \item{\bold{Standardize}}{
Fit the regression to a standardised dataset with function \code{\link[base]{scale}} (only continuous variables are standardised, not categorical variables).
} 

  \item{\bold{Print summary}}{
Provide a summary of the regression with  functions \code{\link[stats]{summary.lm}} , \code{\link[stats]{summary.glm}} or \code{\link[mgcv]{summary.gam}}.
} 

  \item{\bold{Print anova}}{
Provide a summary of the regression with functions \code{\link[stats]{anova.lm}}, \code{\link[stats]{anova.glm}}, \code{\link[mgcv]{anova.gam}}, \code{\link[stats:add1]{drop1}} or \code{\link[car]{Anova}} (latter two type-II ANOVAs only invoced for multiple regression). 
} 

  \item{\bold{add predictions to data frame}}{
Adds the predicted values to the environmental dataset using the model name combined with ".fit" (using the appropriate \code{\link[stats]{predict}} function).
} 

  \item{\bold{Response variable}}{
Type the name of the response variable, or select and double-click from the list that is provided. This variable will be displayed on the left-hand side of the formula (variable ~) and is also the response variable that is plotted in the various result plots. The variable is selected as one of the variables (species) of the community dataset, and is first added to the environmental dataset. When you select the environmental dataset to be the community dataset as well, then you can select variables of the environmental dataset as response variable.
} 

  \item{\bold{Explanatory}}{
Type the right-hand side of the model formula (~ explanatory), or select and double-click for variables and select and click for operators to construct the right-hand side of the model formula.
} 

  \item{\bold{Remove site with name}}{
The name of the site to be removed from the environmental dataset.
} 

  \item{\bold{Plot options}}{
The options provide various functions that can be used to plot regression results of the current model (shown on top of the window; should have been estimated first).
 
Option "diagnostic plots" chooses functions \code{\link[stats]{plot.lm}} or \code{\link[mgcv]{gam.check}} to plot diagnostic plots. For regression trees, the residuals are plotted against the residuals via \code{\link[rpart]{predict.rpart}} and \code{\link[rpart]{residuals.rpart}}. 
 
Option "levene test" chooses function \code{\link[car]{levene.test}} and plots residuals of the selected categorical variable (shown on the right).
 
Option "term plot" chooses functions \code{\link[stats]{termplot}} or \code{\link[mgcv]{plot.gam}} to plot a termplot of the selected categorical variable (shown on the right).
 
Option "effect plot" chooses function \code{\link[effects]{effect}} to plot an effect plot of the selected variable (shown on the right). (The menu option of the R-Commander of models > Graphs plots all the variables).
 
Option "qq plot" chooses function \code{\link[car]{qq.plot}} to plot the residuals from the model.
 
Option "result plot (new)" chooses an appropriate \code{\link[stats]{predict}} function to plot a new plot of the model predictions for the selected variable (shown on the right).
 
Option "result plot (add)" chooses an appropriate \code{\link[stats]{predict}} function to add a new plot of the model predictions for the selected variable (shown on the right)
 
Option "result plot (interpolate)" chooses an appropriate \code{\link[stats]{predict}} function to add a new plot of the model predictions for the selected variable (shown on the right). This model is predicted from a new dataset that only contains 1000 interpolated values for the selected explanatory variable. 
 
Option "cr plot" chooses function \code{\link[car]{cr.plots}} to plot a component + residual plots of the selected variable (shown on the right). (The menu option of the R-Commander of models > Graphs plots all the variables).
 
Option "av plot" chooses function \code{\link[car]{av.plots}} to plot added variable plots of the selected variable (shown on the right). (The menu option of the R-Commander of models > Graphs plots all the variables and has an option of identifying sites with the mouse.)
  
Option "influence plot" chooses function \code{\link[car]{influence.plot}} to plot influence plots. (The menu option of the R-Commander of models > Graphs includes the option of identifying sites with the mouse.) 
 
Option "multcomp" chooses function \code{\link[multcomp]{glht}} to plot simultaneous confidence intervals of the selected categorical variable (shown on the right).
 
Option "rpart" chooses functions \code{\link[rpart]{plot.rpart}} and \code{\link[rpart]{text.rpart}} to plot a dendrogram for the regression tree result.
 
} 

  \item{\bold{Plot variable}}{
Variable of the environmental dataset that is used for some plotting functions.
} 


  \item{\bold{OK}}{
Fit the selected models.
} 

  \item{\bold{Plot}}{
Plot results for the model with name that appears on top. The model options need to apply to the model (e.g. if a GLM method was used to fit the model, this option should also be selected when plotting the results).
} 

  \item{\bold{Cancel}}{
Close the window and do not estimate new regression models.
} 
}
}

\section{Species presence-absence as response}{

The window fits and plots regression models for presence-absence data with a response variable selected from the community dataset and explanatory variables selected from the environmental dataset. The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the fitted regression model after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you can plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{Model options}}{
Select the method of regression analysis.
 
Option "crosstab" calculates a cross-tabulation of the selected response (rescaled as presence-absence) and one selected environmental variable, and estimates a Chi-square test of the contingency table with function \code{\link[stats]{chisq.test}}.
 
Option "binomial model" fits GLMs with binomial variance functions and logit link functions through function \code{\link[stats]{glm}}.
 
Option "quasi-binomial model" fits GLMs with quasi-binomial variance functions and log link functions through function \code{\link[stats]{glm}}.
 
Option "gam model" fits GAMs with binomial variance functions and logit link functions through function \code{\link[mgcv]{gam}}.
 
Option "gam quasi-binomial model" fits GAMs with quasi-binomial variance functions and logit link functions through function \code{\link[mgcv]{gam}}.
 
Option "rpart" fits a regression tree through function \code{\link[rpart]{rpart}}.
 
Option "nnet" fits a forward-feeding artificial neural network through function \code{\link{nnetrandom}}.
 
} 

  \item{\bold{Standardize}}{
Fit the regression to a standardised dataset with function \code{\link[base]{scale}} (only continuous variables are standardised, not categorical variables).
} 

  \item{\bold{Print summary}}{
Provide a summary of the regression with functions \code{\link[stats]{summary.glm}} or \code{\link[mgcv]{summary.gam}}, or use \code{\link[rpart]{summary.rpart}} or \code{\link[nnet:nnet]{summary.nnet}}
} 

  \item{\bold{Print anova}}{
Provide a summary of the regression with  functions \code{\link[stats]{anova.glm}}, \code{\link[mgcv]{anova.gam}}, \code{\link[stats:add1]{drop1}} or \code{\link[car]{Anova}} (latter two type-II ANOVAs only invoced for multiple regression). 
} 

  \item{\bold{add predictions to data frame}}{
Adds the predicted values to the environmental dataset using the model name combined with ".fit" (using the appropriate \code{\link[stats]{predict}} function).
} 

  \item{\bold{Response variable}}{
Type the name of the response variable, or select and double-click from the list that is provided. This variable will be displayed on the left-hand side of the formula (variable >0 ~) and is also the response variable that is plotted in the various result plots. The variable is selected as one of the variables (species) of the community dataset, it will be transformed to presence-absence and is first added to the environmental dataset. When you select the environmental dataset to be the community dataset as well, then you can select variables of the environmental dataset as response variable.
} 

  \item{\bold{Explanatory}}{
Type the right-hand side of the model formula (~ explanatory), or select and double-click for variables and select and click for operators to construct the right-hand side of the model formula.
} 

  \item{\bold{Remove site with name}}{
The name of the site to be removed from the environmental dataset.
} 

  \item{\bold{Plot options}}{
The options provide various functions that can be used to plot regression results of the current model (shown on top of the window; should have been estimated first).
 
Option "tabular" chooses function \code{\link[graphics]{plot}} to plot presence-absence of the response variable against the selected categorical variable (shown on the right).
 
Option "diagnostic plots" chooses functions \code{\link[stats]{plot.lm}} or \code{\link[mgcv]{gam.check}} to plot diagnostic plots. For regression trees and artificial neural networks, the predicted values are plotted against the original presence-absence information.
 
Option "levene test" chooses function \code{\link[car]{levene.test}} and plots residuals of the selected categorical variable (shown on the right).
 
Option "term plot" chooses functions \code{\link[stats]{termplot}} or \code{\link[mgcv]{plot.gam}} to plot a termplot of the selected categorical variable (shown on the right).
 
Option "effect plot" chooses function \code{\link[effects]{effect}} to plot an effect plot of the selected variable (shown on the right). (The menu option of the R-Commander of models > Graphs plots all the variables).
 
Option "qq plot" chooses function \code{\link[car]{qq.plot}} to plot the residuals from the model.
 
Option "result plot (new)" chooses an appropriate \code{\link[stats]{predict}} function to plot a new plot of the model predictions for the selected variable (shown on the right).
 
Option "result plot (add)" chooses an appropriate \code{\link[stats]{predict}} function to add a new plot of the model predictions for the selected variable (shown on the right)
 
Option "result plot (interpolate)" chooses an appropriate \code{\link[stats]{predict}} function to add a new plot of the model predictions for the selected variable (shown on the right). This model is predicted from a new dataset that only contains 1000 interpolated values for the selected explanatory variable. 
 
Option "cr plot" chooses function \code{\link[car]{cr.plots}} to plot a component + residual plots of the selected variable (shown on the right). (The menu option of the R-Commander of models > Graphs plots all the variables.)
 
Option "av plot" chooses function \code{\link[car]{av.plots}} to plot added variable plots of the selected variable (shown on the right). (The menu option of the R-Commander of models > Graphs plots all the variables and has an option of identifying sites with the mouse.)
  
Option "influence plot" chooses function \code{\link[car]{influence.plot}} to plot influence plots. (The menu option of the R-Commander of models > Graphs has an option of identifying sites with the mouse.)
 
Option "multcomp" chooses function \code{\link[multcomp]{glht}} to plot simultaneous confidence intervals of the selected categorical variable (shown on the right).
 
Option "rpart" chooses functions \code{\link[rpart]{plot.rpart}} and \code{\link[rpart]{text.rpart}} to plot a dendrogram for the regression tree result.
 
} 

  \item{\bold{Plot variable}}{
Variable of the environmental dataset that is used for some plotting functions.
} 


  \item{\bold{OK}}{
Fit the selected models.
} 

  \item{\bold{Plot}}{
Plot results for the model with name that appears on top. The model options need to apply to the model (e.g. if a GLM method was used to fit the model, this option should also be selected when plotting the results).
} 

  \item{\bold{Cancel}}{
Close the window and do not estimate new regression models.
} 
}
}

\section{Calculate distance matrix}{

This window calculates a distance matrix from the community dataset and provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new distance matrix that will be calculated after "OK" was clicked.
} 

  \item{\bold{Distance}}{
Ecological distance measure. Passed as argument for "method" for function \code{\link[vegan]{vegdist}}. 
} 

  \item{\bold{Make community dataset)}}{
Make the data frame derived from the new distance matrix the active community data set. This distance matrix can be used directly in the other menus for analysis of ecological distance after selecting the "as.dist" options of these windows.
} 

  \item{\bold{OK}}{
Calculate the distance matrix.
} 

  \item{\bold{Cancel}}{
Close the window and do not calculate a new distance matrix.
}
}
}

\section{Unconstrained ordination}{

The window fits and plots unconstrained ordination models. The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the unconstrained ordination model after "OK" was clicked, or the name of the object that will be plotted when Plot is clicked. In case that you saved a result earlier, then you can plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{Ordination method}}{
Select the method of ordination analysis.
 
Option "PCA" fits a Principal Components Analysis model with function \code{\link[vegan:cca]{rda}}.
 
Option "PCA (prcomp)" fits a Principal Components Analysis model with function \code{\link[stats]{prcomp}}.
 
Option "PCoA" fits a Principal Coordinates Analysis model with function \code{\link[stats]{cmdscale}} using the distance measure selected on the right-hand side (except if the community matrix is interpreted as distance matrix).
 
Option "PCoA (Caillez)" fits a Principal Coordinates Analysis model with function \code{\link[stats]{cmdscale}} using the distance measure selected on the right-hand side (except if the community matrix is interpreted as distance matrix) and setting add=T.
 
Option "CA" fits a Correspondence Analysis (Reciprocal Averaging) model with function \code{\link[vegan]{cca}}.
 
Option "DCA" fits a Detrended Correspondence Analysis model with function \code{\link[vegan]{decorana}}.
 
Option "metaMDS" fits a Non-metric Multidimensional Scaling model with function \code{\link[vegan]{metaMDS}} using the distance measure selected on the right-hand side (except if the community matrix is interpreted as distance matrix).
 
Option "NMS (standard)" fits a Non-metric Multidimensional Scaling model with function \code{\link{NMSrandom}} using the distance measure selected on the right-hand side (except if the community matrix is interpreted as distance matrix).
  
} 

  \item{\bold{Distance}}{
Select the distance measure for the PCoA and NMS methods (other methods have fixed intrinsic distance measures [Euclidean or chi] that can not be changed). 
 
For the methods that provide ordinations based on a distance matrix (PCoA and NMSstandard): passed as argument for "method" for function \code{\link[vegan]{vegdist}} that calculates the distance matrix first.
 
Passed as argument for "distance" for function \code{\link[vegan]{metaMDS}}.
} 

  \item{\bold{PCoA or NMS axes}}{
Select the number of axes to feature in PCoA and NMS results. 
Passed as argument for "k" for functions \code{\link[stats]{cmdscale}}, \code{\link[vegan]{metaMDS}} or \code{\link{NMSrandom}}.
} 

  \item{\bold{NMS permutation}}{
Select the number of permutations for the NMS results. The solution with the lowest stress after all permutations of random starting positions will be provided. 
Passed as argument for "trymax" for function \code{\link[vegan]{metaMDS}} or argument for "perm" for function \code{\link{NMSrandom}}.
} 

  \item{\bold{PCoA or NMS species}}{
Fit species scores to PCoA and NMS results with function \code{\link{add.spec.scores}}. This function adds some other information for PCoA.
} 

  \item{\bold{Model summary}}{
Provide a summary of the ordination with functions \code{\link[vegan:plot.cca]{summary.cca}}, \code{\link[vegan:decorana]{summary.decorana}} orotherwise list the model object.
} 

  \item{\bold{Scaling}}{
Provide the scaling method.
Passed as argument for "scaling" for functions \code{\link[vegan:plot.cca]{summary.cca}}, \code{\link[vegan:decorana]{summary.decorana}} or \code{\link{add.spec.scores}}.
} 

  \item{\bold{as.dist(Community)}}{
Treat the community dataset as a distance matrix.
The community dataset will be used as a distance matrix (via \code{\link[stats:dist]{as.dist}}) for unconstrained ordination methods that use a distance matrix as input (\code{\link[stats]{cmdscale}} and \code{\link{NMSrandom}} for ordination results and via \code{\link[vegan:ordihull]{ordicluster}}, \code{\link[vegan:spantree]{lines.spantree}}, \code{\link{ordicluster2}}, \code{\link{ordinearest}} or \code{}\link{distdisplayed} for plotting options). 
} 

  \item{\bold{Plot method}}{
The options provide various functions that can be used to plot ordination results, or to add information to ordination diagrams.
 
Option "plot" chooses function \code{\link[vegan]{plot.cca}} to plot results from \code{\link[vegan:cca]{rda}}, \code{\link[vegan]{cca}} , \code{\link[vegan]{metaMDS}} or \code{\link[vegan]{decorana}} and function \code{\link[graphics]{plot}} to plot the other ordination results (obtained by function \code{\link[vegan]{scores}}).
 
Option "ordiplot" chooses function \code{\link[vegan]{ordiplot}} to plot ordination results.
 
Option "ordiplot empty" chooses function \code{\link[vegan]{ordiplot}} to plot ordination results, but sites and species will be invisible.
 
Option "identify sites" chooses function \code{\link[vegan:ordiplot]{identify.ordiplot}} to add names of sites to site symbols (circles) created by function \code{\link[vegan]{ordiplot}}. You can choose where the name is added by left-clicking in the quadrant next to the symbol where you want to symbol to be plotted. You can stop identifying sites by right-clicking.
 
Option "identify species" chooses function \code{\link[vegan:ordiplot]{identify.ordiplot}} to add names of species to species symbols (crosses) created by function \code{\link[vegan]{ordiplot}}. You can choose where the name is added by left-clicking in the quadrant next to the symbol where you want to symbol to be plotted. You can stop identifying species by right-clicking.
 
Option "text sites" chooses function \code{\link[vegan:ordiplot]{text.ordiplot}} to add names of all sites to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "text species" chooses function \code{\link[vegan:ordiplot]{text.ordiplot}} to add names of all species to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "points sites" chooses function \code{\link[vegan:ordiplot]{points.ordiplot}} to add symbols for all sites to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "points species" chooses function \code{\link[vegan:ordiplot]{points.ordiplot}} to add symbols for all species to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "origin axes" adds a horizontal and vertical line through the origin of the ordination graph (the origin is the location with coordinates [0,0]).
 
Option "envfit" chooses function \code{\link[vegan]{envfit}} to add information for the variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. 
 
Option "ordihull" chooses function \code{\link[vegan]{ordihull}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordiarrows" chooses function \code{\link[vegan:ordihull]{ordiarrows}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordisegments" chooses function \code{\link[vegan:ordihull]{ordisegments}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordispider" chooses function \code{\link[vegan:ordihull]{ordispider}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordiellipse" chooses function \code{\link[vegan:ordihull]{ordiellipse}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordisurf" chooses function \code{\link[vegan:ordihull]{ordisurf}} to add information for the continuous variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordicluster" chooses function \code{\link[vegan:ordihull]{ordicluster}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix.) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordispantree" chooses function \code{\link[vegan:spantree]{lines.spantree}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordibubble" chooses function \code{\link{ordibubble}} to add information for the continuous variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordisymbol" chooses function \code{\link{ordisymbol}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. Make sure that you click in the graph to show where the legend should be placed!
 
Option "ordivector" chooses function \code{\link{ordivector}} to add information on the selected species of the community dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. You should first make the community dataset the environmental datset to get the list of species on the right-hand side.
 
Option "ordivector interpretation" chooses function \code{\link{ordivector}} to add information on the selected species of the community dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. You should first make the community dataset the environmental datset to get the list of specie son the right-hand side. The function will drop down perpendicular lines from each site to the line connecting the origin and the species position.
 
Option "ordicluster2" chooses function \code{\link{ordicluster2}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordinearest" chooses function \code{\link{ordinearest}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordiequilibriumcircle" chooses function \code{\link{ordiequilibriumcircle}} to plot an equilibrium circle to ordination diagrams created by function \code{\link[vegan]{ordiplot}} from the Principal Components Analysis fitted by \code{\link[vegan:cca]{rda}}.
 
Option "distance displayed" compares the distances between each pair of sites in a distance matrix (with distance measure selected in window above) with distances in ordination diagrams created by function \code{\link[vegan]{ordiplot}} by means of function \code{\link{distdisplayed}}.
 
Option "screeplot.cca" provides a screeplot for PCA results obtained by function \code{\link[vegan:cca]{rda}} by means of function \code{\link[vegan]{screeplot.cca}}.
 
Option "stress" provides a stress plot (Shepard diagram) for NMS results obtained by function \code{\link[vegan]{metaMDS}} by means of function \code{\link[vegan:goodness.metaMDS]{stressplot}}.
 
Option "coenocline" fits coenoclines for all species to the first ordination axis of ordination diagrams created by function \code{\link[vegan]{ordiplot}} by means of function \code{\link{ordicoeno}}.
 } 

  \item{\bold{Plot variable}}{
Variable of the environmental dataset that is used for some plotting functions. For Plot method "ordivector", make the community dataset the environmental dataset first. Some other plot methods may also work with the community dataset as the environmental dataset as well (e.g. "ordibubble", "ordisurf"). Some methods run into problems when the variable has missing observations: in this case, you may need to repeat the ordination analysis after removing sites with missing observations for the variable with the "remove NA" option of the Community dataset menu list.
} 

  \item{\bold{axes}}{
The position of the axes of the ordination result to be plotted in the ordination diagram ("1,2" selects the first two axes of the ordination result).
Passed as argument for "choices" for functions \code{\link[vegan]{plot.cca}}, \code{\link[vegan]{scores}} or \code{\link[vegan]{ordiplot}}.
} 

  \item{\bold{add scores to dataframe}}{
Adds the scores of the sites from the \code{\link[vegan]{ordiplot}} graph to the environmental dataset using the model name combined with ".ax1" and ".ax2".
} 

  \item{\bold{cex}}{
The size of the characters in the resulting plot when "Plot" is clicked.
} 

  \item{\bold{colour}}{
The colour of the resulting plot when "Plot" is clicked.
} 

  \item{\bold{OK}}{
Fit the selected models.
} 

  \item{\bold{Plot}}{
Plot results for the model with name that appears on top. The model options need to apply to the model (e.g. if \code{\link[vegan:cca]{rda}} was used to fit the model, this option should also be selected when plotting the results).
} 

  \item{\bold{Cancel}}{
Close the window and do not fit or plot ordination models.
} 
}
}

\section{Constrained ordination}{

The window fits and plots constrained ordination models to the community dataset, using variables of the environmental dataset to contrain the ordination model (direct gradient analysis, canonical ordination analysis). The menu provides the following options:

\itemize{
  \item{\bold{Save result as}}{
The name for the new object that will save the results from the unconstrained ordination model after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you can plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{Ordination method}}{
Select the method of ordination analysis.
 
Option "RDA" fits a Redundancy Analysis model with function \code{\link[vegan:cca]{rda}}.
  
Option "CCA" fits a Canonical Correspondence Analysis (Reciprocal Averaging) model with function \code{\link[vegan]{cca}}.
 
Option "capscale" fits a scaled Constrained Analysis of Principal Coordinates (distance-based Redundancy Analysis) with function \code{\link[vegan]{capscale}} using the distance measure selected on the right-hand side (except if the community matrix is interpreted as distance matrix).
 
Option "CAPdiscrim" fits a Constrained Analysis of Principal Coordinates (based on discriminant analysis) with function \code{\link{CAPdiscrim}} using the distance measure selected on the right-hand side (except if the community matrix is interpreted as distance matrix) and the categorical variable selected as explanatory variable.
  
Option "prc" fits principal response curves with function \code{\link[vegan]{prc}}. To implement the example provided in the documentation for the \code{\link[vegan]{prc}} function, you need to include the additional steps of defining \code{pyrifos.env <- data.frame(dose,week)} and making this data set the environmental data set.
 
Option "multiconstrained (RDA)" provides the first row of all ANOVA results (\code{\link[vegan]{anova.cca}}) for all possible pairwise combinations of the levels of the first explanatory variable (assumed to be a categorical variable) through function \code{\link{multiconstrained}} with method="rda". (When you change contrast to a particular contrast indicator, you obtain an ordination result that can be analyzed further. For several plotting options, you need to change the community and environmental datasets to "newcommunity" and "newenvdata").
 
Option "multiconstrained (CCA)" provides the first row of all ANOVA results (\code{\link[vegan]{anova.cca}}) for all possible pairwise combinations of the levels of the first explanatory variable (assumed to be a categorical variable) through function \code{\link{multiconstrained}} with method="cca". (When you change contrast to a particular contrast indicator, you obtain an ordination result that can be analyzed further.  For several plotting options, you need to change the community and environmental datasets to "newcommunity" and "newenvdata").
 
Option "multiconstrained (capscale)" provides the first row of all ANOVA results (\code{\link[vegan]{anova.cca}}) for all possible pairwise combinations of the levels of the first explanatory variable (assumed to be a categorical variable) through function \code{\link{multiconstrained}} with method="capscale". (When you change contrast to a particular contrast indicator, you obtain an ordination result that can be analyzed further.  For several plotting options, you need to change the community and environmental datasets to "newcommunity" and "newenvdata").
 
} 

  \item{\bold{Distance}}{
Select the distance measure for the CAP methods (other methods have fixed intrinsic distance measures [Euclidean or chi] that can not be changed). 
Passed as argument for "dist" for function \code{\link[vegan]{capscale}} or \code{\link{CAPdiscrim}}. This argument is ignored by the actual functions if the community dataset is interpreted to be a distance matrix already.
} 

  \item{\bold{Model summary}}{
Provide a summary of the ordination with functions \code{\link[vegan:plot.cca]{summary.cca}} or \code{\link[vegan:prc]{summary.prc}},or otherwise list the model object (\code{\link{CAPdiscrim}}).
} 

  \item{\bold{as.dist(Community)}}{
Treat the community dataset as a distance matrix.
The community dataset will be used as a distance matrix (via \code{\link[stats:dist]{as.dist}}) for constrained ordination methods that can use a distance matrix as input (\code{\link[vegan]{capscale}} or \code{\link{CAPdiscrim}} for ordination results and via \code{\link[vegan:ordihull]{ordicluster}}, \code{\link[vegan:spantree]{lines.spantree}}, \code{\link{ordicluster2}}, \code{\link{ordinearest}} or \code{\link{distdisplayed}} for plotting options).
} 

  \item{\bold{Scaling}}{
Provide the scaling method. This option is not available for function \code{\link{CAPdiscrim}}.
Passed as argument for "scaling" for function \code{\link[vegan:plot.cca]{summary.cca}} or \code{\link[vegan:prc]{summary.prc}}.
} 

  \item{\bold{permutations}}{
Select the number of permutations for testing the significance of the constrained ordination by Monte-Carlo randomization tests. The default of "0" means that no permutation test will be done.
Passed as argument for "permutations" for functions \code{\link[vegan:anova.cca]{permutest.cca}}, \code{\link{CAPdiscrim}} or \code{\link[vegan]{envfit}} (one of the plotting options) or as argument for perm.max for function \code{\link[vegan]{anova.cca}} (which is also called by \code{\link{multiconstrained}}).
} 

  \item{\bold{Explanatory}}{
Type the right-hand side of the model formula (~ explanatory), or select and double-click for variables and select and click for operators to construct the right-hand side of the model formula. It is possible to include conditional variables for partial ordination analysis, except for function \code{\link{CAPdiscrim}} and \code{\link[vegan]{prc}}. For function \code{\link[vegan]{prc}}, the explanatory variables should be separated by a comma and indicate the "treatment" and "time" factors.
} 

  \item{\bold{Plot method}}{
The options provide various functions that can be used to plot ordination results, or to add information to ordination diagrams.
 
Option "plot" chooses function \code{\link[vegan]{plot.cca}} to plot results from \code{\link[vegan:cca]{rda}}, \code{\link[vegan]{cca}} or \code{\link[vegan]{capscale}}, function \code{\link[vegan:prc]{plot.prc}} to plot results from \code{\link[vegan]{prc}} and function \code{\link[graphics]{plot}} to plot the other ordination results (obtained by function \code{\link[vegan]{scores}}).
 
Option "ordiplot" chooses function \code{\link[vegan]{ordiplot}} to plot ordination results.
 
Option "ordiplot empty" chooses function \code{\link[vegan]{ordiplot}} to plot ordination results, but sites and species will be invisible.
 
Option "identify sites" chooses function \code{\link[vegan:ordiplot]{identify.ordiplot}} to add names of sites to site symbols (circles) created by function \code{\link[vegan]{ordiplot}}. You can choose where the name is added by left-clicking in the quadrant next to the symbol where you want to symbol to be plotted. You can stop identifying sites by right-clicking.
 
Option "identify species" chooses function \code{\link[vegan:ordiplot]{identify.ordiplot}} to add names of species to species symbols (crosses) created by function \code{\link[vegan]{ordiplot}}. You can choose where the name is added by left-clicking in the quadrant next to the symbol where you want to symbol to be plotted. You can stop identifying species by right-clicking.
 
Option "identify centroids" chooses function \code{\link[vegan:ordiplot]{identify.ordiplot}} to add names of centroids to centroid symbols (X) created by function \code{\link[vegan]{ordiplot}}. You can choose where the name is added by left-clicking in the quadrant next to the symbol where you want to symbol to be plotted. You can stop identifying species by right-clicking.
 
Option "text sites" chooses function \code{\link[vegan:ordiplot]{text.ordiplot}} to add names of all sites to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "text species" chooses function \code{\link[vegan:ordiplot]{text.ordiplot}} to add names of all species to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "text centroids" chooses function \code{\link[vegan:ordiplot]{text.ordiplot}} to add names of all centroids to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
  
Option "points sites" chooses function \code{\link[vegan:ordiplot]{points.ordiplot}} to add symbols for all sites to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "points species" chooses function \code{\link[vegan:ordiplot]{points.ordiplot}} to add symbols for all species to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "points centroids" chooses function \code{\link[vegan:ordiplot]{points.ordiplot}} to add symbols for all centroids to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "origin axes" adds a horizontal and vertical line through the origin of the ordination graph (the origin is the location with coordinates [0,0]).
 
Option "envfit" chooses function \code{\link[vegan]{envfit}} to add information for the variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. 
 
Option "ordihull" chooses function \code{\link[vegan]{ordihull}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordihull.centroids" chooses function \code{\link{ordihull.centroids}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordiarrows" chooses function \code{\link[vegan:ordihull]{ordiarrows}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordisegments" chooses function \code{\link[vegan:ordihull]{ordisegments}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordispider" chooses function \code{\link[vegan:ordihull]{ordispider}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordispider.centroids" chooses function \code{\link{ordispider.centroids}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordiellipse" chooses function \code{\link[vegan:ordihull]{ordiellipse}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordisurf" chooses function \code{\link[vegan:ordihull]{ordisurf}} to add information for the continuous variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordicluster" chooses function \code{\link[vegan:ordihull]{ordicluster}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordispantree" chooses function \code{\link[vegan:spantree]{lines.spantree}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordibubble" chooses function \code{\link{ordibubble}} to add information for the continuous variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordisymbol" chooses function \code{\link{ordisymbol}} to add information for the categorical variable of the environmental dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. Make sure that you click in the graph to show where the legend should be placed!
 
Option "ordivector" chooses function \code{\link{ordivector}} to add information on the selected species of the community dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. You should first make the community dataset the environmental datset to get the list of species on the right-hand side.
 
Option "ordivector interpretation" chooses function \code{\link{ordivector}} to add information on the selected species of the community dataset selected on the right-hand side to ordination diagrams created by function \code{\link[vegan]{ordiplot}}. You should first make the community dataset the environmental datset to get the list of specie son the right-hand side. The function will drop down perpendicular lines from each site to the line connecting the origin and the species position.
 
Option "ordicluster2" chooses function \code{\link{ordicluster2}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "ordinearest" chooses function \code{\link{ordinearest}} to add information (with distance measure selected in window above - except if the community matrix is interpreted as distance matrix) to ordination diagrams created by function \code{\link[vegan]{ordiplot}}.
 
Option "distance displayed" compares the distances between each pair of sites in a distance matrix (with distance measure selected in window above) with distances in ordination diagrams created by function \code{\link[vegan]{ordiplot}} by means of function \code{\link{distdisplayed}}.
 
Option "coenocline" fits coenoclines for all species to the first ordination axis of ordination diagrams created by function \code{\link[vegan]{ordiplot}} by means of function \code{\link{ordicoeno}}.
 
} 

  \item{\bold{Plot variable}}{
Variable of the environmental dataset that is used for some plotting functions. For Plot method "ordivector", make the community dataset the environmental dataset first. Some other plot methods may also work with the community dataset as the environmental dataset as well (e.g. "ordibubble", "ordisurf"). Some methods run into problems when the variable has missing observations: in this case, you may need to repeat the ordination analysis after removing sites with missing observations for the variable with the "remove NA" option of the Community dataset menu list.
} 

  \item{\bold{axes}}{
The position of the axes of the ordination result to be plotted in the ordination diagram ("1,2" selects the first two axes of the ordination result).
Passed as argument for "choices" for functions \code{\link[vegan]{plot.cca}}, \code{\link[vegan]{scores}} or \code{\link[vegan]{ordiplot}}.
} 

  \item{\bold{add scores to dataframe}}{
Adds the scores of the sites from the \code{\link[vegan]{ordiplot}} graph to the environmental dataset using the model name combined with ".ax1" and ".ax2".
} 

  \item{\bold{cex}}{
The size of the characters in the resulting plot when "Plot" is clicked.
} 

  \item{\bold{colour}}{
The colour of the resulting plot when "Plot" is clicked.
} 

  \item{\bold{OK}}{
Fit the selected models.
} 

  \item{\bold{Plot}}{
Plot results for the model with name that appears on top. The model options need to apply to the model (e.g. if \code{\link[vegan:cca]{rda}} was used to fit the model, this option should also be selected when plotting the results).
} 

  \item{\bold{Cancel}}{
Close the window and do not fit or plot ordination models.
} 
}
}

\section{Clustering}{

This window performs various methods of cluster analysis based on the information of the community dataset. The menu provides the following options:

\itemize{
  \item{\bold{Save cluster as}}{
The name for the new object that will save the results from the cluster analysis after "OK" was clicked, or the name of the object that will be plotted when "Plot" is clicked. In case that you saved a result earlier, then you can plot the result by typing in the name of previous result first in this box.
} 

  \item{\bold{Cluster method}}{
Select the method of ordination analysis.
   
Option "hclust" results in a cluster analysis fitted by function \code{\link[stats]{hclust}}. The distance for the distance matrix derived from the community dataset is selected on the right-hand side.
 
Option "agnes" results in a cluster analysis fitted by function \code{\link[cluster]{agnes}}. The distance for the distance matrix derived from the community dataset is selected on the right-hand side.
 
Option "diana" results in a cluster analysis fitted by function \code{\link[cluster]{diana}}. The distance for the distance matrix derived from the community dataset is selected on the right-hand side.
 
Option "kmeans" results in a cluster analysis fitted by function \code{\link[stats]{kmeans}}. This method is based on the Euclidean distance.
 
Option "cascadeKM" results in a cluster analysis fitted by function \code{\link[vegan]{cascadeKM}}. This method is based on the Euclidean distance as it is based on K-means clustering.
 
Option "pam" results in a cluster analysis fitted by function \code{\link[cluster]{pam}}. The distance for the distance matrix derived from the community dataset is selected on the right-hand side.
 
Option "clara" results in a cluster analysis fitted by function \code{\link[cluster]{clara}}. The distance for the distance matrix derived from the community dataset is selected on the right-hand side.
 
Option "fanny" results in a cluster analysis fitted by function \code{\link[cluster]{fanny}}. The distance for the distance matrix derived from the community dataset is selected on the right-hand side.
  
} 

  \item{\bold{Distance}}{
Ecological distance measure used for the distance matrix.
Passed as argument for "method" for function \code{\link[vegan]{vegdist}}.
} 

  \item{\bold{as.dist(Community)}}{
Treat the community dataset as a distance matrix.
The community dataset will be used as a distance matrix (via \code{\link[stats:dist]{as.dist}}). This option is not available for \code{\link[stats]{kmeans}}).
} 

  \item{\bold{cluster summary}}{
Provide the results of the cluster analysis with \code{\link[cluster]{summary.agnes}}, \code{\link[cluster]{summary.diana}}, \code{\link[cluster]{summary.pam}}, \code{\link[cluster]{summary.clara}} or \code{\link[cluster:print.fanny]{summary.fanny}} or provide results of \code{\link[stats]{hclust}} or \code{\link[stats]{kmeans}}
} 

  \item{\bold{cophenetic correlation}}{
Calculate the correlation of the distances in the distance matrix with the distances in the dendrogram (estimated with function \code{\link[stats]{cophenetic}}) by the Mantel test (\code{\link[vegan]{mantel}}). It only works for hierarchical clustering methods (\code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}} and \code{\link[cluster]{diana}}).
} 

  \item{\bold{clusters}}{
Determine a fixed number of clusters.
  
This number selects the number of clusters to be calculated by the non-hierarchical cluster methods as it is passed as argument for "centers" for function \code{\link[stats]{kmeans}} and argument for "k" for functions \code{\link[cluster]{pam}}, \code{\link[cluster]{clara}} and \code{\link[cluster]{fanny}}.
 
This number selects the number of groups for the partition with the largest number of groups of the cascade as it is passed as argument for "sup.gr" for function \code{\link[vegan]{cascadeKM}} (the argument for "inf.gr" is set to "2").
 
This number selects the number of clusters to be reported for cluster membership for hierarchical clustering methods (\code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}} and \code{\link[cluster]{diana}}) as determined by function \code{\link[stats]{cutree}}: passed as argument for "k" for this function.
 
This number selects the number of rectangles to be plotted on a dendrogram with plotting option of "rectangles": passed as argument for "k" for function \code{\link[stats]{rect.hclust}}.
 
This number selects the number of clusters to be plotted with plotting option of "pruned dendrogram": passed as argument for "k" for function \code{\link[maptree]{clip.clust}}.

} 

  \item{\bold{Save cluster membership}}{
Save the identity of the cluster to which each site belongs into the environmental data set. For hierarchical clustering methods (\code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}} and \code{\link[cluster]{diana}}) as determined by function \code{\link[stats]{cutree}}, with parameter "k" obtained from the box above.
} 

  \item{\bold{Cluster options}}{
Choose the options that are available for some of the hierarchical clustering methods.
   
Options "average", "single", "complete", "ward", "median" and "centroid" can be passed meaningfully as argument for "method" for \code{\link[stats]{hclust}}.
 
Options "average", "single", "complete", "ward" and "weighted" can be passed meaningfully as argument for "method" for \code{\link[cluster]{agnes}}.
 
} 

  \item{\bold{Plot options}}{
Choose the options that are available for plotting hierarchical clustering results (except for "cascadeKM").
 
Option "dendrogram1" selects function \code{\link[stats:hclust]{plot.hclust}}, \code{\link[cluster]{plot.agnes}} or \code{\link[cluster]{plot.diana}} to plot clustering results.
 
Option "dendrogram2" selects function \code{\link[stats:hclust]{plot.hclust}}, \code{\link[cluster]{plot.agnes}} or \code{\link[cluster]{plot.diana}} to plot clustering results with argument hang set to -1. This option will result in each branch of the dendrogram to reach "ground level".
 
Option "rectangles" selects function \code{\link[stats]{rect.hclust}} to plot rectangles around the number of cluster determined by option "clusters" selected above.
 
Option "pruned dendrogram" selects function \code{\link[maptree]{clip.clust}} to prune the cluster to the number of cluster selected by option "clusters" selected above. This option may only work with cluster results obtained by \code{\link[stats:hclust]{plot.hclust}}. 
 
Option "kgs" selects function \code{\link[maptree]{kgs}} as one method of selecting the optimal number of clusters and plots its results.
 
Option "cophenetic" uses function \code{\link[stats]{cophenetic}} to the distance in the dendrogram against the distance of the distance matrix (calculated earlier for the clustering algorithm). A reference line (y=x) is added to the graph.
 
Option "cascadeKM" selects function \code{\link[vegan:cascadeKM]{plot.cascadeKM}} to plot resuls obtained by function \code{\link[vegan]{cascadeKM}}. 
} 

  \item{\bold{cex}}{
The size of the characters in the resulting plot when "Plot" is clicked.
} 

  \item{\bold{colour}}{
The colour of the resulting plot when "Plot" is clicked.
} 

  \item{\bold{OK}}{
Fit the selected models.
} 

  \item{\bold{Plot}}{
Plot results for the cluster with name that appears on top. Plotting will only be meaningfull for hierarchical methods (\code{\link[stats]{hclust}}, \code{\link[cluster]{agnes}} and \code{\link[cluster]{diana}}).
} 

  \item{\bold{Cancel}}{
Close the window and do not analyse or plot clusters..
} 
}
}

\section{Compare distance matrices}{

This window calculates a distance matrix from the community dataset. This distance matrix can be analysed by a Mantel, MRPP or ANOSIM test based on information from the environmental dataset. You can compare two different community datasets if you make one the community dataset and the other one the environmental dataset. The menu provides the following options:

\itemize{
  \item{\bold{Type of test}}{
Selects the type of test to be used.
 
Option "mantel" results in a Mantel test estimated by function \code{\link[vegan]{mantel}}. The distance for distance matrix derived from the community dataset is selected below, the distance to be derived from the environmental dataset is selected on the right-hand side.
 
Option "anosim" results in a ANOSIM test estimated by function \code{\link[vegan]{anosim}} as summarized by \code{\link[vegan:anosim]{summary.anosim}}. The distance measure for the distance matrix derived from the community dataset is selected below, the categorical variable of the environmental dataset is selected at the right-hand side.
 
Option "mrpp" results in a MRPP test estimated by function \code{\link[vegan]{mrpp}}. The distance measure for the distance matrix derived from the community dataset is selected below, the categorical variable of the environmental dataset is selected at the right-hand side.
 
Option "rankindex" results in a series of Mantel tests with a series of distance measures selected by function \code{\link[vegan]{rankindex}} for the community dataset and the Euclidean distance for the environmental dataset (except for datasets that contain factors where \code{\link[cluster]{daisy}} is used).
 
} 

  \item{\bold{Environmental distance}}{
The environmental distance is only used for the test option of "mantel" (test option of "rank index" makes its own choice in between "daisy" or "euclidean" distance). The distance determines the type of distance matrix that is obtained from the environmental data set.
 
Option "daisy" results in function \code{\link[cluster]{daisy}} to be used for providing the distance matrix. This is the only realistic method for environmental datasets that contain categorical variables.
 
The other options are passed as arguments for "method" for function \code{\link[vegan]{vegdist}}.
   
} 

  \item{\bold{Community distance}}{
Ecological distance measure used for the distance matrix obtained from the community data set.
 
Passed as argument for "method" for function \code{\link[vegan]{vegdist}}.
 
For the "rankindex" type of test, a series of distance measures are tested automatically.
 
} 

  \item{\bold{Environmental variable}}{
Selection of the environmental variable(s). Some methods run into problems when the variable has missing observations: in this case, you may need to repeat the ordination analysis after removing sites with missing observations for the variable with the "remove NA" option of the Community dataset menu list.
 
For test option "mantel", when "all" is selected, then the distance matrix is calculated for all variables of the environmental dataset. For environmental datasets with some categorical variables, only environmental distance "daisy" will result in actual distance matrices.
 
For test option "mantel", when a variable is selected, then the distance matrix is only calculated for that variable. In case that the variable is categorical, then the \code{\link[cluster]{daisy}} distance is used automatically.
 
For test option "anosim", the selected environmental variable is passed as argument for "grouping" for function \code{\link[vegan]{anosim}}.
 
For test option "mrpp", the selected environmental variable is passed as argument for "grouping" for function \code{\link[vegan]{mrpp}}.
 
For test option "rankindex", when "all" is selected, then the environmental dataset is passed as argument for "grad" for function \code{\link[vegan]{rankindex}}.
 
For test option "rankindex", the selected variable is passed as argument for "grad" for function \code{\link[vegan]{rankindex}}. 
 
} 

  \item{\bold{as.dist(Community)}}{
Treat the community dataset as a distance matrix.
The community dataset will be used as a distance matrix (via \code{\link[stats:dist]{as.dist}}).
} 


  \item{\bold{Plot results}}{
Plots the distances of the community dataset against the distance of the environmental dataset for test options "mantel", "anosim" and "mrpp". For categorical variables (the only possibility for "anosim" and "mrpp"), environmental distance equals "0" if sites belong to the same group and "1" if they belong to a different group except if they are ordered categorical variables (depending on the results of the \code{\link[cluster]{daisy}} distance; for ordered factors, it is recommended to create a new factor that is unordered and use this variable for the analysis; see \code{\link[base]{factor}}).
} 

  \item{\bold{permutations}}{
Number of permutations. 
Passed as argument for "permutations" for functions \code{\link[vegan]{mantel}}, \code{\link[vegan]{anosim}} and \code{\link[vegan]{mrpp}}.

} 

  \item{\bold{correlation}}{
Correlation method. 
Passed as argument for "method" for function \code{\link[vegan]{mantel}}.
} 

  \item{\bold{OK}}{
Estimate the selected tests.
} 

  \item{\bold{Cancel}}{
Close the window and do estimate a new test.
} 
}
}


\value{
None
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{http://www.worldagroforestry.org/treesandmarkets/tree_diversity_analysis.asp} 
}   

\author{Roeland Kindt (with some help from Jari Oksanen) }

\keyword{multivariate} 



