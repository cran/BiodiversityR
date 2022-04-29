\name{ensemble.concave.hull}
\alias{ensemble.concave.hull}
\alias{ensemble.concave.venn}
\alias{ensemble.concave.union}

\title{
Analysis of Niche Overlap in Environmental Space for Changed Climates via Concave Hulls
}

\description{
Building on methodologies described by Pironon et al. (\doi{10.1038/s41558-019-0585-7}), function \code{\link[BiodiversityR]{ensemble.concave.hull}} constructs two hulls in environmental space for the baseline and a changed (typically a future climate, but possibly also a historical or paleo-climate) for a focal species. Functions \code{\link[BiodiversityR]{ensemble.concave.venn}} and \code{\link[BiodiversityR]{ensemble.concave.union}} create a third hull for candidate accessions that represent different geographies and/or different species. Subsequently overlaps between hulls are investigated. Information is also provided for each accession of the focal species in the novel climate if these are included within the hull of the candidate accessions. 
}

\usage{
ensemble.concave.hull(
        baseline.data,
        change.data,
        complete.cases = TRUE,
        VIF = TRUE, VIF.max = 20, VIF.silent = TRUE,
        method = c("rda", "pca", "prcomp"),
        ax1 = 1, ax2 = 2,
        concavity = 2.5,
        buffer.dist = NA,
        ggplot = TRUE)

ensemble.concave.venn(
    x,
    candidate.data,
    concavity = x$concavity,
    buffer.dist = x$buffer.dist,
    ggplot = TRUE,
    show.candidate.points = TRUE)

ensemble.concave.union(
    x,
    candidate.venns,
    buffer.dist = x$buffer.dist,
    ggplot = TRUE,
    show.candidate.points = TRUE)
}

\arguments{
  \item{baseline.data}{ data.frame with climatic variables for the accessions in the baseline climate.  }
  \item{change.data}{ data.frame with climatic variables for the accessions in the changed (potentially future) climate.  }
  \item{complete.cases}{ Reduce cases with those without missing data via \code{\link[stats]{complete.cases}}.  }
  \item{VIF}{ Select a subset of climatic variables via \code{\link[BiodiversityR]{ensemble.VIF.dataframe}}.  }
  \item{VIF.max}{ Argument setting for \code{\link[BiodiversityR]{ensemble.VIF.dataframe}}.  }
  \item{VIF.silent}{ Argument setting for \code{\link[BiodiversityR]{ensemble.VIF.dataframe}}.  }
  \item{method}{ Method of constructing the hulls; see details.  }  
  \item{ax1}{ Idex for the first ordination axis to be analyzed; see also \code{\link[vegan]{scores}}.  }
  \item{ax2}{ Index for second ordination axis to be analyzed; see also \code{\link[vegan]{scores}}.  }
  \item{concavity}{ A relative measure of concavity used by \code{\link[concaveman]{concaveman}}.  }
  \item{buffer.dist}{ Buffer width used internally by \code{\link[sf]{st_buffer}}.  }
  \item{ggplot}{ Should a ggplot object be included in the output?  }

  \item{x}{ Output similar to those of \code{\link[BiodiversityR]{ensemble.concave.hull}}.  }
  \item{candidate.data}{ data.frame with climatic variables for candidate accessions such as accessions from other geographical areas or other species.  }
  \item{show.candidate.points}{ Should the ggplot object show the locations of the candidate accessions?  }

  \item{candidate.venns}{ list with outputs from the \code{\link[BiodiversityR]{ensemble.concave.venn}} function.  }

}

\details{

Whereas the metholology of Pironon et al. (2019) uses convex hulls, concave hulls can also be used in the methodology provided here. Convex hulls will be obtained by using large values for the \code{concavity} argument (see the description for the \code{\link[concaveman]{concaveman}} function). By using more concave hulls, the influence of outliers on measures of niche overlap can be reduced.

Three methods are available for mapping accessions in environmental space. Methods \code{pca} and \code{prcomp} use principal components analysis, respectively via the \code{\link[vegan]{rda}} and \code{\link[stats]{prcomp}} functions. In both the methods, climatic variables are scaled. As results with \code{pca} are also rescaled via \code{\link[BiodiversityR]{caprescale}}, both methods of \code{pca} and \code{prcomp} should theoretically result in the same configurations.

Method \code{rda} internally uses \code{\link[vegan]{envfit}} to select a subset of climatic variables that are significantly correlated (P <= 0.05, R2 >= 0.50) with the first two axes of a redundancy analysis that uses the climate (baseline vs. changed) as predictor variable.

Candidate accessions are mapped in the environmental space created by \code{\link[BiodiversityR]{ensemble.concave.hull}} via prediction methods available from \code{\link[vegan]{predict.cca}} and \code{\link[stats]{predict.prcomp}}.

Function \code{\link[BiodiversityR]{ensemble.concave.union}} combines candidate hulls obtained from \code{\link[BiodiversityR]{ensemble.concave.venn}}, using \code{\link[sf]{st_union}} internally.

Both \code{\link[BiodiversityR]{ensemble.concave.venn}} and \code{\link[BiodiversityR]{ensemble.concave.union}} return measures of niche overlap based on areas of overlap between the candidate hull and the part of hull for the changed climate that is not covered by the hull for the baseline climate. These functions also indicate for each of the accessions of the focal species in the changed climate whether they occur in a novel climate (novel == TRUE; this information was obtained by \code{\link[BiodiversityR]{ensemble.concave.hull}}) and whether they are inside the hull of the candidate accessions (candidate.in == TRUE).

The optional plot shows the locations of accessions for the changed climate. For \code{\link[BiodiversityR]{ensemble.concave.hull}}, colouring is based on having novel climates (not occurring in the overlap between the two hulls) or not. For the other functions, locations are only shown for accessions with novel climates. Colouring is based on being inside the hull for the candidate accessions or not.

}

\value{
Function \code{ensemble.concave.hull} returns a list with following elements:

- rda.object: result of the ordination method used; 
- method: method used in the function; 
- baseline.hull: polygon for the hull for the baseline climate; 
- baseline.area: area of the baseline hull; 
- change.hull: polygon for the hull for the changed climate; 
- change.area: area of the hull for the changed climate; 
- overlap.hull: polygon for the overlap (intersection) of the baseline and changed hull; 
- overlap.area: area of the overlap hull;  
- novel.hull: polygon for the part of the changed hull that does not cover the baseline hull;  
- change.area: area of the novel hull; 
- buffer.dist: distance used in checking whether accessions are in novel conditions;
- change.points: plotting coordinates and details on novel conditions for accessions of the changed climate;
- baseline.points: plotting coordinates for accessions of the baseline climate
 
}

\author{ Roeland Kindt (World Agroforestry Centre) and Maarten van Zonneveld (World Vegetable Center)}

\references{
Pironon et al. (2019). Potential adaptive strategies for 29 sub-Saharan crops under future climate change. Nat. Clim. Chang. 9: 758-736. \doi{10.1038/s41558-019-0585-7}

van Zonneveld et al. (2022 in prep). Forgotten foods for healthy diets in a warmer Africa.
}

\examples{

\dontrun{
library(ggplot2)
library(sf)
library(concaveman)

data(CucurbitaClim)

alata.data <- CucurbitaClim[CucurbitaClim$species == "Cucurbita_palmata", ]

bioc.names <- paste0("bioc", 1:19)

alata.data2 <- alata.data[alata.data$ADM0_A3 == "USA", ]
alata.base <- alata.data2[alata.data2$climate == "baseline", bioc.names]
alata.fut  <- alata.data2[alata.data2$climate == "future", bioc.names]

conc2.res <- ensemble.concave.hull(baseline.data=alata.base,
                                  change.data=alata.fut,
                                  method="pca",
                                  VIF.max=40,
                                  concavity=2)

plot(conc2.res$ggplot.out)
conc2.res$baseline.area
conc2.res$change.area
conc2.res$novel.area
conc2.res$novel.area / conc2.res$change.area

# Which accessions have novel climates?
summary(conc2.res$change.points)
change.points <- conc2.res$change.points
rownames(change.points[change.points$novel == TRUE, ])
nrow(change.points[change.points$novel == TRUE, ]) / nrow(change.points)

# Analysis via convex hulls
conc100.res <- ensemble.concave.hull(baseline.data=alata.base,
                                  change.data=alata.fut,
                                  method="pca",
                                  concavity=100)

plot(conc100.res$ggplot.out)
conc100.res$baseline.area
conc100.res$change.area
conc100.res$novel.area
conc100.res$novel.area / conc100.res$change.area

# Which accessions have novel climates?
summary(conc100.res$change.points)
change.points <- conc100.res$change.points
rownames(change.points[change.points$novel == TRUE, ])
nrow(change.points[change.points$novel == TRUE, ]) / nrow(change.points)

# Checking niche overlaps with other accessions
# Alternative 1: niche overlap with accessions from Mexico
alata.data2 <- alata.data[alata.data$ADM0_A3 == "MEX", ]
alata.MEX <- alata.data2[alata.data2$climate == "baseline", bioc.names]

venn2.res <- ensemble.concave.venn(conc2.res,
                                   candidate.data=alata.MEX,
                                   concavity=2)
plot(venn2.res$ggplot.out)
table(venn2.res$change.points[ , c("novel", "candidate.in")])


# alternative 1 for convex hulls
venn100.res <- ensemble.concave.venn(conc100.res,
                                   candidate.data=alata.MEX,
                                   concavity=100)
plot(venn100.res$ggplot.out)
table(venn100.res$change.points[ , c("novel", "candidate.in")])

# alternative 2: niche overlap with other species
cucurbita2 <- CucurbitaClim[CucurbitaClim$climate == "baseline", ]
cordata.data <- cucurbita2[cucurbita2$species == "Cucurbita_cordata", bioc.names]
digitata.data <- cucurbita2[cucurbita2$species == "Cucurbita_digitata", bioc.names]

venn.cordata <- ensemble.concave.venn(conc2.res,
                                      candidate.data=cordata.data,
                                      concavity=2)
plot(venn.cordata$ggplot.out)

venn.digitata <- ensemble.concave.venn(conc2.res,
                                      candidate.data=digitata.data,
                                      concavity=2)
plot(venn.digitata$ggplot.out)

# check the union of the two species
spec.res <- vector("list", 2)
spec.res[[1]] <- venn.cordata
spec.res[[2]] <- venn.digitata
union2.res <- ensemble.concave.union(conc2.res,
                                     candidate.venns=spec.res)
table(union2.res$change.points[ , c("novel", "candidate.in")])

# Analysis via convex hulls
venn.digitata <- ensemble.concave.venn(conc100.res,
                                      candidate.data=digitata.data,
                                      concavity=100)

venn.cordata <- ensemble.concave.venn(conc100.res,
                                      candidate.data=cordata.data,
                                      concavity=100)
spec.res <- vector("list", 2)
spec.res[[1]] <- venn.cordata
spec.res[[2]] <- venn.digitata

union100.res <- ensemble.concave.union(conc100.res,
                                     candidate.venns=spec.res)
plot(union100.res$ggplot.out)
table(union100.res$change.points[ , c("novel", "candidate.in")])

}

}

