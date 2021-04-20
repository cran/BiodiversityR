\name{sites.long}
\alias{sites.long}
\alias{species.long}
\alias{centroids.long}
\alias{vectorfit.long}
\alias{ordisurfgrid.long}
\alias{ordiellipse.long}
\alias{pvclust.long}
\alias{axis.long}
\alias{accumcomp.long}
\alias{renyicomp.long}
\alias{renyi.long}

\title{Helper Functions to Prepare Plotting of Accumulation, Diversity Profile and Ordiplot Results via ggplot2}

\description{
These functions organize outputs from \code{\link[vegan]{ordiplot}}, \code{\link[BiodiversityR]{accumcomp}} and \code{\link[BiodiversityR]{renyicomp}} so these can be plotted with \code{\link[ggplot2]{ggplot}}.
}

\usage{
sites.long(x, env.data = NULL)

species.long(x, spec.data = NULL)

centroids.long(y, grouping, FUN = mean, centroids.only = FALSE)

vectorfit.long(z)

ordisurfgrid.long(z)

ordiellipse.long(z, grouping.name = "Grouping")

pvclust.long(cl, cl.data)

axis.long(w, choices = c(1, 2), cmdscale.model=FALSE, CAPdiscrim.model=FALSE)

accumcomp.long(x, ci = 2, label.freq = 1)

renyicomp.long(x, label.freq = 1)

renyi.long(x, env.data=NULL, label.freq = 1)
}

\arguments{
  \item{x}{ Result of \code{\link[vegan]{ordiplot}}, \code{\link[BiodiversityR]{accumcomp}} or \code{\link[BiodiversityR]{renyicomp}} }
  \item{env.data}{ Environmental descriptors for each site. }
  \item{spec.data}{ Descriptors for each species. }
  \item{y}{ Result of function {\link[BiodiversityR]{sites.long}}. }
  \item{grouping}{ Variable defining the centroids }
  \item{FUN}{ A function to compute the summary statistics which can be applied to all data subsets, as in \code{\link[stats]{aggregate}} }
  \item{centroids.only}{ Return the coordinates for the centroids }
  \item{z}{ Result of \code{\link[vegan]{envfit}}, \code{\link[vegan]{ordisurf}} or \code{\link[vegan]{ordiellipse}}  }
  \item{grouping.name}{ Name for the categorical variable, expected as the factor used in the earlier ordiellipse call. }
  \item{cl}{ Result of \code{\link[pvclust]{pvclust}}  }
  \item{cl.data}{ Result of \code{\link[vegan]{ordicluster}}  }  
  \item{w}{ Ordination object from which the \code{\link[vegan]{ordiplot}} was obtained, expected to be fitted in vegan.  }
  \item{choices}{ Ordination axes selected, as in \code{\link[vegan]{ordiplot}} }
  \item{cmdscale.model}{ Use \code{TRUE} is the model was fitted via \code{\link[stats]{cmdscale}} }
  \item{CAPdiscrim.model}{ Use \code{TRUE} is the model was fitted via \code{\link[BiodiversityR]{CAPdiscrim}} }
  \item{ci}{ Multiplier for confidence intervals as in \code{\link[vegan]{specaccum}}. In case 'NA' is provided, then the multiplier is calculated via \code{\link[stats]{qt}}. }
  \item{label.freq}{ Frequency of labels along the x-axis (count between labels). }
}

\details{
Examples for ordiplot results are shown below.
   
Function \code{pvclust.long} combines data from \code{\link[pvclust]{pvclust}} with coordinates of nodes and branches from \code{\link[vegan]{ordicluster}}. The variable of \code{prune} allows to remove higher levels of nodes and branches in the clustering hierarchy in a similar way as argument \code{prune} for \code{\link[vegan]{ordicluster}} - see examples.
   
See also section: see examples for species accumulation curves and Renyi diversity profiles
}


\value{
These functions produce data.frames that can subsequentially be plotted via \code{ggplot} methods.
}

\references{
Kindt, R. & Coe, R. (2005)
Tree diversity analysis: A manual and
  software for common statistical methods for ecological and
  biodiversity studies.
 
  \url{https://www.worldagroforestry.org/output/tree-diversity-analysis} 

\url{https://rpubs.com/Roeland-KINDT}

}   

\author{Roeland Kindt}

\seealso{
\code{\link{accumcomp}}, \code{\link{renyicomp}}
}

\examples{

\dontrun{
# ggplot2 plotting method
library(ggplot2)
library(ggforce)
library(concaveman)
library(ggrepel)
library(ggsci)
library(dplyr)

library(vegan)
data(dune)
data(dune.env)

attach(dune)
attach(dune.env)

Ordination.model1 <- capscale(dune ~ Management, data=dune.env, 
  distance='kulczynski', sqrt.dist=F, add='cailliez')

plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling='species')

# obtain 'long' data from the ordiplot object
sites1 <- sites.long(plot1, env.data=dune.env)
species1 <- species.long(plot1)
centroids1 <- centroids.long(sites1, Management, FUN=median)
centroids2 <- centroids.long(sites1, Management, FUN=median, centroids.only=TRUE)

# information on percentage variation from the fitted ordination
axislabs <- axis.long(Ordination.model1, choices=c(1 , 2))

# change the theme
# possibly need for extrafont::loadfonts(device="win") to have Arial
# as alternative, use library(ggThemeAssist)
BioR.theme <- theme(
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line("gray25"),
        text = element_text(size = 12, family="Arial"),
        axis.text = element_text(size = 10, colour = "gray25"),
        axis.title = element_text(size = 14, colour = "gray25"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.key = element_blank())

# no species scores
# centroids calculated directly via the centroids.long function 
plotgg1 <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axislabs[1, "label"]) +
    ylab(axislabs[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_mark_hull(data=sites1, aes(x=axis1, y=axis2, colour = Management), 
        concavity = 0.1, alpha=0.8, size=0.2, show.legend=FALSE) +
    geom_point(data=sites1, aes(x=axis1, y=axis2, colour=Management, shape=Management), 
        size=5) +
#    geom_segment(data=species1, aes(x=0, y=0, xend=axis1*2, yend=axis2*2), 
#        size=1.2, arrow=arrow()) +
#    geom_label_repel(data=species1, aes(x=axis1*2, y=axis2*2, label=labels)) +
    geom_point(data=centroids.long(sites1, grouping=Management, centroids.only=TRUE), 
        aes(x=axis1c, y=axis2c, colour=Centroid, shape=Centroid), size=10, show.legend=FALSE) +
    geom_segment(data=centroids.long(sites1, grouping=Management), 
        aes(x=axis1c, y=axis2c, xend=axis1, yend=axis2, colour=Management), 
        size=1, show.legend=FALSE) +
    BioR.theme +
    ggsci::scale_colour_npg() +
    coord_fixed(ratio=1)

plotgg1

# select species to plot based on goodness of fit
spec.envfit <- envfit(plot1, env=dune)
spec.data1 <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species2 <- species.long(plot1, spec.data=spec.data1)
species2 <- species2[species2$r > 0.6, ]

# after_scale introduced in ggplot2 3.3.0
plotgg2 <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axislabs[1, "label"]) +
    ylab(axislabs[2, "label"]) +   
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_mark_ellipse(data=sites1, aes(x=axis1, y=axis2, 
        colour=Management, fill=after_scale(alpha(colour, 0.2))), 
        expand=0, size=0.2, show.legend=TRUE) +
    geom_point(data=sites1, aes(x=axis1, y=axis2, colour=Management, shape=Management), 
        size=5) +
    geom_segment(data=centroids.long(sites1, grouping=Management), 
        aes(x=axis1c, y=axis2c, xend=axis1, yend=axis2, colour=Management), 
        size=1, show.legend=FALSE) +
    geom_segment(data=species2, aes(x=0, y=0, xend=axis1*2, yend=axis2*2), 
        size=1.2, arrow=arrow()) +
    geom_label_repel(data=species2, aes(x=axis1*2, y=axis2*2, label=labels)) +
    BioR.theme +
    ggsci::scale_colour_npg() +
    coord_fixed(ratio=1)

plotgg2

# Add contour and vector for a continuous explanatory variable
Ordination.model2 <- capscale(dune ~ Management, data=dune.env, 
  distance='kulczynski', sqrt.dist=F, add='cailliez')

plot2 <- ordiplot(Ordination.model2, choices=c(1,2), scaling='species')

sites2 <- sites.long(plot2, env.data=dune.env)
axislabs <- axis.long(Ordination.model2, choices=c(1 , 2))

dune.envfit <- envfit(plot2, dune.env)
vectors2 <- vectorfit.long(dune.envfit)

A1.surface <- ordisurf(plot2, y=A1)
A1.surface
A1.grid <- ordisurfgrid.long(A1.surface)

plotgg3 <- ggplot() + 
    geom_contour_filled(data=A1.grid, aes(x=x, y=y, z=z)) +
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axislabs[1, "label"]) +
    ylab(axislabs[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_point(data=sites2, aes(x=axis1, y=axis2, size=A1), shape=21, colour="black", fill="red") +
    geom_segment(data=subset(vectors2, vector=A1), aes(x=0, y=0, xend=axis1*2, yend=axis2*2), 
        size=1.2, arrow=arrow()) +
    geom_label_repel(data=subset(vectors2, vector=A1), aes(x=axis1*2, y=axis2*2, 
        label=vector), size=5) +
    BioR.theme +
    scale_fill_viridis_d() +
    scale_size(range=c(1, 20)) +
    labs(fill="A1") +
    coord_fixed(ratio=1)

plotgg3

# after_stat introduced in ggplot2 3.3.0
plotgg4 <- ggplot() + 
    geom_contour(data=A1.grid, aes(x=x, y=y, z=z, colour=factor(after_stat(level))), size=2) +
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axislabs[1, "label"]) +
    ylab(axislabs[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_point(data=sites2, aes(x=axis1, y=axis2, size=A1), shape=21, colour="black", fill="red") +
    geom_label_repel(data=sites2, aes(x=axis1, y=axis2, label=labels), 
        colour='black', size=4) +
    BioR.theme +
    scale_colour_viridis_d() +
    scale_size(range=c(1, 20)) +
    labs(colour="A1") +
    coord_fixed(ratio=1)

plotgg4

# example of Voronoi segmentation
plotgg5 <- ggplot(data=sites2, aes(x=axis1, y=axis2)) +
    geom_voronoi_tile(aes(fill = Management, group=-1L), max.radius=0.2) +
    geom_voronoi_segment(colour="grey50") +
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axislabs[1, "label"]) +
    ylab(axislabs[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
    geom_point() +
    BioR.theme +
    ggsci::scale_colour_npg() +
    coord_fixed(ratio=1)

plotgg5

# adding ellipse via ordiellipse

plot3 <- ordiplot(Ordination.model1, choices=c(1,2), scaling='species')
axislabs <- axis.long(Ordination.model1, choices=c(1 , 2))

Management.ellipses <- ordiellipse(plot3, groups=Management, display="sites", kind="sd")
Management.ellipses.long2 <- ordiellipse.long(Management.ellipses, grouping.name="Management")

plotgg6 <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axislabs[1, "label"]) +
    ylab(axislabs[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    geom_polygon(data=Management.ellipses.long2, 
                   aes(x=axis1, y=axis2, colour=Management, 
                       fill=after_scale(alpha(colour, 0.2))), 
              size=0.2, show.legend=FALSE) +
    geom_point(data=sites1, aes(x=axis1, y=axis2, colour=Management, shape=Management), 
        size=5) +
    geom_segment(data=centroids.long(sites1, grouping=Management), 
        aes(x=axis1c, y=axis2c, xend=axis1, yend=axis2, colour=Management), 
        size=1, show.legend=FALSE) +
    BioR.theme +
    ggsci::scale_colour_npg() +
    coord_fixed(ratio=1)

plotgg6

# adding cluster results via pvclust.long

library(pvclust)
# transformation as pvclust works with Euclidean distance
dune.Hellinger <- disttransform(dune, method='hellinger')
dune.pv <- pvclust(t(dune.Hellinger), 
                   method.hclust="mcquitty",
                   method.dist="euclidean",
                   nboot=1000)

plot(dune.pv)
pvrect(dune.pv, alpha=0.89, pv="au")

# Model fitted earlier
plot1 <- ordiplot(Ordination.model1, choices=c(1,2), scaling='species')
cl.data1 <- ordicluster(plot1, cluster=as.hclust(dune.pv$hclust))

sites1 <- sites.long(plot1, env.data=dune.env)
axislabs <- axis.long(Ordination.model1, choices=c(1 , 2))

cl.data1 <- ordicluster(plot2, cluster=as.hclust(dune.pv$hclust))
pvlong <- pvclust.long(dune.pv, cl.data1)

# as in example for ordicluster, prune higher level hierarchies
plotgg7 <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab(axislabs[1, "label"]) +
    ylab(axislabs[2, "label"]) +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +   
    geom_segment(data=subset(pvlong$segments, pvlong$segments$prune > 3),
               aes(x=x1, y=y1, xend=x2, yend=y2, colour=au>=0.89, 
                   size=au),
               show.legend=TRUE) +
    geom_point(data=subset(pvlong$nodes, pvlong$nodes$prune > 3), 
               aes(x=x, y=y, fill=au>=0.89), 
               shape=21, size=2, colour="black") +
    geom_point(data=sites1, 
               aes(x=axis1, y=axis2, shape=Management), 
               colour="darkolivegreen4", alpha=0.9, size=5) +
    geom_text(data=sites1,
              aes(x=axis1, y=axis2, label=labels)) +
    BioR.theme +
    ggsci::scale_colour_npg() +
    scale_size(range=c(0.3, 2)) +
    scale_shape_manual(values=c(15, 16, 17, 18)) +
    guides(shape = guide_legend(override.aes = list(linetype = 0))) +
    coord_fixed(ratio=1)

plotgg7

} # dontrun
}

\keyword{multivariate} 


