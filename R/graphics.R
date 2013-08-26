# 
# @filename: graphics.R 
# 
# @date: 28.03.2013
#
# @author: Stephan Gade <stephan.gade@germanbreastgroup.de>
#
# various plotting functions
#


#'
#' Align ggplo2 plots vertically
#'
#' @param ... ggplo2 plots
#' @param heights the heights of the single plots. If NULL an equal height of 1 is asusmed for all plots. Default to NULL
#' 
#' @return arranged plots
#'
#' @import ggplot2
#' @import gtable
#' @export
AlignPlots <- function(..., heights=NULL) {

    ## get plots
    plots <- list(...)

    ## get tables from all plots
    gtables <- lapply(plots, function(xx)  ggplot_gtable(ggplot_build(xx)))

    ## get maximal width
    width.max <- do.call(unit.pmax, lapply(gtables, function(xx) xx$widths[2:3]))

    ## set maximal with for all objects
    gtables <- lapply(gtables, 
                      function(xx) {
                          xx$widths[2:3] <- width.max
                          xx
                      })

    ## set height to unuiform height of 1 if not given
    ## and create unit classes
    if (is.null(heights)) {
        heights <- rep(1, length(plots))
    }
    heights <- unit(heights, rep("null", length(heights)))

    ret <- do.call(arrangeGrob, c(gtables, list(nrow=length(plots), ncol=1, clip=T, heights=heights)))

    return(ret)
}



#'
#' Kaplan Meier Plot from Matt Cooper (http://mcfromnz.wordpress.com/2012/05/05/kaplan-meier-survival-plot-with-at-risk-table-by-sub-groups/)
#'
#' @author Matt Cooper
#' \url{http://statbandit.wordpress.com/2011/03/08/an-enhanced-kaplan-meier-plot/}
#' @param sfit a \code{\link[survival]{survfit}} object
#' @param sdiff a \code{\link[survival]{survdiff}} object, used for compuation of exact p-values. If no survdiff object is given computation of p-values is based on the survfit object
#' @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#' @param returns logical: if \code{TRUE}, return an arrangeGrob object
#' @param plot logical: if \code{TRUE}, plot the graph
#' @param color logical: if \code{TRUE}, use a colored scheme, otherwise black & white
#' @param xlabs x-axis label
#' @param ylabs y-axis label
#' @param xlims the limit of x axis
#' @param ylims the limit on the y axis
#' @param ystratalabs The strata labels. \code{Default = levels(summary(sfit)$strata)}
#' @param ystrataname The legend name. Default = "Strata"
#' @param timeby numeric: control the granularity along the time-axis
#' @param main plot title
#' @param pval logical: add the pvalue to the plot?
#' @param pval.size size of the pvalue in the plot. Default to 3
#' @param probs numeric: vector with survail probabilities that will be marked by dotted lines, dfault to NULL
#' @param subs default to NULL
#'
#' @import survival
#' @import gridExtra
#' @import plyr
#' @export
#'
ggkm <- function(sfit,
                 sdiff=NULL,
                 table = TRUE,
                 returns = FALSE,
                 plot=T,
                 color = TRUE,
                 xlabs = "Time",
                 ylabs = "Survival Probability",
                 xlims = c(0,max(sfit$time)),
                 ylims = c(0,1),
                 ystratalabs = NULL,
                 ystrataname = NULL,
                 timeby = 100,
                 main = "Kaplan-Meier Plot",
                 pval = TRUE,
                 pval.size=4,
                 probs = NULL,
                 subs = NULL,
                 ...) {

    #############
    # libraries #
    #############

    #     require(ggplot2)
    #     require(survival)
    #     require(gridExtra)

    #################################
    # sorting the use of subsetting #
    #################################

    times <- seq(0, max(sfit$time), by = timeby)

    if(is.null(subs)){
        subs1 <- 1:length(levels(summary(sfit)$strata))
        subs2 <- 1:length(summary(sfit,censored=T)$strata)
        subs3 <- 1:length(summary(sfit,times = times,extend = TRUE)$strata)
    } else{
        for(i in 1:length(subs)){
            if(i==1){
                ssvar <- paste("(?=.*\\b=",subs[i],sep="")
            }
            if(i==length(subs)){
                ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],"\\b)",sep="")
            }
            if(!i %in% c(1, length(subs))){
                ssvar <- paste(ssvar,"\\b)(?=.*\\b=",subs[i],sep="")
            }
            if(i==1 & i==length(subs)){
                ssvar <- paste("(?=.*\\b=",subs[i],"\\b)",sep="")
            }
        }
        subs1 <- which(regexpr(ssvar,levels(summary(sfit)$strata), perl=T)!=-1)
        subs2 <- which(regexpr(ssvar,summary(sfit,censored=T)$strata, perl=T)!=-1)
        subs3 <- which(regexpr(ssvar,summary(sfit,times = times,extend = TRUE)$strata, perl=T)!=-1)
    }

    if( !is.null(subs) ) pval <- FALSE

    ##################################
    # data manipulation pre-plotting #
    ##################################

    if(is.null(ystratalabs)) ystratalabs <- as.character(sub("group=*","",names(sfit$strata))) #[subs1]
    if(is.null(ystrataname)) ystrataname <- "Strata"
    m <- max(nchar(ystratalabs))
    times <- seq(0, max(sfit$time), by = timeby)

    .df <- data.frame(                      # data to be used in the survival plot
                      time = sfit$time[subs2],
                      n.risk = sfit$n.risk[subs2],
                      n.event = sfit$n.event[subs2],
                      surv = sfit$surv[subs2],
                      strata = factor(summary(sfit, censored = T)$strata[subs2]),
                      upper = sfit$upper[subs2],
                      lower = sfit$lower[subs2]
                      )

    levels(.df$strata) <- ystratalabs       # final changes to data for survival plot
    zeros <- data.frame(time = 0, surv = 1,
                        strata = factor(ystratalabs, levels=levels(.df$strata)),
                        upper = 1, lower = 1)
    .df <- rbind.fill(zeros, .df)
    d <- length(levels(.df$strata))

    ###################################
    # specifying plot parameteres etc #
    ###################################

    ## use own environment
    ## if probs is given prob is evaluated in this environment and can be found
    .e <- environment()

    ## assemble plot
    p <- ggplot( .df, aes(time, surv), environment=.e)
    if(color) {
        p <- p + 
        geom_step(aes(colour = strata), size = 0.7) +
        theme_grey()+
        labs(colour = ystrataname)
        if (!is.null(probs)) {
            for (prob in probs)
                p <- p + geom_hline(aes(yintercept=prob), colour="#990000", linetype="dashed")
        }
    } else {
        p <- p + 
        geom_step(aes(linetype = strata), size = 0.7) +
        theme_bw()+
        labs(linetype = ystrataname)
        if (!is.null(probs)) {
            for (prob in probs)
                p <- p + geom_hline(aes(yintercept=prob), linetype="dashed")
        }
    }
    p <- p + 
    theme(legend.position = c(ifelse(m < 10, ylims.28, .35),ifelse(d < 4, max(ylims[1]+.25, 0.8), max(ylims[1]+.35, 0.8)))) +    # MOVE LEGEND HERE [first is x dim, second is y dim]
    theme(legend.key = element_rect(colour = NA))+
    theme(axis.title.x = element_text(vjust = 0.5)) +
    theme(axis.title.y = element_text(vjust=1, angle=90)) + 
    scale_x_continuous(xlabs, breaks = times, limits = xlims) +
    scale_y_continuous(ylabs, limits = ylims) +
    theme(panel.grid.minor = element_blank()) +
    #theme(plot.margin = unit(c(0, 1, .5,ifelse(m < 10, 1.5, 4.5)),"lines")) +
    theme(plot.margin = unit(c(1,1,1.5,0.5), "lines")) + 
    ggtitle(main)


    #####################
    # p-value placement #
    #####################a

    if(pval) {
        if (is.null(sdiff)) {
            sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
        }
        pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
        pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 2)))
        p <- p + annotate("text",x = 0.8 * max(sfit$time),y = (ylims[1] + 0.1),label = pvaltxt, size=pval.size)
    }

    ###################################################
    # Create table graphic to include at-risk numbers #
    ###################################################

    if(table) {
        risk.data <- data.frame(
                                strata = factor(summary(sfit,times = times,extend = TRUE)$strata[subs3], levels=names(sfit$strata)),
                                time = summary(sfit,times = times,extend = TRUE)$time[subs3],
                                n.risk = summary(sfit,times = times,extend = TRUE)$n.risk[subs3]
                                )
        risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))

        data.table <- ggplot(risk.data,aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
        geom_text(size = 4) + theme_bw() +
        scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = rev(ystratalabs)) +
        scale_x_continuous("Numbers at risk", limits = xlims) +
        theme(axis.title.x = element_text(size = 11, vjust = 1),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),axis.text.x = element_blank(),
              axis.ticks = element_blank(),axis.text.y = element_text(face = "bold", size=11, hjust = 1))

        data.table <- data.table +
        theme(legend.position = "none") + xlab(NULL) + ylab(NULL) + 
        theme(plot.margin = unit(c(0,1,0.5,0.5),"lines"))


        #######################
        # Plotting the graphs #
        #######################

        a <- AlignPlots(p, data.table, heights=c(2, .5))
        
        if (plot) {
            print(a)
        }

        if(returns) {
            ret <- list(km=a, main=p, table=data.table)
            attr(ret, "class") <- "ggkm"
            return(ret)
        }
    } else {

        if(plot) print(p)
        if(returns) return(p)
    }
}

#'
#' Print function for ggkm plots
#'
#' A helper function for seamless call of print to results of a ggkm function
#'
#' @param x object of class 'ggkm'
#'
#' @import ggplot2
#' @export
#'
print.ggkm <- function(x) {

    print(x$km)
}




#'
#' @title SplitGraphs
#' Splits list of ggplot2 plots into several grids od define size and plots them
#'
#' @param plot.list the list with gglot2 graphs
#' @param nrow, the number of rows for one grid, default to 3
#' @param ncol, the number of columns of one grid, default to 2
#' @param ... additional arguments to grid.arrange
#'
#' @import ggplot2
#' @import gridExtra
#'
#' @export
#'
#' @author Stephan Gade
#'
SplitGraphs <- function(plot.list, nrow=3, ncol=2,...) {
  
  ## get maximal number of plots in one grid
  num <- nrow*ncol
  
  ## get index
  index <- SplitVec(1:length(plot.list), ceiling(length(plot.list)/num))
  
  ## for every entry in index create a plot
  for (i in 1:length(index)) {
    
    args.list <- c(plot.list[index[[i]]], list(nrow=nrow, ncol=ncol, ...))
    do.call(grid.arrange, args.list)
  }
}



#'
#' Creates a survival data frame 
#'
#' The functions creates an data frame needed for \link{qplot_survival}
#'
#' @author R. Saccilotto
#'
#' @param s.survfit a survival fit object
#'
#' @return a data frame with survival data
#'
#' @export
#'
CreateSurvivalFrame <- function(f.survfit){

    # initialise frame variable
    f.frame <- NULL

    # check if more then one strata
    if(length(names(f.survfit$strata)) == 0){

        # create data.frame with data from survfit
        f.frame <- data.frame(time=f.survfit$time, 
                              n.risk=f.survfit$n.risk, 
                              n.event=f.survfit$n.event, 
                              n.censor = f.survfit$n.censor, 
                              surv=f.survfit$surv, 
                              upper=f.survfit$upper, 
                              lower=f.survfit$lower)

        # add first row to dataset (start at 1)
        f.frame <- rbind(data.frame(time=c(0, f.frame$time[1]), 
                                    n.risk=c(f.survfit$n, f.survfit$n), 
                                    n.event=c(0,0),
                                    n.censor=c(0,0), 
                                    surv=c(1,1), 
                                    upper=c(1,1), 
                                    lower=c(1,1)),
                         f.frame)

    } else { 	 	 
        
        # create vector for strata identification
        f.strata <- NULL
        
        for (f.i in 1:length(f.survfit$strata)){
            # add vector for one strata according to number of rows of strata
            f.strata <- c(f.strata, rep(names(f.survfit$strata)[f.i], f.survfit$strata[f.i]))
        } 	 	 
        
        # create data.frame with data from survfit (create column for strata)
        f.frame <- data.frame(time=f.survfit$time, 
                              n.risk=f.survfit$n.risk, 
                              n.event=f.survfit$n.event, 
                              n.censor = f.survfit$n.censor, 
                              surv=f.survfit$surv, 
                              upper=f.survfit$upper, 
                              lower=f.survfit$lower, 
                              strata=factor(f.strata))

        # create first two rows (start at 1) for each strata
        for (f.i in 1:length(f.survfit$strata)) {

            # take only subset for this strata from data
            f.subset <- subset(f.frame, strata==names(f.survfit$strata)[f.i])

            # create first two rows (time: 0, time of first event)
            # add first two rows to dataset
            f.frame <- rbind(data.frame(time=c(0, f.subset$time[1]), 
                                        n.risk=rep(f.survfit[f.i]$n, 2), 
                                        n.event=c(0,0),
                                        n.censor=c(0,0), 
                                        surv=c(1,1), 
                                        upper=c(1,1), 
                                        lower=c(1,1), 
                                        strata=rep(names(f.survfit$strata)[f.i],2)),
                             f.frame)  

        }

        # reorder data
        f.frame <- f.frame[order(f.frame$strata, f.frame$time), ]

        # rename row.names
        rownames(f.frame) <- NULL
    } 	 

    # return frame
    return(f.frame)
}

#'
#' Survival (KM) plot for ggplot2
#'
#' @author R. Saccilotto
#'
#' @param f.frame a survival data frame constructed by \link{CreateSurvivalFrame}
#' @param f.ci wheter to draw confidence interval: TRUE, FALSE or 'default'. With 'default' it is dependent wheter a strata is given: no with strata, yes without
#' @param f.shape the shape of censoring indicators
#' 
#' @import ggplot2
#' @export 
#' 
qplot_survival <- function(f.frame, f.CI="default", f.shape=3){
  
  # use different plotting commands dependig whether or not strata's are given
  if("strata" %in% names(f.frame) == FALSE){
    
    # confidence intervals are drawn if not specified otherwise
    if (f.CI=="default" | f.CI==TRUE ) {
      
      # create plot with 4 layers (first 3 layers only events, last layer only censored)
      # hint: censoring data for multiple censoring events at timepoint are overplotted
      
      # (unlike in plot.survfit in survival package)
      
      ggplot(data=f.frame) + 
        geom_step(aes(x=time, y=surv), direction="hv") + 
        geom_step(aes(x=time,y=upper), directions="hv", linetype=2) + geom_step(aes(x=time,y=lower), direction="hv", linetype=2) +
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
      
    } else {
      
        # create plot without confidence intervalls
      ggplot(data=f.frame) + 
        geom_step(aes(x=time, y=surv), direction="hv") +
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)	
      
    }	
  } else { 	 	 
    
    if (f.CI=="default" | f.CI==FALSE) {
      
      # without CI
      ggplot(data=f.frame, aes(group=strata, colour=strata)) + 
        geom_step(aes(x=time, y=surv), direction="hv") + 
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
      
    } else {
      
      # with CI (hint: use alpha for CI)
      ggplot(data=f.frame, aes(colour=strata, group=strata)) + 
        geom_step(aes(x=time, y=surv),direction="hv") + 
        geom_step(aes(x=time, y=upper), directions="hv", linetype=2, alpha=0.5) +
        geom_step(aes(x=time,y=lower), direction="hv", linetype=2, alpha=0.5) +
        geom_point(data=subset(f.frame, n.censor==1), aes(x=time, y=surv), shape=f.shape)
    }
  }
}
