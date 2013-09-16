#
# @filename: tools.R
#
# @date:     16 Jul 2012, 05:24pm
#
# @author:   Stephan Gade <stephan.gade@germanbreastgroup.de>
#


#'
#' Helper function to replace factors in coeffecients by factor: level
#'
#' @param object model object 
#'
ReplaceCoefNames <- function(object) {

    ## get coefficients names
    coef.names <- rownames(summary(object)$coefficients)

    if ("coxph" %in% class(object)) {

        ## get the coefficients which are factors
        ttt <- attr(object$terms, "dataClass")
        facs <- names(ttt)[ttt=="factor"]    

        ## replace factor co-variates with factor: level
        for (i in facs) {
            m <- object$assign[[i]]
            coef.names[m] <- sub(i, paste(i,": ", sep=""), coef.names[m])
        }

        return(coef.names)

    } else {
        return(coef.names)
    }

}

#'
#' Function to extract odds-ratios, confidence intervals, and p values as given by the 'summary' function for 
#' 'glm' and related objects
#'
#' @param object a 'glm' object
#' @param level the confidence level (has to be between 0 and 1). Default to 0.95
#' @param ... further arguments to \code{confint}
#'
#' @return matrix with odds-ratios (exponential of coefficients), confidence intervals, and p values
#'
#' @export
#'
ExtractOR <- function(object, level=0.95, ...) {

    ## get the summary and therewith the p.values
    coefs <- summary(object)$coefficients

    ## get odds ratios, simply the exp of coefs
    or <- signif(exp(coef(object)),3)

    ## get CI intervals
    ci <- signif(exp(confint(object, level=level, ...)),3)
    
    ## assemble data frame
    ret <- cbind(or,
                 ci,
                 signif(coefs[,4,drop=F], 2))
    
    ## get level as percentage value
    lev <- paste(format(level*100), "%", sep="")

    ## set col and rownames
    colnames(ret) <- c("OR", paste(lev, "lower"), paste(lev, "upper"), "P")
    rownames(ret) <- rownames(coefs)

    return(ret)
}



#'
#' 'coxph' and related objects
#'
#' @param object a 'coxph' object
#' @param level the confidence level (has to be between 0 and 1). Default to 0.95
#' @param ... further arguments to \code{confint}
#'
#' @return matrix with hazard-ratios (exponential of coefficients), confidence intervals, and p values
#'
#' @export
#'
ExtractHR <- function(object, level=0.95, ...) {

    ## get the summary and therewith the p.values
    coefs <- summary(object)$coefficients

    ## get odds ratios, simply the exp of coefs
    hr <- signif(exp(coef(object)),3)

    ## get CI intervals
    ci <- signif(exp(confint(object, level=level, ...)),3)
    #ci <- paste("[", apply(signif(exp(confint(object)), 3),1, paste, collapse=", "), "]", sep="")
   
    ## assemble data frame
    ret <- cbind(hr,
                 ci,
                 signif(coefs[,5,drop=F], 2))

    ## get level as percentage value
    lev <- paste(format(level*100), "%", sep="")

    ## set col and rownames
    colnames(ret) <- c("HR", paste(lev, "lower"), paste(lev, "upper"), "P")
    rownames(ret) <- rownames(coefs)

    return(ret)
}

#'
#' Prints output of Cox model 
#'
#' The function is wrapper around the \link{texreg} function of the \link{texreg} package which replaces the names of coefficients
#'
#' @param object a coxph model
#' @param ... further argument to texreg
#'
#' @import texreg
#' @export
#'
PrintHR <- function(object, ...) {

    texreg(object, custom.names=ReplaceCoefNames(object), ... )
}


#'
#' function to split a vector in pieces, needed for e.g. cross validation
#' 
#' @param x the vector to split
#' @param n the number of splits
#' 
#' @export
#'
SplitVec <- function(x,n) split(x, factor(sort(rank(x)%%n)))



#'
#' Helper function to read clinical annotation from GEO Matrix Series files
#' 
#' @param file the name of the matrix series file
#' @param start.line the number of the line holding the names of the samples
#' @param end.line the line after (!) the last line to read
#' @param pattern.search a named character vector giving the search pattern for the different clinical parameters, names should be the names of the pararmeters
#' @param pattern.extract character vector holding the pattern for extraction the desired parameter value, set to "" if no extraction is need, vector is recycled to same length as  pattern.search
#' @param drop.first.col boolean, indicates if the first column of the file contains line descriptors and thus should be dropped, default to TRUE
#' @param ... further arguments to sub
#'
#' @return a data frame with the parameters specified in pattern.search
#'
#' @export
#' 
ParseMatrixFile <- function(file, start.line, end.line, pattern.search, pattern.extract="^.*:\\s*", drop.first.col=T, ...) {

    ## read file 
    anno <- read.delim(file, as.is=T, header=F, skip=(start.line-1), nrows=(end-start))

    ## clearning annotation, if needed, drop the first column
    if (drop.first.col) {
        anno <- anno[,-1]
    }

    ## set colnames, the sample names
    colnames(anno) <- anno[1,]
    anno <- anno[-1,]


    ## extract clinical variables
    ## bases on the search pattern

    # all samples
    samples <- colnames(anno)

    # recycle extraxt pattern
    pattern.extract=rep(pattern.extract, length.out=length(pattern.search))

    ## get glinical parmeters
    clin.par <- matrix("", ncol=length(pattern.search), nrow=length(samples))
    colnames(clin.par) <- names(pattern.search)
    rownames(clin.par) <- samples
    for (i  in seq_along(pattern.search)) {
        ttt <- unlist(apply(anno, 2, 
                            function(xx) {
                                as.character(sub(pattern.extract[i], "\\1", xx[grep(pattern.search[i], xx)], ...))
                            }))
        clin.par[,i] <- ttt[samples]
    }

    return(as.data.frame(clin.par, stringsAsFactors=F))

}



#'
#' Normalize Quantiles with another data matrix as template to get the quantiles
#'
#' @param x: the data matrix or vector we want to normalize
#' @param template: a data matrix or vector providing the quantiles
#'
#' @return the data matrix/vector normalized to the same quantiles like the template
#'
#' @export
#'
NormalizeQuantilesTemplate <- function(x, template) {

    # perform some checks

    if(!is.matrix(template) && !is.vector(template)){
	stop("The template has to be either a vector or a matrix!")
    }
    else if(is.matrix(template)) {
	template <- template[,1]
    }

    if(!is.matrix(x) && !is.vector(x)){
	stop("x has to be either a vector or a matrix!")
    }
    else if(is.vector(x)) {
	dim(x) <- c(length(x),1)
    }

   
    # quantile normalization works columnwise
    res <- apply(x, 2, function(xx) {
		
		# get a "mask" for the normalization, this is the sorted template vector
		mask <- sort(template)	
		# order the current colukn ascending
		o <- order(xx)

		# use the order to access xx and set the order xx to the mask
		# note! the actual order of xx wasn't change
		xx[o] <- mask

		# return the column
		return(xx)

	    })


    rownames(res) <- rownames(x)
    colnames(res) <- colnames(x)

    return(res)
}


#'
#' Function to summarize probes of an expression matrix based on common gene identifiers (e.g gene symbols)
#'
#' @param eset the data matrix we want to summarize, the rownames have to be identiefiers in the first column of mappings
#' @param mapping a character vector with the IDs we want to use for summarization
#' 
#' @return the expression matrix with summarized probes
#'
#' @export
#' 
SummarizeProbes <- function(eset, mapping, FUN=median,...) {

    # check if mapping has the same length as eset has rows
    if(length(mapping)!=nrow(eset)) {
        stop("Mapping must have the same length as eset rows!")
    }

    # create new expression ratio matrix
    ttt <- tapply(eset[,1], as.character(mapping), FUN, ...)
    
    exprValues <- matrix(0, length(ttt), ncol(eset))
    colnames(exprValues) <- colnames(eset)

    # summarize the single probes using the median 
    for (i in 1:ncol(exprValues)) { 

        # use FUN to summarize the probes of one gene
        exprValues[,i] <- tapply(eset[,i], as.character(mapping), FUN, ...)

    }
    rownames(exprValues) <- names(ttt)

    return(exprValues)
}



#'
#' A convenient wrapper for limma (moderated t-test)
#'
#' @param data the data matrix, rows corresponds to features (lipds, genes, etc.), columns to samples (patients)
#' @param design the design matrix, in case of a simple t-test only two columns with the two groups 
#' @param contrasts a character vector specifying the contrast to be tested. When using not using eBayes and toptable currently only one contrast is allowed.
#' @param ebayes boolean indicating wheter the resulting fit should be moderated with an empirical bayes approach. If not, the method boils down to a simple linear model (anova in case of groups). Default to TRUE.
#' @param ordinary.F boolean indicating wheter p-values are calculated, similar to the \code{\link{anova.lm}} function, using the ordinary F statistic and degrees of freedom. This is useful if limma is used to compute many t-test or ordinary ANOVAs in a vectorized fashion. Otherwise, a more sophistacated approach from topTable is used to compute the F-tests.
#' @param ... further arguments to \code{\link{toptable}} 
#'
#' @return a data frame with the test results, one for every row, as returned by \code{\link{toptable}}.
#' @seealso \code{\link{model.matrix}}
#' @seealso \code{\link{makeContrasts}}
#'
#' @import limma
#' @export
#'
#' @examples
#' data <- cbind(matrix(rnorm(25, mean=1), nrow=5),
#'               matrix(rnorm(25, mean=3), nrow=5))
#' groups <- factor(rep(c("T", "N"), each=5), levels=c("T", "N"))
#' design <- model.matrix(~0+groups)
#' colnames(design) <- levels(groups)
#' contrasts <- "T-N"
#' # simple vectorized ANOVA with limma
#' res <- PerformLimma(data=data, design=design, contrasts=contrasts, ebayes=FALSE, ordinary.F=TRUE)
#' 
PerformLimma <- function(data, design, contrasts, ebayes=TRUE, ordinary.F=FALSE, ...) {

    # first fit using the data matrix and the design matrix
    fit <- lmFit(data, design)

    
    # make contrasts and using the contrast matrix for the second fit
    if(!is.null(contrasts)) {
        contr <- makeContrasts(contrasts=contrasts, levels=design)
        fit <- contrasts.fit(fit, contr)
    }

    # case 1: ebayes
    if(ebayes)  {

        # moderate variance with empirical bayes
        fit <- eBayes(fit)

        # return toptable with all features
        res <- toptable(fit, number=Inf, ...)


    } else {


        # compute ordinary t-statistics
        fit$t <- fit$coef/fit$stdev.unscaled/fit$sigma 
        fit$t.p.value <- 2*pt(-abs(fit$t), df=fit$df.residual)

        # get F statistic and corresponding p-values
        # try classifyTestsF from limma
        # this code is from the limma function eBayes
        F.stat <- classifyTestsF(fit, fstat.only=T)
        fit$F <- as.vector(F.stat)
        df1 <- attr(F.stat, "df1")
        df2 <- attr(F.stat, "df2")

        if (ordinary.F) {

            if(is.null(contrasts)) {
                warning("No contrasts were given! In case of ANOVA the F statistic might be misleading!")
            }


            if(!is.null(contrasts) && length(contrasts)>1) {
                stop("When not using eBayes and toptable, only one contrast is allowed!")
            }

            # in case of ordinary F statistics
            # use df.residual to calculate p-values for the F statistic
            fit$F.p.value <- pf(fit$F, df1, fit$df.residual, lower.tail = F)

        } else {
            # otherwise use the df from limma's classifyTestsF (FStat)
            # this is the same method used in the toptable function
            if (df2[1] > 1e+06) {
                fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)
            } else {
                fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)
            }
        }

        ## assemble result data frame
        
        # first IDs, if available take rownames of data
        if(is.null(rownames(data))) {
            res <- data.frame(ID=1:nrow(data))
        } else {
            res <- data.frame(ID=rownames(data))
        }
        

        # append logFC
        colnames(fit$coef) <- paste("logFC(", colnames(fit$coef), ")", sep="")
        res <- cbind(res, fit$coef)

        # append t-statistic and corresponding p.values
        colnames(fit$t.p.value) <- paste("p(", colnames(fit$t), ")", sep="")
        colnames(fit$t) <- paste("t(", colnames(fit$t), ")", sep="")
        res <- cbind(res, fit$t, fit$t.p.value)

        # append F statstic and p-values
        res <- cbind(res, F=fit$F, F.p=fit$F.p.value)
    }


    # return the result
    return(res)
}



#'
#' Helper function to calculate Logrank Test p values
#'
#' @note code is from the print.survdiff function in the survival package by Terry Therneau
#'
#' @param x a survfit objejct
#'
#' @return the logrank test p-value
#' 
#' @import survival
#' @export
#' 
LogRankP <- function(x) {

  if (is.matrix(x$obs)) {
    etmp <- apply(x$exp, 1, sum)
  }
  else {
    etmp <- x$exp
  }
  
  df <- (sum(1 * (etmp > 0))) - 1
  p <- 1 - pchisq(x$chisq, df)

  return(p)    
}

#'
#' Helper function to calculate Cox Model p values
#'
#' @note code is from the print.coxph function in the survival package by Terry Therneau
#'
#' @param x a coxph objejct
#'
#' @return the p-values from the cox model
#' 
#' @import survival
#' @export
#' 
CoxP <- function(x) {

    logtest <- -2 * (x$loglik[1] - x$loglik[2])

    if (is.null(x$df)) {
        df <- sum(!is.na(x$coefficients))
    }
    else {
        df <- round(sum(x$df), 2)
    }

    p <- 1 - pchisq(logtest, df)

    return(p)
}


#'
#' Helper function to make Hmisc's latex easier
#'
#' @param data frame or matrix to be printed
#' @param file the file to be written, default to "" which means no file but stdout
#' @param booktabs boolean, whether to use booktabs package or not, default to TRUE
#' @param ctable boolean, whether to use ctable package or not, default to FALSE
#' @param rowlabel label for first column with rownames, default to ""
#' @param n.cgroup see Hmisc::latex, default to NULL. If NULL but cgroup is not will be set to ncol(x)
#' @param cgroup see Hmisc::latex
#' @param ... further arguments to Hmisc::latex
#'
#' @import Hmisc
#' @export
#'
PrintLatex <- function(x, file="", booktabs=T, ctable=F, rowlabel="", n.cgroup=NULL, cgroup=NULL, n.ngroup=NULL, ngroup=NULL, ...) {

    ## sanitize col and row names
    ## not done automatically by Hmisc::latex
    colnames(x) <- sani(colnames(x))
    rownames(x) <- sani(rownames(x))

    ## set appropriate n.cgroup if not given but cgroup
    if (is.null(n.cgroup) && !is.null(cgroup)) {
        if (ncol(x) %% length(cgroup) != 0) {
            stop('ncol(x) is not a multiply of length(cgroup). n.cgroup must be givern!')
        }
        n.cgroup <- rep((ncol(x) %/% length(cgroup)), length(cgroup))
    }

    ## set appropriate n.ngroup if not given but ngroup
    if (is.null(n.ngroup) && !is.null(ngroup)) {
        if (nrow(x) %% length(ngroup) != 0) {
            stop('nrow(x) is not a multiply of length(ngroup). n.ngroup must be givern!')
        }
        n.ngroup <- rep((nrow(x) %/% length(ngroup)), length(ngroup))
    }

    Hmisc::latex(x, file=file, booktabs=booktabs, ctable=ctable, rowlabel=rowlabel, n.cgroup=n.cgroup, cgroup=cgroup, ...)
}


#'
#' CIndexTest
#'
#' Tests for a variable the Improvement in C-Index of a Cox model
#'
#' @param formula model formula with a Surv object on the left side
#' @param var the name of variable to test with C-Index. Will be added to the formula if not included
#' @param data the data frame for model fitting. Must contain all variables in the formula
#' @param N number of permutations for the test, default to 10,000
#' @param cores number of cores to use for permutation
#'
#' @return list with the original C-Index, the permutation C-Indices and the p-value (one sided, alternative="higher")
#'
CIndexTest <- function(formula, var, data, N=10000, seed=123, cores=2, ...) {


    require(survival)
    require(parallel)

    ## crate data frame with needed variables
    mf <- match.call()
    m <- match(c("formula", "data", "subset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")

    ## check if var is already in the formula
    if (!var %in% attr(terms(formula), "term.labels")) {
        cat("CIndexTest: Adding", var, "to model formula ... ")
        formula <- update(formula, as.formula(paste("~.+", var, sep="")))
        mf[[2L]] <- formula
        cat("!\n")
    }

    ## eval mf, check if all variables are there
    mf <- eval(mf, parent.frame())

    ## fit original Cox model and get original C-Index
    cox.ori <- coxph(formula=formula, data=data, ...)
    c.ori   <- summary(cox.ori)$concordance[1]

    ## get permutation
    set.seed(seed)
    perms <- lapply(1:N, function(xx) sample(x=nrow(data), size=nrow(data), replace=F))

    
    ## parallel computation of c-indices
    c.perms <- mclapply(perms, 
                        function(perm) {

                            ## get data with permutated variable
                            data.perm <- data
                            data.perm[[var]] <- data[perm,var]

                            ## get cox model and concordance
                            cox.perm <- coxph(formula=formula, data=data.perm, ...)

                            ## return concordance index
                            return(summary(cox.perm)$concordance[1])

                        }, mc.cores=cores)
    c.perms <- unlist(c.perms)

    ## get cdf and p value one sided
    cdf <- ecdf(c.perms)
    p.value <- 1-cdf(c.ori)

    ## get return object
    ret <- list(c.ori=c.ori, c.perms=c.perms, p.value=p.value, model=formula, var=var, N=N)
    class(ret) <- "citest"
    return(ret)
}

#'
#' Plot function for 'citest' class
#'
plot.citest <- function(x)  {

    ## load libraries
    require(ggplot2)

    data.plot <- data.frame(c.perms=x$c.perms)
    pl <- ggplot(data.plot, aes(x=c.perms)) + geom_histogram() + geom_vline(xintercep=x$c.ori, colour="red") + xlab("C-Index (permutations)") 
    pl <- pl +  annotate("text", label=paste("C-Index: ", signif(x$c.ori, 3)), x=x$c.ori, y=1000)

    print(pl)

    invisible(pl)

}


#'
#' Print function for 'citest' class
#'
#' @param x object of class 'citest'
#'
#' @return x (invisible)
#'
print.citest <- function(x) {

    ## print out information
    cat("\n       Concordance index test for", x$var, "in Cox model\n        ") 
    print(x$model) 
    cat("\n\n")
    cat("  C-Index for full model:", x$c.ori,"\n")
    cat("  p-value:", x$p.value,"\n")
    cat("  number of permutations:", x$N, "\n\n")

    ## return x
    invisible(x)
}
