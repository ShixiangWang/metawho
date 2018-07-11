##'@author Shixiang Wang <w_shixiang@163.com>
##'@rdname log_ES
##'
##'@title log transformation for effect size estimation according to confidence level and distribution
##'
##'@description A variety of different outcome measures which used in meta-analysis,
##'as input, are in the form of log, such as Hazard Ratio. This function
##'is used to do log transformation to calculate effect size and
##'standard error. Then the result can be easier used for model fit.
##'
##'@param data  a dataframe include effect size, ci.lb, ci.ub in the first, second and third columns.
##'@param distribution a character specify distribution. N for normal, t for student and F for F distribution. Default is N for normal distribution.
##'@param conf_level a number specify confidence level, default is 0.05.
##'@param df1 a number specify degree of freedom for t distribution or F distribution (the first one).
##'@param df2 a number specify the second degree of freedom for F distribution.
##'@param var default is FALSE. If TRUE, the sampling variance will be computed.
##'
##'@examples hr <- c(3.12, 1.15) # specify hazard ratios (hr)
##'ci.lb <- c(2.2, 1.03) # specify lower bound for hr confidence intervals
##'ci.ub <- c(4.1, 2.6) # specify upper bound for hr confidence intervals
##' dat = data.frame(hr=hr, ci.lb=ci.lb, ci.ub=ci.ub)
##' log_ES(dat)
##'
##'@export
##'
log_ES = function(data, distribution=c("N", "t", "F"), conf_level=0.05, df1=Inf, df2=NULL, var=FALSE){

    distribution = match.arg(distribution)
    message("Please check the effect size, low/upper confidence interval bounds should in 1,2,3rd columns of input data.")
    hr = data[,1]
    ci.lb = data[, 2]
    ci.ub = data[, 3]

    # ci.lb <- yi - qnorm(level/2, lower.tail = FALSE) * sqrt(vi)
    # ci.ub <- yi + qnorm(level/2, lower.tail = FALSE) * sqrt(vi)

    if(distribution == "N"){
        conf_q = qnorm(1 - conf_level/2)
        yi  = log(hr)
        sei  = (log(ci.ub) - log(ci.lb)) / (2*conf_q)

    }else if (distribution == "t"){
        conf_q = qt(1 - conf_level/2, df = df1)
        yi  = log(hr)
        sei  = (log(ci.ub) - log(ci.lb)) / (2*conf_q)
    }else {
        conf_q = qf(1 - conf_level/2, df = df1, df2 = df2)
        yi  = log(hr)
        sei  = (log(ci.ub) - log(ci.lb)) / (2*conf_q)
    }

    if (!var){
        dat = data.frame(yi=yi, sei=sei)
        dat = cbind.data.frame(data, dat)
    }else{
        vi = sei^2
        dat = data.frame(yi=yi, sei=sei, vi=vi)
        dat = cbind.data.frame(data, dat)
    }


    return(dat)
}


##' @author Shixiang Wang <w_shixiang@163.com>
##' @rdname easy_forest
##'
##' @title easy forest plots
##'
##' @description A wrapper of forest function in metafor package.
##' More easier use for subgroup analysis.
##'
##' @details This is a generic function build on forest function in metafor
##' package. Even though the forest function is flexible, it is not very easy
##' generate plots readily for publication, especially in subgroup analysis.
##'
##' @export
##'
easy_forest = function(){
    UseMethod("easy_forest")
}


##' @rdname easy_forest
##'
##' @title Function to easily create forest plot for object of class "rma"
##'
##' @inheritParams metafor::forest.rma
##' @export
##'
easy_forest.rma = function(x){}

rma(yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i, x2i, t1i, t2i,
    m1i, m2i, sd1i, sd2i, xi, mi, ri, ti, sdi, r2i, ni, mods,
    measure="GEN", intercept=TRUE, data, slab, subset,
    add=1/2, to="only0", drop00=FALSE, vtype="LS",
    method="REML", weighted=TRUE, test="z",
    level=95, digits=4, btt, tau2, verbose=FALSE, control, ...)

dat = escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
             data = dat.bcg, append = TRUE)
res = rma(yi, vi, data = dat)
forest(x, annotate=TRUE, addfit=TRUE, addcred=FALSE, showweights=FALSE,
       xlim, alim, clim, ylim, at, steps=5,
       level=x$level, refline=0, digits=2L, width,
       xlab, slab, mlab, ilab, ilab.xpos, ilab.pos,
       order, transf, atransf, targs, rows,
       efac=1, pch=15, psize, col, border, lty,
       cex, cex.lab, cex.axis, annosym, ...)

forest(res, slab = paste(dat$author, dat$year, sep = ", "),
       xlim = c(-16, 6), at = log(c(0.05, 0.25, 1, 4)), atransf = exp,
       ilab = cbind(dat$tpos, dat$tneg, dat$cpos, dat$cneg),
       ilab.xpos = c(-9.5, -8, -6, -4.5), cex = 0.75)
opar = par(cex = 0.75, font = 2)
text(c(-9.5, -8, -6, -4.5), 15, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75, -5.25), 16, c("Vaccinated", "Control"))
text(-16, 15, "Author(s) and Year", pos = 4)
text(6, 15, "Relative Risk [95% CI]", pos = 2)
par(opar)

##' @rdname easy_forest
##'
##' @title forest plot for object of class "rma", default
##'
##' @inheritParams metafor::forest.rma
##' @export
easy_forest.default = function (x, annotate = TRUE, addfit = TRUE, addcred = FALSE,
          showweights = FALSE, xlim, alim, clim, ylim, at, steps = 5,
          level = x$level, refline = 0, digits = 2L, width, xlab, slab,
          mlab, ilab, ilab.xpos, ilab.pos, order, transf, atransf,
          targs, rows, efac = 1, pch = 15, psize, col, border, lty,
          cex, cex.lab, cex.axis, annosym, ...)
{
    if (!inherits(x, "rma"))
        stop("Argument 'x' must be an object of class \"rma\".")
    if (inherits(x, "rma.ls"))
        stop("Method not yet implemented for objects of class \"rma.ls\". Sorry!")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail",
                              "na.pass")))
        stop("Unknown 'na.action' specified under options().")
    if (missing(transf))
        transf <- FALSE
    if (missing(atransf))
        atransf <- FALSE
    transf.char <- deparse(substitute(transf))
    atransf.char <- deparse(substitute(atransf))
    if (is.function(transf) && is.function(atransf))
        stop("Use either 'transf' or 'atransf' to specify a transformation (not both).")
    if (missing(targs))
        targs <- NULL
    if (missing(at))
        at <- NULL
    if (missing(ilab))
        ilab <- NULL
    if (missing(ilab.xpos))
        ilab.xpos <- NULL
    if (missing(ilab.pos))
        ilab.pos <- NULL
    if (missing(order))
        order <- NULL
    if (missing(psize))
        psize <- NULL
    if (missing(cex))
        cex <- NULL
    if (missing(cex.lab))
        cex.lab <- NULL
    if (missing(cex.axis))
        cex.axis <- NULL
    if (x$int.only) {
        if (missing(col)) {
            col <- c("black", "gray50")
        }
        else {
            if (length(col) == 1L)
                col <- c(col, "gray50")
        }
        if (missing(border))
            border <- "black"
    }
    else {
        if (missing(col))
            col <- "gray"
        if (missing(border))
            border <- "gray"
    }
    if (missing(lty)) {
        lty <- c("solid", "dotted", "solid")
    }
    else {
        if (length(lty) == 1L)
            lty <- c(lty, "dotted", "solid")
        if (length(lty) == 2L)
            lty <- c(lty, "solid")
    }
    if (length(efac) == 1L)
        efac <- rep(efac, 3)
    if (length(efac) == 2L)
        efac <- c(efac[1], efac[1], efac[2])
    if (missing(annosym))
        annosym <- c(" [", ", ", "]")
    if (length(annosym) != 3)
        stop("Argument 'annosym' must be a vector of length 3.")
    measure <- x$measure
    if (inherits(x, "rma.glmm") && showweights)
        stop("Option 'showweights=TRUE' currently not possible for 'rma.glmm' objects. Sorry!")
    if (length(digits) == 1L)
        digits <- c(digits, digits)
    level <- ifelse(level > 1, (100 - level)/100, ifelse(level >
                                                             0.5, 1 - level, level))
    yi <- x$yi.f
    vi <- x$vi.f
    X <- x$X.f
    k <- length(yi)
    if (missing(slab)) {
        if (x$slab.null) {
            slab <- paste("Study", x$slab)
        }
        else {
            slab <- x$slab
        }
    }
    else {
        if (length(slab) == 1 && is.na(slab))
            slab <- rep("", k)
    }
    if (length(yi) != length(slab))
        stop("Number of outcomes does not correspond to the length of the 'slab' argument.")
    if (is.null(dim(ilab)))
        ilab <- cbind(ilab)
    if (length(pch) == 1L)
        pch <- rep(pch, k)
    if (length(pch) != length(yi))
        stop("Number of outcomes does not correspond to the length of the 'pch' argument.")
    options(na.action = "na.pass")
    if (x$int.only) {
        pred <- fitted(x)
        pred.ci.lb <- rep(NA_real_, k)
        pred.ci.ub <- rep(NA_real_, k)
    }
    else {
        temp <- predict(x, level = level)
        pred <- temp$pred
        if (addcred) {
            pred.ci.lb <- temp$cr.lb
            pred.ci.ub <- temp$cr.ub
        }
        else {
            pred.ci.lb <- temp$ci.lb
            pred.ci.ub <- temp$ci.ub
        }
    }
    if (inherits(x, "rma.glmm")) {
        weights <- NULL
    }
    else {
        weights <- weights(x)
    }
    options(na.action = na.act)
    if (!is.null(psize)) {
        if (length(psize) == 1L)
            psize <- rep(psize, k)
        if (length(psize) != length(yi))
            stop("Number of outcomes does not correspond to the length of the 'psize' argument.")
    }
    if (!is.null(order)) {
        if (is.character(order)) {
            if (length(order) != 1)
                stop("Incorrect length of 'order' argument.")
            if (order == "obs")
                sort.vec <- order(yi)
            if (order == "fit")
                sort.vec <- order(pred)
            if (order == "prec")
                sort.vec <- order(vi, yi)
            if (order == "resid")
                sort.vec <- order(yi - pred, yi)
            if (order == "rstandard")
                sort.vec <- order(rstandard(x)$z, yi)
            if (order == "abs.resid")
                sort.vec <- order(abs(yi - pred), yi)
            if (order == "abs.rstandard")
                sort.vec <- order(abs(rstandard(x)$z), yi)
        }
        else {
            sort.vec <- order
        }
        yi <- yi[sort.vec]
        vi <- vi[sort.vec]
        X <- X[sort.vec, , drop = FALSE]
        slab <- slab[sort.vec]
        ilab <- ilab[sort.vec, , drop = FALSE]
        pred <- pred[sort.vec]
        pred.ci.lb <- pred.ci.lb[sort.vec]
        pred.ci.ub <- pred.ci.ub[sort.vec]
        weights <- weights[sort.vec]
        pch <- pch[sort.vec]
        psize <- psize[sort.vec]
    }
    k <- length(yi)
    if (missing(rows)) {
        rows <- k:1
    }
    else {
        if (length(rows) == 1L) {
            rows <- rows:(rows - k + 1)
        }
    }
    if (length(rows) != length(yi))
        stop("Number of outcomes does not correspond to the length of the 'rows' argument.")
    yi <- yi[k:1]
    vi <- vi[k:1]
    X <- X[k:1, , drop = FALSE]
    slab <- slab[k:1]
    ilab <- ilab[k:1, , drop = FALSE]
    pred <- pred[k:1]
    pred.ci.lb <- pred.ci.lb[k:1]
    pred.ci.ub <- pred.ci.ub[k:1]
    weights <- weights[k:1]
    pch <- pch[k:1]
    psize <- psize[k:1]
    rows <- rows[k:1]
    yiviX.na <- is.na(yi) | is.na(vi) | apply(is.na(X), 1, any)
    if (any(yiviX.na)) {
        not.na <- !yiviX.na
        if (na.act == "na.omit") {
            yi <- yi[not.na]
            vi <- vi[not.na]
            X <- X[not.na, , drop = FALSE]
            slab <- slab[not.na]
            ilab <- ilab[not.na, , drop = FALSE]
            pred <- pred[not.na]
            pred.ci.lb <- pred.ci.lb[not.na]
            pred.ci.ub <- pred.ci.ub[not.na]
            weights <- weights[not.na]
            pch <- pch[not.na]
            psize <- psize[not.na]
            rows.new <- rows
            rows.na <- rows[!not.na]
            for (j in seq_len(length(rows.na))) {
                rows.new[rows >= rows.na[j]] <- rows.new[rows >=
                                                             rows.na[j]] - 1
            }
            rows <- rows.new[not.na]
        }
        if (na.act == "na.fail")
            stop("Missing values in results.")
    }
    k <- length(yi)
    ci.lb <- yi - qnorm(level/2, lower.tail = FALSE) * sqrt(vi)
    ci.ub <- yi + qnorm(level/2, lower.tail = FALSE) * sqrt(vi)
    if (is.function(transf)) {
        if (is.null(targs)) {
            yi <- sapply(yi, transf)
            ci.lb <- sapply(ci.lb, transf)
            ci.ub <- sapply(ci.ub, transf)
            pred <- sapply(pred, transf)
            pred.ci.lb <- sapply(pred.ci.lb, transf)
            pred.ci.ub <- sapply(pred.ci.ub, transf)
        }
        else {
            yi <- sapply(yi, transf, targs)
            ci.lb <- sapply(ci.lb, transf, targs)
            ci.ub <- sapply(ci.ub, transf, targs)
            pred <- sapply(pred, transf, targs)
            pred.ci.lb <- sapply(pred.ci.lb, transf, targs)
            pred.ci.ub <- sapply(pred.ci.ub, transf, targs)
        }
    }
    tmp <- .psort(ci.lb, ci.ub)
    ci.lb <- tmp[, 1]
    ci.ub <- tmp[, 2]
    tmp <- .psort(pred.ci.lb, pred.ci.ub)
    pred.ci.lb <- tmp[, 1]
    pred.ci.ub <- tmp[, 2]
    if (!missing(clim)) {
        clim <- sort(clim)
        if (length(clim) != 2L)
            stop("Argument 'clim' must be of length 2.")
        ci.lb[ci.lb < clim[1]] <- clim[1]
        ci.ub[ci.ub > clim[2]] <- clim[2]
        pred.ci.lb[pred.ci.lb < clim[1]] <- clim[1]
        pred.ci.ub[pred.ci.ub > clim[2]] <- clim[2]
    }
    if (is.null(psize)) {
        if (is.null(weights)) {
            if (any(vi <= 0, na.rm = TRUE)) {
                psize <- rep(1, k)
            }
            else {
                wi <- 1/sqrt(vi)
                psize <- wi/sum(wi, na.rm = TRUE)
                psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,
                                                                 na.rm = TRUE) - min(psize, na.rm = TRUE))
                psize <- (psize * 1) + 0.5
                if (all(is.na(psize)))
                    psize <- rep(1, k)
            }
        }
        else {
            wi <- weights
            psize <- wi/sum(wi, na.rm = TRUE)
            psize <- (psize - min(psize, na.rm = TRUE))/(max(psize,
                                                             na.rm = TRUE) - min(psize, na.rm = TRUE))
            psize <- (psize * 1) + 0.5
            if (all(is.na(psize)))
                psize <- rep(1, k)
        }
    }
    rng <- max(ci.ub, na.rm = TRUE) - min(ci.lb, na.rm = TRUE)
    if (annotate) {
        if (showweights) {
            plot.multp.l <- 2
            plot.multp.r <- 2
        }
        else {
            plot.multp.l <- 1.2
            plot.multp.r <- 1.2
        }
    }
    else {
        plot.multp.l <- 1.2
        plot.multp.r <- 0.4
    }
    if (missing(xlim)) {
        xlim <- c(min(ci.lb, na.rm = TRUE) - rng * plot.multp.l,
                  max(ci.ub, na.rm = TRUE) + rng * plot.multp.r)
        xlim <- round(xlim, digits[2])
    }
    alim.spec <- TRUE
    if (missing(alim)) {
        if (is.null(at)) {
            alim <- range(pretty(x = c(min(ci.lb, na.rm = TRUE),
                                       max(ci.ub, na.rm = TRUE)), n = steps - 1))
            alim.spec <- FALSE
        }
        else {
            alim <- range(at)
        }
    }
    alim <- sort(alim)
    xlim <- sort(xlim)
    if (xlim[1] > min(yi, na.rm = TRUE)) {
        xlim[1] <- min(yi, na.rm = TRUE)
    }
    if (xlim[2] < max(yi, na.rm = TRUE)) {
        xlim[2] <- max(yi, na.rm = TRUE)
    }
    if (alim[1] < xlim[1]) {
        xlim[1] <- alim[1]
    }
    if (alim[2] > xlim[2]) {
        xlim[2] <- alim[2]
    }
    if (missing(ylim)) {
        if (x$int.only && addfit) {
            ylim <- c(-1.5, k + 3)
        }
        else {
            ylim <- c(0.5, k + 3)
        }
    }
    else {
        ylim <- sort(ylim)
    }
    if (is.null(at)) {
        if (alim.spec) {
            at <- seq(from = alim[1], to = alim[2], length.out = steps)
        }
        else {
            at <- pretty(x = c(min(ci.lb, na.rm = TRUE), max(ci.ub,
                                                             na.rm = TRUE)), n = steps - 1)
        }
    }
    else {
        at[at < alim[1]] <- alim[1]
        at[at > alim[2]] <- alim[2]
        at <- unique(at)
    }
    at.lab <- at
    if (is.function(atransf)) {
        if (is.null(targs)) {
            at.lab <- formatC(sapply(at.lab, atransf), digits = digits[2],
                              format = "f", drop0trailing = ifelse(class(digits) ==
                                                                       "integer", TRUE, FALSE))
        }
        else {
            at.lab <- formatC(sapply(at.lab, atransf, targs),
                              digits = digits[2], format = "f", drop0trailing = ifelse(class(digits) ==
                                                                                           "integer", TRUE, FALSE))
        }
    }
    else {
        at.lab <- formatC(at.lab, digits = digits[2], format = "f",
                          drop0trailing = ifelse(class(digits) == "integer",
                                                 TRUE, FALSE))
    }
    par.mar <- par("mar")
    par.mar.adj <- par.mar - c(0, 3, 1, 1)
    par.mar.adj[par.mar.adj < 0] <- 0
    par(mar = par.mar.adj)
    on.exit(par(mar = par.mar))
    plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
         yaxt = "n", xaxt = "n", xaxs = "i", bty = "n", ...)
    abline(h = ylim[2] - 2, lty = lty[3], ...)
    if (is.numeric(refline))
        segments(refline, ylim[1] - 5, refline, ylim[2] - 2,
                 lty = "dotted", ...)
    par.usr <- par("usr")
    height <- par.usr[4] - par.usr[3]
    if (is.null(cex)) {
        lheight <- strheight("O")
        cex.adj <- ifelse(k * lheight > height * 0.8, height/(1.25 *
                                                                  k * lheight), 1)
    }
    if (is.null(cex)) {
        cex <- par("cex") * cex.adj
    }
    else {
        if (is.null(cex.lab))
            cex.lab <- cex
        if (is.null(cex.axis))
            cex.axis <- cex
    }
    if (is.null(cex.lab))
        cex.lab <- par("cex") * cex.adj
    if (is.null(cex.axis))
        cex.axis <- par("cex") * cex.adj
    if (addfit && !x$int.only) {
        for (i in seq_len(k)) {
            if (is.na(pred[i]))
                next
            polygon(x = c(max(pred.ci.lb[i], alim[1]), pred[i],
                          min(pred.ci.ub[i], alim[2]), pred[i]), y = c(rows[i],
                                                                       rows[i] + (height/100) * cex * efac[3], rows[i],
                                                                       rows[i] - (height/100) * cex * efac[3]), col = col,
                    border = border, ...)
        }
    }
    if (addfit && x$int.only) {
        if (inherits(x, "rma.mv") && x$withG && x$tau2s > 1) {
            if (!is.logical(addcred)) {
                if (length(addcred) == 1)
                    addcred <- c(addcred, addcred)
                temp <- predict(x, level = level, tau2.levels = addcred[1],
                                gamma2.levels = addcred[2])
                addcred <- TRUE
            }
            else {
                if (addcred) {
                    stop("Need to specify the level of the inner factor(s) via the 'addcred' argument.")
                }
                else {
                    temp <- predict(x, level = level, tau2.levels = 1,
                                    gamma2.levels = 1)
                }
            }
        }
        else {
            temp <- predict(x, level = level)
        }
        beta <- temp$pred
        beta.ci.lb <- temp$ci.lb
        beta.ci.ub <- temp$ci.ub
        beta.cr.lb <- temp$cr.lb
        beta.cr.ub <- temp$cr.ub
        if (is.function(transf)) {
            if (is.null(targs)) {
                beta <- sapply(beta, transf)
                beta.ci.lb <- sapply(beta.ci.lb, transf)
                beta.ci.ub <- sapply(beta.ci.ub, transf)
                beta.cr.lb <- sapply(beta.cr.lb, transf)
                beta.cr.ub <- sapply(beta.cr.ub, transf)
            }
            else {
                beta <- sapply(beta, transf, targs)
                beta.ci.lb <- sapply(beta.ci.lb, transf, targs)
                beta.ci.ub <- sapply(beta.ci.ub, transf, targs)
                beta.cr.lb <- sapply(beta.cr.lb, transf, targs)
                beta.cr.ub <- sapply(beta.cr.ub, transf, targs)
            }
        }
        tmp <- .psort(beta.ci.lb, beta.ci.ub)
        beta.ci.lb <- tmp[, 1]
        beta.ci.ub <- tmp[, 2]
        tmp <- .psort(beta.cr.lb, beta.cr.ub)
        beta.cr.lb <- tmp[, 1]
        beta.cr.ub <- tmp[, 2]
        if (!missing(clim)) {
            beta.ci.lb[beta.ci.lb < clim[1]] <- clim[1]
            beta.ci.ub[beta.ci.ub > clim[2]] <- clim[2]
            beta.cr.lb[beta.cr.lb < clim[1]] <- clim[1]
            beta.cr.ub[beta.cr.ub > clim[2]] <- clim[2]
        }
        if (x$method != "FE" && addcred) {
            segments(max(beta.cr.lb, alim[1]), -1, min(beta.cr.ub,
                                                       alim[2]), -1, lty = lty[2], col = col[2], ...)
            if (beta.cr.lb >= alim[1]) {
                segments(beta.cr.lb, -1 - (height/150) * cex *
                             efac[1], beta.cr.lb, -1 + (height/150) * cex *
                             efac[1], col = col[2], ...)
            }
            else {
                polygon(x = c(alim[1], alim[1] + (1.4/100) *
                                  cex * (xlim[2] - xlim[1]), alim[1] + (1.4/100) *
                                  cex * (xlim[2] - xlim[1]), alim[1]), y = c(-1,
                                                                             -1 + (height/150) * cex * efac[2], -1 - (height/150) *
                                                                                 cex * efac[2], -1), col = col[2], border = col[2],
                        ...)
            }
            if (beta.cr.ub <= alim[2]) {
                segments(beta.cr.ub, -1 - (height/150) * cex *
                             efac[1], beta.cr.ub, -1 + (height/150) * cex *
                             efac[1], col = col[2], ...)
            }
            else {
                polygon(x = c(alim[2], alim[2] - (1.4/100) *
                                  cex * (xlim[2] - xlim[1]), alim[2] - (1.4/100) *
                                  cex * (xlim[2] - xlim[1]), alim[2]), y = c(-1,
                                                                             -1 + (height/150) * cex * efac[2], -1 - (height/150) *
                                                                                 cex * efac[2], -1), col = col[2], border = col[2],
                        ...)
            }
        }
        polygon(x = c(beta.ci.lb, beta, beta.ci.ub, beta), y = c(-1,
                                                                 -1 + (height/100) * cex * efac[3], -1, -1 - (height/100) *
                                                                     cex * efac[3]), col = col[1], border = border,
                ...)
        if (missing(mlab))
            mlab <- ifelse((x$method == "FE"), "FE Model", "RE Model")
        text(xlim[1], -1, mlab, pos = 4, cex = cex, ...)
    }
    axis(side = 1, at = at, labels = at.lab, cex.axis = cex.axis,
         ...)
    if (missing(xlab))
        xlab <- .setlab(measure, transf.char, atransf.char, gentype = 1)
    mtext(xlab, side = 1, at = min(at) + (max(at) - min(at))/2,
          line = par("mgp")[1] - 0.5, cex = cex.lab, ...)
    for (i in seq_len(k)) {
        if (is.na(yi[i]) || is.na(vi[i]))
            next
        if (ci.lb[i] >= alim[2]) {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[2]), y = c(rows[i],
                                                                   rows[i] + (height/150) * cex * efac[2], rows[i] -
                                                                       (height/150) * cex * efac[2], rows[i]), col = "black",
                    ...)
            next
        }
        if (ci.ub[i] <= alim[1]) {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[1]), y = c(rows[i],
                                                                   rows[i] + (height/150) * cex * efac[2], rows[i] -
                                                                       (height/150) * cex * efac[2], rows[i]), col = "black",
                    ...)
            next
        }
        segments(max(ci.lb[i], alim[1]), rows[i], min(ci.ub[i],
                                                      alim[2]), rows[i], lty = lty[1], ...)
        if (ci.lb[i] >= alim[1]) {
            segments(ci.lb[i], rows[i] - (height/150) * cex *
                         efac[1], ci.lb[i], rows[i] + (height/150) * cex *
                         efac[1], ...)
        }
        else {
            polygon(x = c(alim[1], alim[1] + (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[1] + (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[1]), y = c(rows[i],
                                                                   rows[i] + (height/150) * cex * efac[2], rows[i] -
                                                                       (height/150) * cex * efac[2], rows[i]), col = "black",
                    ...)
        }
        if (ci.ub[i] <= alim[2]) {
            segments(ci.ub[i], rows[i] - (height/150) * cex *
                         efac[1], ci.ub[i], rows[i] + (height/150) * cex *
                         efac[1], ...)
        }
        else {
            polygon(x = c(alim[2], alim[2] - (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[2] - (1.4/100) * cex *
                              (xlim[2] - xlim[1]), alim[2]), y = c(rows[i],
                                                                   rows[i] + (height/150) * cex * efac[2], rows[i] -
                                                                       (height/150) * cex * efac[2], rows[i]), col = "black",
                    ...)
        }
    }
    text(xlim[1], rows, slab, pos = 4, cex = cex, ...)
    if (!is.null(ilab)) {
        if (is.null(ilab.xpos))
            stop("Must specify 'ilab.xpos' argument when adding information with 'ilab'.")
        if (length(ilab.xpos) != ncol(ilab))
            stop(paste0("Number of 'ilab' columns (", ncol(ilab),
                        ") does not match length of 'ilab.xpos' argument (",
                        length(ilab.xpos), ")."))
        if (!is.null(ilab.pos) && length(ilab.pos) == 1)
            ilab.pos <- rep(ilab.pos, ncol(ilab))
        for (l in seq_len(ncol(ilab))) {
            text(ilab.xpos[l], rows, ilab[, l], pos = ilab.pos[l],
                 cex = cex, ...)
        }
    }
    if (annotate) {
        if (is.function(atransf)) {
            if (is.null(targs)) {
                if (addfit && x$int.only) {
                    annotext <- cbind(sapply(c(yi, beta), atransf),
                                      sapply(c(ci.lb, beta.ci.lb), atransf), sapply(c(ci.ub,
                                                                                      beta.ci.ub), atransf))
                }
                else {
                    annotext <- cbind(sapply(yi, atransf), sapply(ci.lb,
                                                                  atransf), sapply(ci.ub, atransf))
                }
            }
            else {
                if (addfit && x$int.only) {
                    annotext <- cbind(sapply(c(yi, beta), atransf,
                                             targs), sapply(c(ci.lb, beta.ci.lb), atransf,
                                                            targs), sapply(c(ci.ub, beta.ci.ub), atransf,
                                                                           targs))
                }
                else {
                    annotext <- cbind(sapply(yi, atransf, targs),
                                      sapply(ci.lb, atransf, targs), sapply(ci.ub,
                                                                            atransf, targs))
                }
            }
            tmp <- .psort(annotext[, 2:3])
            annotext[, 2:3] <- tmp
        }
        else {
            if (addfit && x$int.only) {
                annotext <- cbind(c(yi, beta), c(ci.lb, beta.ci.lb),
                                  c(ci.ub, beta.ci.ub))
            }
            else {
                annotext <- cbind(yi, ci.lb, ci.ub)
            }
        }
        if (showweights) {
            if (addfit && x$int.only) {
                annotext <- cbind(c(weights, 100), annotext)
            }
            else {
                annotext <- cbind(weights, annotext)
            }
        }
        annotext <- formatC(annotext, format = "f", digits = digits[1])
        if (missing(width)) {
            width <- apply(annotext, 2, function(x) max(nchar(x)))
        }
        else {
            if (length(width) == 1L)
                width <- rep(width, ncol(annotext))
        }
        for (j in seq_len(ncol(annotext))) {
            annotext[, j] <- formatC(annotext[, j], width = width[j])
        }
        if (showweights) {
            annotext <- cbind(annotext[, 1], "%   ", annotext[,
                                                              2], annosym[1], annotext[, 3], annosym[2], annotext[,
                                                                                                                  4], annosym[3])
        }
        else {
            annotext <- cbind(annotext[, 1], annosym[1], annotext[,
                                                                  2], annosym[2], annotext[, 3], annosym[3])
        }
        annotext <- apply(annotext, 1, paste, collapse = "")
        if (addfit && x$int.only) {
            text(x = xlim[2], c(rows, -1), labels = annotext,
                 pos = 2, cex = cex, ...)
        }
        else {
            text(x = xlim[2], rows, labels = annotext, pos = 2,
                 cex = cex, ...)
        }
    }
    for (i in seq_len(k)) {
        if (is.na(yi[i]))
            next
        if (yi[i] >= alim[1] && yi[i] <= alim[2])
            points(yi[i], rows[i], pch = pch[i], cex = cex *
                       psize[i], ...)
    }
    if (x$int.only && addfit)
        abline(h = 0, lty = lty[3], ...)
    res <- list(xlim = par("usr")[1:2], alim = alim, at = at,
                ylim = ylim, rows = rows, cex = cex, cex.lab = cex.lab,
                cex.axis = cex.axis)
    invisible(res)
}


