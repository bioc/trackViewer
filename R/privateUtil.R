convertFont <- function(){
        return(36*min(c(convertWidth(unit(1, "snpc"), "inches", valueOnly=TRUE), 
                        convertHeight(unit(1, "snpc"), "inches", valueOnly=TRUE))
                      ))
}

is_color <- function(x) {
  vapply(x, FUN=function(X) {
    tryCatch(is.matrix(col2rgb(X)), 
             error = function(e) FALSE)
  }, FUN.VALUE = logical(1L))
}

inRange <- function(x, scale){
  x>=scale[1] & x<=scale[2]
}

hrScale <- function(r){
    if(any(r<0))
      stop("'r' must be positive")
    if(length(r)!=2)
        stop("'r' must be a vector with two numbers")
    x <- abs(r[2] - r[1])
    suffix <- c("bp", "K", "M", "G")#NOT USE: "T", "P", "E", "Z", "Y"
    base <- 1000L
    n <- length(suffix)
    for(i in 1:n){
        if(x >= base*10){
            if(i < n)
                x <- x/base
        }else break
    }
    x <- round(x=x/10, digits=0)
    if(i==1) x <- 10^round(log10(x), digits=0) / 2
    return(list(scale=x, unit=suffix[i], range=r))
}

convertNum2HumanNum <- function(x){
  x <- x[1]
  stopifnot(x>0)
  suffix <- c("bp", "K", "M", "G")
  base <- 1000L
  n <- length(suffix)
  for(i in 1:n){
    if(x >= base){
      if(i < n)
        x <- x/base
    }else break
  }
  return(paste(x, suffix[i]))
}

locateScale <- function(x, y, maxY, scale){
    suffix <- c(1, 1000, 1000000, 1000000000)
    names(suffix) <- c("bp", "K", "M", "G")
    rg <- scale[[3]]
    scale <- as.numeric(scale[[1]] * suffix[scale[[2]]])
    scale.10 <- round(x=scale/10, digits=1)
    threshold <- maxY/2
    start <- c()
    end <- c()
    START <- TRUE
    defaultOut <- c(start=mean(rg)-scale, end=mean(rg)+scale)
    if(length(x)==0) return(defaultOut)
    y <- y[order(x)]
    x <- x[order(x)]
    for(i in seq_along(x)){
        if(y[i]>threshold){
            if(START){
                START <- FALSE
                start <- c(start, x[i])
            }
        }else{
            if(!START){
                START <- TRUE
                end <- c(end, x[i])
            }
        }
    }
    if(length(start)==length(end)+1) 
        end <- c(end, x[i])
    if(length(start)==0){
        mid <- sum(rg)/2
        return(c(start=mid-scale, end=mid+scale))
    }
    ir <- IRanges(start=c(x[1],start, x[length(x)])-scale.10, 
                  end=c(x[1], end, x[length(x)])+scale.10)
    ir <- gaps(ir)
    if(length(ir)==0) return(defaultOut)
    ir <- ir[order(width(ir), decreasing=TRUE)]
    mid <- sum(c(start(ir)[1], end(ir)[1]))/2
    return(c(start=mid-scale, end=mid+scale))
}

orderedGR <- function(gr=GRanges()){
    if(length(gr)>0){
        gr[order(as.character(seqnames(gr)), start(gr))]
    }else{
        gr
    }
}

condenceGRs <- function(gr=GRanges(), FUN=sum){
    .gr <- reduce(gr, min.gapwidth=0, with.revmap=TRUE)
    scores <- score(gr)
    .gr$score <- sapply(.gr$revmap, function(.id) FUN(scores[.id]))
    .gr$revmap <- NULL
    .gr
}

disjoinGRs <- function(gr=GRanges(), FUN=sum){
    if(length(gr)<1) return(gr)
    .gr <- disjoin(gr)
    ol <- findOverlaps(.gr, gr)
    s <- tapply(score(gr[subjectHits(ol)]), queryHits(ol), FUN=FUN)
    .gr$score <- 0
    .gr$score[as.numeric(names(s))] <- s
    .gr
}

filterTracks <- function(tl, chrom, from, to, st){
    for(i in seq_along(tl)){
        if(tl[[i]]@type %in% c("data", "scSeq")){
            if(tl[[i]]@format=="WIG") {
                tl[[i]] <- parseWIG(tl[[i]], chrom, from, to)
            }
            dat <- tl[[i]]@dat
            dat <- disjoinGRs(dat)
            tl[[i]]@dat <- dat[end(dat)>=from &
                                   start(dat)<=to &
                                   seqnames(dat)==chrom]
            if(length(tl[[i]]@dat2)>0){
                dat2 <- tl[[i]]@dat2
                dat2 <- disjoinGRs(dat2)
                tl[[i]]@dat2 <- dat2[end(dat2)>=from &
                                        start(dat2)<=to &
                                        seqnames(dat2)==chrom]
            }
            if(st %in% c("+", "-")) {
                dat <- tl[[i]]@dat
                tl[[i]]@dat <- dat[strand(dat)==st]
                if(length(tl[[i]]@dat2)>0){
                    dat2 <- tl[[i]]@dat2
                    tl[[i]]@dat2 <- dat2[strand(dat2)==st]
                }
            }
        }else{
          if(tl[[i]]@type=="interactionData"){
            ## dat, dat2 are paired, or with target
            getKeep <- function(dat, dat2){
              keep <- ((end(dat)>=from & start(dat)<=to) |
                         (end(dat2)>=from & start(dat2)<=to)) & 
                seqnames(dat)==chrom & seqnames(dat2)==chrom
              ## remove duplicates
              idx1 <- paste(as.character(seqnames(dat)), start(dat), end(dat),
                            as.character(seqnames(dat2)), start(dat2), end(dat2))
              idx2 <- paste(as.character(seqnames(dat2)), start(dat2), end(dat2),
                            as.character(seqnames(dat)), start(dat), end(dat))
              idx <- ifelse(start(dat)<start(dat2), idx1, idx2)
              keep <- keep & (!duplicated(idx))
            }
            .dat <- tl[[i]]@dat
            .dat2 <- tl[[i]]@dat2
            .dat1target <- NULL
            .dat2target <- NULL
            if(length(.dat$target)==length(.dat) && is(.dat$target, "GRanges")){
              .dat1target <- .dat$target
              if(length(.dat2)!=0){
                if(length(.dat2$target)==length(.dat2) && is(.dat2$target, "GRanges")){
                  names(.dat2) <- NULL
                  .dat2target <- .dat2$target
                }
              }
            }else{
              .dat1target <- .dat2
            }
            keep2 <- keep1 <- getKeep(.dat, .dat1target)
            tl[[i]]@dat <- tl[[i]]@dat[keep1]
            if(length(.dat2target)){
              keep2 <- getKeep(.dat2, .dat2target)
            }
            tl[[i]]@dat2 <- tl[[i]]@dat2[keep2]
          }else{
            if(tl[[i]]@type=="lollipopData"){
              dat <- tl[[i]]@dat
            }else{
              dat <- range(unname(tl[[i]]@dat))
            }
            dat <- dat[end(dat)>=from &
                         start(dat)<=to &
                         seqnames(dat)==chrom]
            dat2 <- tl[[i]]@dat2
            if(length(dat2)>0){
              dat2 <- dat2[end(dat2)>=from &
                             start(dat2)<=to &
                             seqnames(dat2)==chrom]
            }
            if(tl[[i]]@type=="lollipopData"){
              tl[[i]]@dat <- dat
              tl[[i]]@dat2 <- dat2
            }
            if(length(dat)==0 && length(dat2)==0)
              tl[[i]]@style@height <- 0
          }
        }
    }
    tl
}

is_null_na <- function(.ele){
  if(is.null(.ele)) return(TRUE)
  is.na(.ele)
}
getYlim <- function(tl, op){
    yscales <- mapply(tl, op, FUN=function(.ele, .op){
        ylim <- .ele@style@ylim
        if(length(ylim)!=2){
            if(.ele@type %in% c("data", "lollipopData", "scSeq")){
                if(length(.ele@dat)>0){
                    ylim <- unique(round(range(.ele@dat$score)))
                }else{
                    ylim <- c(0, 0)
                }
                if(length(.ele@dat2)>0 && is_null_na(.op)[1]){
                    ylim2 <- unique(round(range(.ele@dat2$score)))
                    ylim <- c(ylim, -1*ylim2)
                }
                ylim <- range(c(0, ylim))
            }else{
              if(.ele@type == "interactionData"){
                ## max interaction height
                ylim <- c(0, 1)
              }else{
                ylim <- c(0, 0)
              }
            }            
        }
        ylim
    }, SIMPLIFY = FALSE)
    
    yscaleR <- range(unlist(yscales))
    if(diff(yscaleR)==0) yscaleR <- c(0, 1)
    yscales <- lapply(yscales, function(.ele){
        if(diff(.ele)==0){
            if(all(.ele>=0)){
                .ele <- c(0, yscaleR[2])
            }else{
                .ele <- c(yscaleR[1], 0)
            }
        }
        if(.ele[1]>.ele[2]) .ele[1] <- 0
        .ele
    })
    names(yscales) <- names(tl)
    yscales
}

getYheight <- function(tl){
    yHeights <- sapply(tl, function(.ele){
        yh <- .ele@style@height
        if(length(yh)==0) yh <- -1
        yh[1]
    })
    noY <- yHeights == -1
    yHeightsT <- sum(yHeights[!noY])
    if(yHeightsT>1.001)
        stop("total heights of data tracks is greater than 1.")
    if(length(yHeights[noY]) > 0){
        yHeights[noY] <- 
            (1 - yHeightsT) / length(yHeights[noY])
    }
    names(yHeights) <- names(tl)
    yHeights
}

drawXaxis <- function(xscale, style){
    scale <- hrScale(xscale)
    suffix <- c(1, 1000, 1000000, 1000000000)
    names(suffix) <- c("bp", "K", "M", "G")
    interval <- scale$scale * suffix[scale$unit]
    start <- ceiling(xscale[1]/interval)
    end <- floor(xscale[2]/interval)
    if(style@flip) xscale <- rev(xscale)
    if(length(style@xat)>0){
      at <- style@xat
      if(length(style@xlabel)==length(style@xat)){
        label <- style@xlabel
      }else{
        label <- at
      }
      at <- rescale(at, from=xscale)
    }else{
      label <- interval * start:end
      at <- rescale(label, from=xscale)
    }
    gp <- style@xgp
    class(gp) <- "gpar"
    rot <- ifelse(style@xlas %in% c(0, 1), 0, 90)
    grid.xaxis(at=at, label=label, gp=gp,
               edits = gEdit(gPath="labels", rot=rot),
               draw=style@xaxis)
}
putGeneYlab <- function(curViewStyle, style, name, height, xscale, rang, withlollipop=FALSE){
    gap <- (xscale[2] - xscale[1])/100
    strand <- unique(as.character(strand(rang)))
    if(length(strand)>1){
      strand <- names(sort(table(strand), decreasing = TRUE))[1]
      strand(rang) <- strand
      rang <- range(rang)
    }
    just <- style@ylabpos=="upstream"
    if(curViewStyle@flip){
        if(strand=="+"){
            x <- ifelse(just, start(rang) + gap, end(rang) - gap)
            just <- ifelse(just, "left", "right")
        }else{
            x <- ifelse(just, end(rang) - gap, start(rang) + gap)
            just <- ifelse(just, "right", "left")
        }
    }else{
        if(strand=="+"){
            x <- ifelse(just, start(rang) - gap, end(rang) + gap)
            just <- ifelse(just, "right", "left")
        }else{
            x <- ifelse(just, end(rang) + gap, start(rang) - gap)
            just <- ifelse(just, "left", "right")
        }
    }
    
    gp <- style@ylabgp
    class(gp) <- "gpar"
    if(is.null(gp$cex)) gp$cex <- optFontSize("z", curViewStyle, 
                                              height=height)
    pushViewport(viewport(x=curViewStyle@margin[2], y=0, 
                          height=1, 
                          width=1-curViewStyle@margin[2]-curViewStyle@margin[4], 
                          clip="off",
                          just=c(0,0), 
                          xscale=xscale))
    grid.text(x=x, y=ifelse(withlollipop, .25, .5), label=name, rot=0, just=just, gp=gp, 
              default.units="native")
    popViewport()
}

putYlab <- function(curViewStyle, style, name, yHeightBottom, yHeightTop, height, yscale){
    ##c("left", "right", "topleft", "bottomleft", "topright", "bottomright", "abovebaseline", "underbaseline")
    vp <- switch(style@ylabpos,
                 left=viewport(x=curViewStyle@margin[2]*.4, 
                               width=curViewStyle@margin[2]*.8,
                               just="center"),
                 right=viewport(x=1 - curViewStyle@margin[4]*.6, 
                                width=curViewStyle@margin[4]*.8,
                                just="center"),
                 topleft=viewport(y=1-yHeightTop, height=yHeightTop, just="bottom"),
                 bottomleft=viewport(y=0, height=yHeightBottom, just="bottom"),
                 topright=viewport(y=1-yHeightTop, height=yHeightTop, just="bottom"),
                 bottomright=viewport(y=0, height=yHeightBottom, just="bottom"),
                 abovebaseline=viewport(y=0.5, height = 1, just="center",yscale=yscale),
                 underbaseline=viewport(y=0.5, height = 1, just="center",yscale=yscale),
                 viewport(x=curViewStyle@margin[2]*.4, 
                          width=curViewStyle@margin[2]*.8))
    just <- switch(style@ylabpos,
                   left="center",
                   right="center",
                   topleft=c(0, 0),
                   topright=c(1, 0),
                   bottomleft=c(0, 1),
                   bottomright=c(1, 1),
                   abovebaseline=c(0, 0),
                   underbaseline=c(0, 1),
                   "center"
        )
    x <- switch(style@ylabpos,
                left=.5,
                right=.5,
                topleft=curViewStyle@margin[2],
                topright=1-curViewStyle@margin[4],
                bottomleft=curViewStyle@margin[2],
                bottomright=1-curViewStyle@margin[4],
                abovebaseline=curViewStyle@margin[2],
                underbaseline=curViewStyle@margin[2],
                .5
        )
    y <- switch(style@ylabpos,
                left=.5,
                right=.5,
                topleft=.1,
                topright=.1,
                bottomleft=.9,
                bottomright=.9,
                abovebaseline=unit(0, "native")+unit(2, "points"),
                underbaseline=unit(0, "native")-unit(2, "points"),
                .5
        )
    pushViewport(vp)
    rot <- ifelse(style@ylablas %in% c(0, 3), 90, 0)
    gp <- style@ylabgp
    class(gp) <- "gpar"
    if(style@ylabpos %in% c("topleft", "topright", "bottomleft", "bottomright", "abovebaseline", "underbaseline")){
        rot <- 0
        curHeight <- ifelse(grepl("top",style@ylabpos), yHeightTop, yHeightBottom)
        if(is.null(gp$cex)) gp$cex <- optFontSize("z", curViewStyle, 
                                                  height=curHeight)
    }else{
        if(curViewStyle@autolas){
            curWidth <- ifelse(style@ylabpos=="left", 
                               curViewStyle@margin[2],
                               curViewStyle@margin[4])
            rot <- ifelse(convertHeight(unit(height, "npc"), 
                                        unitTo="inches", valueOnly=TRUE)
                          > convertWidth(unit(curWidth, "npc"), 
                                         unitTo="inches", valueOnly=TRUE),
                          90, 0)
            if(is.null(gp$cex)) gp$cex <- optFontSize("z", curViewStyle, 
                                                      height=height)
        }
    }
    grid.text(x=x, y=y, label=name, rot=rot, just=just, gp=gp)

    popViewport()
}

drawYaxis <- function(ylim, yaxisStyle, curViewStyle, heatlegends=list()){
    if(yaxisStyle@main){
        vp <- viewport(x=curViewStyle@margin[2] * .2, 
                       width=curViewStyle@margin[2] *.4,
                       y=0, height=1, clip="off",
                       just=c(0,0), yscale=ylim)
    }else{
        vp <- viewport(x=curViewStyle@margin[2], 
                       width=1 - curViewStyle@margin[2] - curViewStyle@margin[4],
                       y=0, height=1, clip="off", just=c(0,0), yscale=ylim)
    }
    gp <- yaxisStyle@gp
    class(gp) <- "gpar"
    pushViewport(vp)
    if(length(heatlegends)){
      width <- ifelse(yaxisStyle@main, .3,
                      curViewStyle@margin[4]/4)
      x <- ifelse(yaxisStyle@main, 1.15, 1+curViewStyle@margin[4]/8)
      grid.raster(rev(heatlegends$crp), x=x, width=width, height=1)
      # map breaks to ylim
      if(heatlegends$userdefinedbreaks){
        at <- seq(ylim[1], ylim[2], length.out=length(heatlegends$breaks))
        label <- heatlegends$breaks
      }else{
        at <- ylim
        label <- range(heatlegends$breaks)
      }
      tryCatch({
        dif <- rle(diff(label))$lengths
        w <- c(1, 1+cumsum(dif))
        at <- at[w]
        label <- label[w]
      }, error=function(.e){}, warning=function(.w){})
      tryCatch({
        if(any(label%%1!=0)){
          label <- formatC(label, format = "fg")
        }
      }, error=function(.e){}, warning=function(.w){})
      grid.yaxis(at=at, label = label,
                 main=FALSE, gp=gp,
                 draw=yaxisStyle@draw)
    }else{
      at <- unique(ylim)
      at <- at[order(at)]
      label <- yaxisStyle@label
      if(label){
        label <- at
      }
      grid.yaxis(at=at, label=label, main=FALSE, gp=gp,
                 draw=yaxisStyle@draw)
    }
    popViewport()
}

drawXscale <- function(scale){
    gp <- scale@gp
    class(gp) <- "gpar"
    line1 <- abs(as.numeric(convertY(unit(1, "line"), 
                                 scale@from@unit)))
    if(length(gp$cex)>0){
        if(is.numeric(gp$cex)){
            line1 <- line1 * gp$cex
        }
    }
    if(scale@from@y<0){
        scale@from@y <- scale@from@y - line1
        scale@to@y <- scale@to@y - line1
    }
    grid.segments(x0=scale@from@x, 
                  y0=scale@from@y,
                  x1=scale@to@x,
                  y1=scale@to@y,
                  default.units=scale@from@unit,
                  gp=gp)
    grid.segments(x0=scale@from@x,
                  y0=scale@from@y,
                  x1=scale@from@x,
                  y1=scale@from@y + 0.25 * line1,
                  default.units=scale@from@unit,
                  gp=gp)
    grid.segments(x0=scale@to@x,
                  y0=scale@to@y,
                  x1=scale@to@x,
                  y1=scale@to@y + 0.25 * line1,
                  default.units=scale@to@unit,
                  gp=gp)
    grid.text(label=scale@label, 
              x=(scale@from@x + scale@to@x)/2,
              y=(scale@from@y + scale@to@y)/2 + line1 * 0.2,
              gp=gp,
              just="bottom",
              default.units=scale@from@unit)
}

maxStringWidth <- function(labels, spaces="WW", cex){
    max(as.numeric(convertX(stringWidth(paste0(labels, spaces)), "line"))*cex)
}

getColNum <- function(labels, spaces="WW", cex){
    ncol <- floor(as.numeric(convertX(unit(1, "npc"), "line")) / 
                      maxStringWidth(labels, spaces=spaces, cex) / 
              as.numeric(convertX(stringWidth("W"), "line")))
    nrow <- ceiling(length(labels) / ncol)
    ncol <- ceiling(length(labels) / nrow)
    ncol
}

checkLegendPosition <- function(legendPosition){
  positions <- c('top', 'left', 'right')
  if(!is.list(legendPosition)){
    legendPosition <- list(
      position = match.arg(legendPosition, choices = positions)
    )
  }else{
    legendPosition$position <- 
      match.arg(legendPosition$position,
                choices = positions)
  }
  return(legendPosition)
}
############### handle legend ####################
## set the legend as a list, 
## if all the legend for different tracks is same
## set draw legend for last track later
handleLegend <- function(legend, len, dat){
  if(length(legend)>0){
    if(is.character(legend)){
      if(!missing(dat)){
        if(is.list(dat)){
          if(length(legend)<length(dat)){
            legend <- rep(legend, length(dat))[seq_along(dat)]
          }
        }else{
          dat <- list(dat)
          legend <- legend[1]
        }
        para <- c("shape", "color", "border", "alpha")
        preset <- list("circle", "white", "black", 1)
        shapeMap <- c("circle"=21, "square"=22, "diamond"=23, 
                      "triangle_point_up"=24, "triangle_point_down"=25)
        names(preset) <- para
        legend <- mapply(function(.legend, .dat){
          coln <- colnames(mcols(.dat))
          if(.legend %in% coln){
            labels <- mcols(.dat)[, .legend]
            gp <- lapply(para, function(.ele){
              if(.ele %in% coln){
                mcols(.dat)[, .ele]
              }else{
                rep(preset[[.ele]], length(.dat))
              }
            })
            names(gp) <- para
            names(gp)[names(gp)=="color"] <- "fill"
            if(is.list(gp[["shape"]])){
              warning('multiple shape for one legend is not accept. Using default')
              gp[["shape"]] <- shapeMap[vapply(gp[["shape"]], FUN=function(.ele) .ele[1], FUN.VALUE = character(1L))]
            }else{
              gp[["shape"]] <- shapeMap[gp[["shape"]]]
            }
            gp <- as.data.frame(gp, stringsAsFactors=FALSE)
            gp <- cbind(labels=labels, gp)
            names(gp)[names(gp)=="shape"] <- "pch"
            gp <- gp[!duplicated(gp[, "labels"]), ]
            gp <- gp[order(gp[, "labels"]), ]
            gp <- as.list(gp)
          }
        }, legend, dat, SIMPLIFY = FALSE)
      }
    }
    if(!is.list(legend)){
      tmp <- legend
      legend <- vector(mode = "list", length = len)
      legend[[len]] <- tmp
      rm(tmp)
    }else{
      if(length(legend)==1){
        tmp <- legend[[1]]
        legend <- vector(mode = "list", length = len)
        legend[[len]] <- tmp
        rm(tmp)
      }else{
        if("labels" %in% names(legend)){
          tmp <- legend
          legend <- vector(mode = "list", length = len)
          legend[[len]] <- tmp
          rm(tmp)
        }else{
          if(length(legend)<len){
            length(legend) <- len
          }
        }
      }
    }
  }
  return(legend)
}

legendNchar <- function(legend){
  n <- 0
  if(length(legend)>0){
    if(is.list(legend)){
      thisLabels <- legend[["labels"]]
    }else{
      thisLabels <- names(legend)
    }
    if(length(thisLabels)>0){
      n <- nchar(thisLabels)
    }
  }
  return(n)
}

getMaxYlimNchar <- function(SNP.gr, minVal, types){
  if(!is.list(SNP.gr)){
    SNP.gr <- list(SNP.gr)
  }
  m <- lapply(SNP.gr[types!='pie'], function(.ele){
    .e <- .ele$score
    if(length(.e)>0){
      return(max(.e[!is.infinite(.e)], na.rm=TRUE))
    }else{
      return(0)
    }
  })
  m <- max(c(1, unlist(m, recursive = TRUE)))
  m <- nchar(as.character(round(m))) + 3.5
  if(m<minVal){
    return(minVal)
  }else{
    return(min(10, m))
  }
}
################ handle ranges #####################
## if !missing(ranges) set ranges as feature ranges
handleRanges <- function(ranges, SNP.gr, features, len){
  if(length(ranges)>0){
    stopifnot(inherits(ranges, c("GRanges", "GRangesList", "list")))
    if(is(ranges, "GRanges")){
      if(length(ranges)==1){
        ranges <- split(rep(ranges, len)[seq.int(len)],
                        seq.int(len))
      }else{
        ranges <- split(rep(ranges, len),
                        rep(seq.int(len), each=len))[seq.int(len)]
      }
    }else{## GRangesList
      if(length(ranges)!=len){
        ranges <- rep(ranges, seq.int(len))[seq.int(len)]
      }
    }
    stopifnot(length(ranges)==len)
  }else{
    if(is(features, "GRanges")){
      ranges <- split(range(unname(features), ignore.strand=TRUE)[rep(1, len)],
                      seq.int(len))
    }else{
      if(length(features)!=len){
        stop("if both SNP.gr and features is GRangesList,",
             " the lengthes of them should be identical.")
      }
      ranges <- GRangesList(lapply(features, function(.ele){
        range(unname(.ele), ignore.strand=TRUE)}))
    }
  }
  return(ranges)
}

##cut all SNP.gr by the range
cutSNP <- function(SNP.gr, ranges, len){
  if(is(ranges, "GRanges")){
    for(i in seq.int(len)){
      range <- ranges[i]
      stopifnot(all(width(SNP.gr[[i]])==1))
      SNP.gr[[i]] <- subsetByOverlaps(SNP.gr[[i]], range, ignore.strand=FALSE)
    }
  }else{
    if(inherits(ranges, c("GRangesList", "list"))){
      for(i in seq.int(len)){
        range <- ranges[[i]]
        stopifnot(all(width(SNP.gr[[i]])==1))
        SNP.gr[[i]] <- subsetByOverlaps(SNP.gr[[i]], range, ignore.strand=FALSE)
      }
    }
  }
  return(SNP.gr)
}

## multiple transcripts in one gene could be separated by featureLayerID
setFeatureLayerID <- function(feature, range){
  feature <- feature[end(feature)>=start(range) & 
                       start(feature)<=end(range)]
  if(length(feature$featureLayerID)!=length(feature)){
    feature$featureLayerID <- rep("1", length(feature))
    feature$featureLayerID <- as.character(feature$featureLayerID)
    start(feature)[start(feature)<start(range)] <- start(range)
    end(feature)[end(feature)>end(range)] <- end(range)
  }
  return(feature)
}

## bottomblank, the transcripts legend height
plotFeatureLegend <- function(feature, LINEH, range, xaxis, xaxis.gp, label_on_feature=FALSE){
  if(length(xaxis)>1 || as.logical(xaxis[1])){
    xaxisSpace <- 2
    if(is.numeric(xaxis.gp$cex)) xaxisSpace <- 2*xaxis.gp$cex
  }else{
    xaxisSpace <- 0
  }
  if(length(names(feature))>0 & !label_on_feature ){ ## features legend
    feature.s <- feature[!duplicated(names(feature))]
    cex <- if(length(unlist(feature.s$cex))==length(feature.s)) 
      unlist(feature.s$cex) else 1
    ncol <- getColNum(names(feature.s), cex=cex)
    featureLegendSpace <- max(ceiling(length(names(feature.s)) / ncol) * cex + 1 )
    pushViewport(viewport(x=.5, y=featureLegendSpace*LINEH/2, 
                          width=1,
                          height=featureLegendSpace*LINEH,
                          xscale=c(start(range), end(range))))
    color <- if(length(unlist(feature.s$color))==length(feature.s)) 
      unlist(feature.s$color) else "black"
    fill <- if(length(unlist(feature.s$fill))==length(feature.s)) 
      unlist(feature.s$fill) else "black"
    pch <- if(length(unlist(feature.s$pch))==length(feature.s)) 
      unlist(feature.s$pch) else 22
    grid.legend(label=names(feature.s), ncol=ncol,
                byrow=TRUE, vgap=unit(.2, "lines"),
                hgap=unit(.5, "lines"),
                pch=pch,
                gp=gpar(col=color, fill=fill, cex=cex))
    popViewport()
  }else{
    featureLegendSpace <- 1
  }
  bottomblank <- (xaxisSpace + featureLegendSpace) * LINEH
  return(bottomblank)
}

plot_grid_xaxis <- function(xaxis, gp=gpar(col="black")){
  ## axis, should be in the bottom of transcripts
  if(length(xaxis)==1 && as.logical(xaxis)) {
    grid.xaxis(gp=gp)
  }
  if(length(xaxis)>1 && is.numeric(xaxis)){
    xaxisLabel <- names(xaxis)
    if(length(xaxisLabel)!=length(xaxis)) xaxisLabel <- TRUE
    grid.xaxis(at=xaxis, label=xaxisLabel, gp=gp)
  }
}

getYXratio <- function(){
  as.numeric(convertHeight(unit(1, 'snpc'), 'npc'))/
    as.numeric(convertWidth(unit(1, "snpc"), "npc"))
}
#"circle", "square", "diamond", "triangle_point_up", "star", or "triangle point_down"
grid.circle1 <- function(x = 0.5, y = 0.5, r = 0.5, 
                         default.units = "npc", name = NULL, 
                         gp = gpar(), draw = TRUE, vp = NULL){
  fill <- gp$fill
  col <- gp$col
  lwd <- if(length(gp$lwd)>0) gp$lwd else 1
  alpha <- gp$alpha
  if(is.null(fill)) fill <- "white"
  twopi <- 2 * pi
  ratio.yx <- getYXratio()
  t2xy <- function(t) {
    t2p <- twopi * t + pi/2
    list(x = r * cos(t2p)/ratio.yx, y = r * sin(t2p))
  }
  P <- t2xy(seq.int(0, 1, length.out = 100))
  invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                         gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
}

grid.square <- function(x = 0.5, y = 0.5, r = 0.5, 
                        default.units = "npc", name = NULL, 
                        gp = gpar(), draw = TRUE, vp = NULL){
  fill <- gp$fill
  col <- gp$col
  lwd <- if(length(gp$lwd)>0) gp$lwd else 1
  alpha <- gp$alpha
  if(is.null(fill)) fill <- "white"
  ratio.yx <- getYXratio()
  invisible(grid.rect(unit(x,"npc"), unit(y, "npc"), 
                      width = unit(r*2/ratio.yx, "npc"), 
                      height = unit(r*2, "npc"),
                      gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
}

grid.diamond <- function(x = 0.5, y = 0.5, r = 0.5, 
                        default.units = "npc", name = NULL, 
                        gp = gpar(), draw = TRUE, vp = NULL){
  fill <- gp$fill
  col <- gp$col
  lwd <- if(length(gp$lwd)>0) gp$lwd else 1
  alpha <- gp$alpha
  if(is.null(fill)) fill <- "white"
  ratio.yx <- getYXratio()
  P <- 
    list(x = c(0, r/ratio.yx, 0, -r/ratio.yx), 
         y = c(-r, 0, r, 0))
  invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                         gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
}

grid.triangle_point_up <- function(x = 0.5, y = 0.5, r = 0.5, 
                         default.units = "npc", name = NULL, 
                         gp = gpar(), draw = TRUE, vp = NULL){
  fill <- gp$fill
  col <- gp$col
  lwd <- if(length(gp$lwd)>0) gp$lwd else 1
  alpha <- gp$alpha
  if(is.null(fill)) fill <- "white"
  ratio.yx <- getYXratio()
  P <- 
    list(x = c(-r/ratio.yx, r/ratio.yx, 0, -r/ratio.yx), 
         y = c(-r, -r, r, -r))
  invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                         gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
}

grid.triangle_point_down <- function(x = 0.5, y = 0.5, r = 0.5, 
                                   default.units = "npc", name = NULL, 
                                   gp = gpar(), draw = TRUE, vp = NULL){
  fill <- gp$fill
  col <- gp$col
  lwd <- if(length(gp$lwd)>0) gp$lwd else 1
  alpha <- gp$alpha
  if(is.null(fill)) fill <- "white"
  ratio.yx <- getYXratio()
  P <- 
    list(x = c(-r/ratio.yx, r/ratio.yx, 0, -r/ratio.yx), 
         y = c(r, r, -r, r))
  invisible(grid.polygon(unit(P$x+x,"npc"), unit(P$y+y, "npc"), 
                         gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
}


grid.star <- function(x = 0.5, y = 0.5, r = 0.5, 
                      default.units = "npc", name = NULL, 
                      gp = gpar(), draw = TRUE, vp = NULL){
  fill <- gp$fill
  col <- gp$col
  lwd <- if(length(gp$lwd)>0) gp$lwd else 1
  alpha <- gp$alpha
  if(is.null(fill)) fill <- "white"
  ratio.yx <- getYXratio()
  i <- 1:11
  angle <- 180
  alpha <- 2*pi / 10
  r <- r * (i %% 2 + 1)/2
  omega <- alpha * i + angle * pi /180
  invisible(grid.polygon(unit(r*sin(omega)/ratio.yx+x,"npc"), 
                         unit(r*cos(omega)+y, "npc"), 
                         gp=gpar(col = col, fill = fill, lwd=lwd, alpha=alpha)))
}

cleanDataMcols <- function(this.dat.mcols, type){
  this.dat.mcols <- 
    this.dat.mcols[, 
                   !colnames(this.dat.mcols) %in% 
                     c("color", "fill", "lwd", "id", 
                       "cex", "dashline.col", 
                       "id.col", "stack.factor", "SNPsideID",
                       "shape", "alpha"), 
                   drop=FALSE]
  if(type!="pie.stack"){
    this.dat.mcols <- 
      this.dat.mcols[, !colnames(this.dat.mcols) %in% 
                       c("stack.factor.order", 
                         "stack.factor.first"), 
                     drop=FALSE]
  }
  this.dat.mcols <- 
    this.dat.mcols[, !grepl("^label.parameter",
                            colnames(this.dat.mcols)), 
                   drop=FALSE]
  this.dat.mcols <- 
    this.dat.mcols[, !grepl("^node.label",
                            colnames(this.dat.mcols)), 
                   drop=FALSE]
  return(this.dat.mcols)
}
