#' Class \code{"trackViewerStyle"}
#' @description An object of class \code{"trackViewerStyle"} 
#'              represents track viewer style.
#' @aliases trackViewerStyle
#' @rdname trackViewerStyle-class
#' @slot margin \code{"numeric"}, specify the bottom, left, top and right margin.
#' @slot xlas \code{"numeric"}, label direction of x-axis mark. It should 
#' be a integer 0-3. See \code{\link[graphics]{par}:las}
#' @slot xgp A \code{"list"}, object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of x-axis. For y-axis, see \code{\link{yaxisStyle}}
#' @slot xaxis \code{"logical"}, draw x-axis or not
#' @slot xat \code{"numeric"}, the values will be passed to grid.xaxis as
#' 'at' parameter.
#' @slot xlabel \code{"character"}, the values will be passed to grid.xaxis as
#' 'label' parameter.
#' @slot autolas \code{"logical"} automatic determine y label direction
#' @slot flip \code{"logical"} flip the x-axis or not, default FALSE
#' @import methods
#' @exportClass trackViewerStyle
#' @examples
#' tvs <- trackViewerStyle()
#' setTrackViewerStyleParam(tvs, "xaxis", TRUE)
#' 
setClass("trackViewerStyle", 
         representation(
             margin="numeric",
             xlas="numeric",
             xgp="list",
             xaxis="logical",
             xat="numeric",
             xlabel="character",
             autolas="logical",
             flip="logical"),
         prototype(
             margin=c(.01, .05, 0, 0),
             xlas=0,
             xgp=list(),
             xaxis=FALSE,
             xat=numeric(0L),
             xlabel=character(0L),
             autolas=FALSE,
             flip=FALSE
             ),
         validity=function(object){
             if(!(object@xlas %in% 0:3))
                 return("xlas should be numeric in {0,1,2,3}. See ?par")
             if(any(object@margin<0 | object@margin>.8))
                 return("margin could not greater than .8 or smaller than 0")
             return(TRUE)
         }
)

#' @rdname trackViewerStyle-class
#' @param \dots Each argument in \dots becomes an slot in the new trackViewerStyle.
#' @export

trackViewerStyle <- function(...){
    new("trackViewerStyle", ...)
}

#' @rdname trackViewerStyle-class
#' @param tvs An object of \code{trackViewerStyle}.
#' @param attr the name of slot to be changed.
#' @param value values to be assigned.
#' @exportMethod setTrackViewerStyleParam
#' @aliases setTrackViewerStyleParam
#' @aliases setTrackViewerStyleParam,trackViewerStyle,character-method
#' 
setGeneric("setTrackViewerStyleParam", function(tvs, attr, value) 
    standardGeneric("setTrackViewerStyleParam"))
#' @rdname trackViewerStyle-class
#' @aliases setTrackViewerStyleParam,trackViewerStyle,character,ANY-method
#' @details 
#' \code{setTrackViewerStyleParam} changes one slot of a
#' \code{\link{trackViewerStyle}} object. Unlike \code{setTrackStyleParam} and
#' its siblings, this style is not attached to a single track: it controls the
#' figure-wide x-axis shared by every track in a plot (drawn once, beneath the
#' bottom-most track), so this function is typically called on the
#' \code{trackViewerStyle} object passed to \code{viewTracks}, not on a track
#' itself.
#'
#' Accepted values of \code{attr}, all slots of \code{\link{trackViewerStyle}}:
#' \describe{
#'   \item{\code{margin}}{\code{"numeric"} of length 4, the bottom, left, top
#'   and right margins of the whole plot, as fractions of the device (each
#'   value must be between 0 and 0.8). Default \code{c(.01, .05, 0, 0)}.}
#'   \item{\code{xlas}}{\code{"numeric"} in \{0, 1, 2, 3\}, the rotation of the
#'   x-axis tick labels. See \code{\link[graphics]{par}:las}. Ignored when
#'   \code{autolas} is \code{TRUE}.}
#'   \item{\code{xgp}}{\code{"list"} of graphical parameters for the x-axis,
#'   converted to \code{\link[grid]{gpar}}, e.g. \code{list(cex=.8, col="black")}.
#'   For the equivalent y-axis settings, see \code{\link{yaxisStyle}} and
#'   \code{\link{setTrackYaxisParam}}, which are per-track rather than
#'   figure-wide.}
#'   \item{\code{xaxis}}{\code{"logical"}, whether the shared x-axis is drawn at
#'   all. Default \code{FALSE}.}
#'   \item{\code{xat}}{\code{"numeric"} vector of genomic positions at which
#'   x-axis tick marks are drawn, passed to \code{grid.xaxis} as \code{at}. If
#'   left at its default (\code{numeric(0)}), positions are chosen
#'   automatically from the plotted range.}
#'   \item{\code{xlabel}}{\code{"character"} vector of labels for the ticks in
#'   \code{xat}, passed to \code{grid.xaxis} as \code{label}. Should be the
#'   same length as \code{xat} when both are set explicitly; left at its
#'   default, labels are generated automatically from the tick positions.}
#'   \item{\code{autolas}}{\code{"logical"}, whether the y-axis label direction
#'   is chosen automatically rather than following \code{xlas}. Default
#'   \code{FALSE}.}
#'   \item{\code{flip}}{\code{"logical"}, whether the x-axis (and the whole
#'   plot) is flipped, e.g. to display a feature on the minus strand
#'   5' to 3' left to right. Default \code{FALSE}.}
#' }
#'
#' Any other value of \code{attr} is rejected. Values are assigned with
#' \code{check = TRUE}, so the validity method of
#' \code{\link{trackViewerStyle}} runs on every call — an \code{xlas} outside
#' 0:3 or a \code{margin} outside [0, 0.8] fails immediately with a descriptive
#' message rather than surfacing later as a plotting error.
#' @section Assignment in the calling frame:
#' Like \code{\link{setTrackStyleParam}}, this method assigns the updated
#' object back to the variable passed as \code{tvs} in the caller's
#' environment, so \code{setTrackViewerStyleParam(tvs, "xaxis", TRUE)} is
#' sufficient on its own and the result need not be re-assigned. It is also
#' returned invisibly for use in contexts, such as \code{lapply}, where that
#' calling-frame assignment has no visible target — see
#' \code{\link{setTrackStyleParam}} for a worked example of that case.
#' @return An object of class \code{\link{trackViewerStyle}} with the modified
#' slot, returned invisibly. The variable supplied as \code{tvs} is updated as
#' a side effect.
#' @examples
#' tvs <- trackViewerStyle()
#'
#' ## turn on the shared x-axis and rotate its labels
#' setTrackViewerStyleParam(tvs, "xaxis", TRUE)
#' setTrackViewerStyleParam(tvs, "xlas", 2)
#'
#' ## choose explicit tick positions and labels
#' setTrackViewerStyleParam(tvs, "xat", c(122929000, 122929500, 122930000))
#' setTrackViewerStyleParam(tvs, "xlabel", c("122,929,000", "122,929,500", "122,930,000"))
#'
#' ## tighten the plot margins
#' setTrackViewerStyleParam(tvs, "margin", c(.02, .08, .01, .01))
#'
#' ## flip the plot, e.g. for a minus-strand gene
#' setTrackViewerStyleParam(tvs, "flip", TRUE)
#'
#' \dontrun{
#' ## rejected: xlas must be one of 0:3
#' setTrackViewerStyleParam(tvs, "xlas", 5)
#' ## rejected: not a slot of trackViewerStyle
#' setTrackViewerStyleParam(tvs, "color", "red")
#' }
#' @seealso \code{\link{trackViewerStyle}}, \code{\link{setTrackStyleParam}},
#' \code{\link{setTrackXscaleParam}}, \code{\link{setTrackYaxisParam}}
setMethod("setTrackViewerStyleParam", 
          signature(tvs="trackViewerStyle", attr="character", value="ANY"),
          function(tvs, attr, value){
              if(!attr %in% c("margin", "xlas", "xgp", "xat", "xlabel", "xaxis", "autolas", "flip"))
                  stop("attr must be a slot name of trackViewerStyle")
              x <- tvs
              slot(x, attr, check = TRUE) <- value
              eval.parent(substitute(tvs <- x))
              return(invisible(x))
          })

#' Class \code{"pos"}
#' @description An object of class \code{"pos"} represents a point location
#' @rdname pos-class
#' @aliases pos
#' @slot x A \code{\link{numeric}} value, indicates the x position
#' @slot y A \code{\link{numeric}} value, indicates the y position
#' @slot unit \code{"character"} apecifying the units for the corresponding
#' numeric values. See \code{\link[grid]{unit}}
#' @exportClass pos
#' 
setClass("pos", representation(x="numeric", y="numeric", unit="character"),
         prototype(x=0.5, y=0.5, unit="npc"),
         validity=function(object){
             if(!object@unit %in% c("npc", "cm", "inches", "mm", "points", "picas",
                                    "bigpts", "dida", "cicero", "scaledpts",
                                    "lines", "char", "native", "snpc",
                                    "strwidth", "strheight", "grobwidth",
                                    "grobheight"))
                 return("unit must be units of grid unit. See ?unit")
             return(TRUE)
         })

#' Class \code{"xscale"}
#' @description An object of class \code{"xscale"} represents x-scale style.
#' @rdname xscale-class
#' @aliases xscale
#' @slot from A \code{\link{pos}} class, indicates the start point 
#' postion of x-scale.
#' @slot to A \code{\link{pos}} class, indicates the end point 
#' postion of x-scale.
#' @slot label \code{"character"} the label of x-scale
#' @slot gp A \code{"list"} object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of x-scale.
#' @slot draw A \code{"logical"} value indicating whether the x-scale
#' should be draw.
#' @exportClass xscale

setClass("xscale",
         representation(from="pos",
                        to="pos",
                        label="character",
                        gp="list",
                        draw="logical"),
         prototype(draw=FALSE, gp=list())
         )

#' Class \code{"yaxisStyle"}
#' @description An object of class \code{"yaxisStyle"} represents y-axis style.
#' @rdname yaxisStyle-class
#' @aliases yaxisStyle
#' @slot at \code{"numeric"} vector of y-value locations for the tick marks
#' @slot label \code{"logical"} value indicating whether to draw the 
#' labels on the tick marks.
#' @slot gp A \code{"list"} object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of y-axis.
#' @slot draw A \code{"logical"} value indicating whether the y-axis
#' should be draw.
#' @slot main A \code{"logical"} value indicating whether the y-axis
#' should be draw in left (TRUE) or right (FALSE).
#' @exportClass yaxisStyle

setClass("yaxisStyle",
         representation(at="numeric",
                        label="logical",
                        gp="list",
                        draw="logical",
                        main="logical"),
         prototype(draw=TRUE, label=FALSE, gp=list(), main=TRUE)
         )

#' Class \code{"trackStyle"}
#' @description An object of class \code{"trackStyle"} represents track style.
#' @rdname trackStyle-class
#' @aliases trackStyle
#' @slot tracktype \code{"character"} track type, could be peak, 
#' line, histogram, or annotation. 
#' Default is "peak". "annotation" is used to mark the peak regions.
#' For interaction data,
#' it could be "heatmap" or "link".
#' @slot color \code{"character"} track color. If the track has dat and dat2 slot,
#' it should have two values.
#' @slot NAcolor \code{"character"} NA color for interactionData.
#' @slot breaks \code{"numeric"} breaks for color keys of interactionData.
#' @slot height \code{"numeric"} track height. It should be a value between 0 and 1
#' @slot marginTop \code{"numeric"} track top margin
#' @slot marginBottom \code{"numeric"} track bottom margin
#' @slot xscale object of \code{\link{xscale}}, describe the details of x-scale
#' @slot yaxis object of \code{\link{yaxisStyle}}, describe the details of y-axis
#' @slot ylim \code{"numeric"} y-axis range
#' @slot ylabpos \code{"character"}, ylable postion, ylabpos should 
#' be 'left', 'right', 'topleft', 'bottomleft', 'topright', 'bottomright',
#' 'abovebaseline', 'underbaseline', or 'none'.
#' For gene type track, it also could be 'upstream' or 'downstream'
#' @slot ylablas \code{"numeric"} y lable direction. It should 
#' be a integer 0-3. See \code{\link[graphics]{par}:las}
#' @slot ylabgp A \code{"list"} object, It will convert to an object of 
#' class \code{\link[grid]{gpar}}. This is basically a list of graphical 
#' parameter settings of y-label.
#' @slot ysplit A \code{"numeric"} to split y plot region for interaction data.
#' Default is 0.5, which will split the region half to half. It will only work
#' for back to back plot.
#' @exportClass trackStyle
#'

setClass("trackStyle",
         representation(tracktype="character",
                        color="character",
                        NAcolor="character",
                        breaks="numeric",
                        height="numeric",
                        marginTop="numeric",
                        marginBottom="numeric",
                        xscale="xscale",
                        yaxis="yaxisStyle",
                        ylim="numeric",
                        ylabpos="character",
                        ylablas="numeric",
                        ylabgp="list",
                        ysplit="numeric"
                        ),
         prototype(
             marginTop=0,
             marginBottom=0.05,
             color=c("black","black"),
             NAcolor="white",
             breaks=0,
             tracktype="peak",
             ylabpos="left",
             ylablas=0,
             ylabgp=list(),
             ysplit=0.5
         ),
         validity=function(object){
             if(!object@ylabpos %in% c("left", "right", "topleft", "bottomleft", 
                                       "topright", "bottomright", "upstream", "downstream",
                                       "abovebaseline", "underbaseline", "none"))
                 return("ylabpos should be 'left', 'right', 'topleft', 'bottomleft', 
                        'topright', 'bottomright', 'upstream', 'downstream', 
                        'abovebaseline', 'underbaseline' or 'none'.")
             if(!(object@ylablas %in% 0:3))
                return("ylas should be numeric in {0,1,2,3}. See ?par")
             if(!all(object@tracktype %in% c("peak", "annotation", "line", "histogram",
                                             "heatmap", "link")))
                 return("tracktype must be on of peak, annotation, line, histogram, heatmap or link")
             return(TRUE)
         }
)

#' Class \code{"track"}
#' @description An object of class \code{"track"} bundles the genomic data to be
#' plotted together with the style used to draw it. A \code{track} is the basic
#' unit accepted by \code{\link{viewTracks}}; several of them are collected into
#' a \code{\link{trackList}} for a multi-panel figure.
#' @rdname trackStyle-class
#' @aliases track
#' @slot dat Object of class \code{\link[GenomicRanges:GRanges-class]{GRanges}},
#' the primary data of the track. The metadata columns it must carry depend on
#' \code{type} (see \sQuote{Details}): coverage-like tracks
#' (\code{"data"}, \code{"scSeq"}), \code{"lollipopData"} and
#' \code{"interactionData"} require a \code{score} column, whereas
#' \code{"gene"} and \code{"transcript"} tracks require a \code{feature} column.
#' @slot dat2 Object of class \code{\link[GenomicRanges:GRanges-class]{GRanges}},
#' the optional second data set of the track; leave it as a zero-length
#' \code{GRanges} when unused. Its meaning depends on \code{type}:
#' \itemize{
#'   \item for \code{"data"}, \code{"scSeq"} and \code{"lollipopData"} tracks it
#'   is drawn back-to-back with \code{dat}: \code{dat} above the baseline and
#'   \code{dat2} below it as \code{-1 * score}, which is the usual way to
#'   compare two samples or the two strands of one sample in a single panel;
#'   \item for \code{"interactionData"} it holds the second anchor of each
#'   interaction, so it must either be the same length as \code{dat} (anchors
#'   matched element-wise) or carry a \code{target} metadata column.
#' }
#' When supplied it must contain a \code{score} column, and for lollipop data
#' every range must have width 1, exactly as for \code{dat}.
#' @slot type The type of track, one of:
#' \describe{
#'   \item{\code{"data"}}{continuous or interval scores — coverage, signal,
#'   peaks. Drawn as a peak, line or histogram depending on
#'   \code{style@@tracktype}.}
#'   \item{\code{"gene"}}{gene models. \code{dat} must have a \code{feature}
#'   column describing each range (for example \code{"exon"}, \code{"CDS"},
#'   \code{"utr5"}, \code{"utr3"}).}
#'   \item{\code{"transcript"}}{transcript models; same \code{feature}
#'   requirement as \code{"gene"}, but isoforms are drawn on separate rows.}
#'   \item{\code{"scSeq"}}{single-cell signal. Validated like \code{"data"},
#'   i.e. a \code{format} and a numeric \code{score} are required.}
#'   \item{\code{"lollipopData"}}{point features such as variants or
#'   modification sites. Every range in \code{dat} (and \code{dat2}) must have
#'   width 1; ranges wider than one base are rejected.}
#'   \item{\code{"interactionData"}}{chromatin interactions, for example from
#'   Hi-C or 4C. Drawn as a heatmap or as arcs according to
#'   \code{style@@tracktype}.}
#' }
#' @slot format The format the data was imported from: \code{"BED"},
#' \code{"bedGraph"}, \code{"WIG"}, \code{"BigWig"} or \code{"BAM"}. It is
#' required, and must be a single string, for \code{"data"} and \code{"scSeq"}
#' tracks; other types ignore it. For every format except \code{"WIG"} the
#' \code{score} column must be \code{numeric} or \code{integer}. \code{"WIG"} is
#' the exception: because a WIG file has no fixed step boundaries until it is
#' plotted, \code{score} is stored as a
#' \code{\link[IRanges:AtomicList-class]{CompressedCharacterList}} and is
#' converted on the fly — see \code{\link{importScore}}.
#' @slot style Object of class \code{\link{trackStyle}} controlling color,
#' height, margins, y-axis, x-scale and the drawing mode. Do not assign to its
#' slots directly; use \code{\link{setTrackStyleParam}},
#' \code{\link{setTrackXscaleParam}} and \code{\link{setTrackYaxisParam}}, which
#' validate the value and update the object in the calling frame.
#' @slot name unused yet
#' @details
#' The validity method enforces the type-specific requirements summarised above;
#' the most common errors are a missing \code{score} or \code{feature} metadata
#' column, a lollipop range wider than one base, and a \code{dat2} whose length
#' does not match \code{dat} for interaction data.
#'
#' Note that \code{dat} and \code{dat2} are plain \code{GRanges} objects, so the
#' seqlevel style of the track must match that of the other tracks and of the
#' plotting range. \code{seqlevelsStyle} and its replacement method are
#' defined for \code{track} and rename both \code{dat} and \code{dat2} at once.
#' @exportClass track
#' @examples
#' extdata <- system.file("extdata", package="trackViewer",
#' mustWork=TRUE)
#' fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED")
#' setTrackStyleParam(fox2, "color", c("red","green"))
#' setTrackXscaleParam(fox2, "gp", list(cex=.5))
#' setTrackYaxisParam(fox2, "gp", list(col="blue"))
#' fox2$dat <- GRanges(score=numeric(0))
#'
#' ## inspect and change the style
#' fox2$type
#' fox2$style$tracktype
#' setTrackStyleParam(fox2, "tracktype", "histogram")
#' setTrackStyleParam(fox2, "height", .2)
#'
#' @seealso Please try to use \code{\link{importScore}} and \code{\link{importBam}} to 
#' generate the object.
setClass("track", representation(dat="GRanges",
                                 dat2="GRanges",
                                 type="character",
                                 format="character",
                                 style="trackStyle",
                                 name="character"),
         validity=function(object){
             if(!object@type %in% 
                c("data", "gene", "transcript", "scSeq",
                  "lollipopData", "interactionData"))
                 return("type must be 'data', 'transcript', 'gene', 'scSeq',
                        'lollipopData', 'interactionData'")
             if(object@type %in% c("data", "scSeq")){
                 if(!length(object@format)==1){
                   return("format must be one of \"BED\", 
                            \"bedGraph\", \"WIG\", \"BigWig\"")
                 }
                 if(!object@format %in% 
                    c("BED", "bedGraph", "WIG", "BigWig", "BAM"))
                     return("format must be one of \"BED\", 
                            \"bedGraph\", \"WIG\", \"BigWig\"")
                 if(is.null(object@dat$score))
                     return("dat should contain score metadata.")
                 if(length(object@dat2)>0){
                     if(is.null(object@dat2$score))
                         return("dat2 should contain score metadata.")
                 }
                 if(object@format!="WIG"){
                     if(!inherits(object@dat$score, c("numeric", "integer")))
                         return("class of score metadata should be numeric")
                 }else{
                     if(!is(object@dat$score, "CompressedCharacterList"))
                         return("Please try ?imortScore for WIG files")
                 }
             }else{
               if(object@type=="lollipopData"){
                 if(is.null(object@dat$score))
                   return("dat should contain score metadata.")
                 if(!all(width(object@dat)==1)){
                   return("Width for lollipopData must be 1")
                 }
                 if(length(object@dat2)>0){
                   if(is.null(object@dat2$score))
                     return("dat2 should contain score metadata.")
                   if(!all(width(object@dat2)==1)){
                     return("Width for lollipopData must be 1")
                   }
                 }
               }else{
                 if(object@type=="interactionData"){
                   if(is.null(object@dat$score))
                     return("dat should contain score metadata.")
                   if(length(object@dat2)!=length(object@dat)){
                       if(length(object@dat$target)!=length(object@dat)){
                           return("dat2 should be same length of dat.") 
                       }else{
                           if(length(object@dat2)>0){
                               if(length(object@dat2$target)!=
                                  length(object@dat2)){
                                   return(paste("dat2 does not contain target",
                                                "metadata."))
                               }
                           }
                       }
                   }
                 }else{
                   if(is.null(mcols(object@dat)$feature))
                     return("The metadata of dat must contain colnumn 'feature'") 
                   if(length(object@dat2)>0){
                     if(is.null(object@dat2$score))
                       return("dat2 should contain score metadata.")
                     if(!all(width(object@dat2)==1)){
                       return("Width for lollipop data must be 1")
                     }
                   }
                 }
               }
             }
             return(TRUE)
         })

#' Method seqlevels
#' @rdname trackStyle-class
#' @exportMethod seqlevels
#' @aliases seqlevels,track-method
setMethod("seqlevels", "track", 
          function(x){ seqlevels(x@dat) })
#' Method seqlevelsStyle
#' @rdname trackStyle-class
#' @exportMethod seqlevelsStyle
#' @aliases seqlevelsStyle,track-method
setMethod("seqlevelsStyle", "track", 
          function(x){ seqlevelsStyle(x@dat) })
#' Method seqlevelsStyle<-
#' @rdname trackStyle-class
#' @exportMethod seqlevelsStyle<-
#' @aliases seqlevelsStyle<-,track-method
setReplaceMethod("seqlevelsStyle", "track", 
          function(x, value){ 
              seqlevelsStyle(x@dat) <- value
              seqlevelsStyle(x@dat2) <- value
              return(x)
})

#' @rdname trackStyle-class
#' @param object an object of trackStyle.
#' @exportMethod show
#' 
#' @aliases show,track-method
setMethod("show", "track", function(object){
    cat("This is an object of track\n", "slot name:", object@name, "\n", 
        "slot type:", object@type, "\n", "slot format:", object@format, "\n")
    cat("slot dat:\n")
    show(object@dat)
    cat("slot dat2:\n")
    show(object@dat2)
    cat("slot style: try object$style to see details.\n")
})
#' Method $
#' @rdname trackStyle-class
#' @param x an object of trackStyle or track
#' @param name slot name of trackStyle or track
#' @exportMethod $
#' @aliases $,track-method
#' @aliases $,trackStyle-method
setMethod("$", "track", function(x, name) slot(x, name))
setMethod("$", "trackStyle", function(x, name) slot(x, name))
#' Method $<-
#' @rdname trackStyle-class
#' @exportMethod $<-
#' @aliases $<-,track-method
#' @aliases $<-,trackStyle-method
setReplaceMethod("$", "track", 
                 function(x, name, value){
                     slot(x, name, check = TRUE) <- value
                     x
                 })
setReplaceMethod("$", "trackStyle", 
                 function(x, name, value){
                   slot(x, name, check = TRUE) <- value
                   x
                 })

#' @rdname trackStyle-class
#' @importFrom utils .DollarNames
#' @method .DollarNames track
#' @param pattern A regular expression. Only matching names are returned.
#' @export 
.DollarNames.track <- function(x, pattern=""){
  grep(pattern, slotNames(x), value = TRUE)
}
#' @rdname trackStyle-class
#' @importFrom utils .DollarNames
#' @method .DollarNames trackStyle
#' @export 
.DollarNames.trackStyle <- function(x, pattern=""){
  grep(pattern, slotNames(x), value = TRUE)
}

#' Method setTrackStyleParam
#' @rdname trackStyle-class
#' @param ts An object of \code{track}.
#' @param attr the name of slot of \code{\link{trackStyle}} object to be changed.
#' @param value values to be assigned.
#' @exportMethod setTrackStyleParam
#' @aliases setTrackStyleParam
#' @aliases setTrackStyleParam,track,character-method
setGeneric("setTrackStyleParam", function(ts, attr, value) 
    standardGeneric("setTrackStyleParam"))
#' 
#' @rdname trackStyle-class
#' @aliases setTrackStyleParam,track,character,ANY-method
#' @details
#' \code{setTrackStyleParam} changes one styling slot of the
#' \code{\link{trackStyle}} object held in \code{ts@@style}. It accepts only the
#' slots that are plain values; the two compound slots have their own setters,
#' \code{\link{setTrackXscaleParam}} for \code{xscale} and
#' \code{\link{setTrackYaxisParam}} for \code{yaxis}. Passing any other name,
#' including \code{"xscale"} or \code{"yaxis"}, raises an error.
#'
#' Accepted values of \code{attr}:
#' \describe{
#'   \item{\code{tracktype}}{\code{"character"}. How the data is drawn. For
#'   \code{"data"} tracks use \code{"peak"} (the default), \code{"line"},
#'   \code{"histogram"} or \code{"annotation"}, the last marking peak regions
#'   rather than plotting their scores. For \code{"interactionData"} tracks use
#'   \code{"heatmap"} or \code{"link"}.}
#'   \item{\code{color}}{\code{"character"} vector of colors. Give two values
#'   when the track has both \code{dat} and \code{dat2}, the first for the
#'   positive side and the second for the negative side of the baseline.}
#'   \item{\code{height}}{\code{"numeric"} between 0 and 1, the fraction of the
#'   figure given to this track. Note that \code{\link{trackList}} rewrites this
#'   slot for every element it is given, so set it after building the list, not
#'   before.}
#'   \item{\code{marginTop}, \code{marginBottom}}{\code{"numeric"}, space left
#'   above and below the drawing region of the track, as a fraction of its
#'   height. Defaults are 0 and 0.05.}
#'   \item{\code{ylim}}{\code{"numeric"} of length 2 fixing the y-axis range.
#'   Useful for making two tracks directly comparable, since the range is
#'   otherwise taken from the data in the current view.}
#'   \item{\code{ylabpos}}{\code{"character"}, one of \code{"left"},
#'   \code{"right"}, \code{"topleft"}, \code{"bottomleft"}, \code{"topright"},
#'   \code{"bottomright"}, \code{"abovebaseline"}, \code{"underbaseline"} or
#'   \code{"none"}; for gene-type tracks also \code{"upstream"} or
#'   \code{"downstream"}.}
#'   \item{\code{ylablas}}{\code{"numeric"} in \{0, 1, 2, 3\}, the rotation of
#'   the y label. See \code{\link[graphics]{par}:las}.}
#'   \item{\code{ylabgp}}{\code{"list"} of graphical parameters for the y label,
#'   converted to \code{\link[grid]{gpar}}, e.g. \code{list(cex=.8, col="gray30")}.}
#'   \item{\code{breaks}}{\code{"numeric"} breaks for the color key of
#'   interaction data.}
#'   \item{\code{NAcolor}}{\code{"character"}, the color used for missing cells
#'   of an interaction heatmap. Default \code{"white"}.}
#'   \item{\code{ysplit}}{\code{"numeric"}, where to split the y region for
#'   back-to-back interaction plots. Default 0.5 splits it evenly; it has no
#'   effect on other layouts.}
#' }
#'
#' The value is assigned with \code{check = TRUE}, so the validity method of
#' \code{\link{trackStyle}} runs on every call and an out-of-range value (an
#' unknown \code{tracktype}, an \code{ylablas} outside 0:3, an unrecognized
#' \code{ylabpos}) fails immediately rather than at plotting time.
#' @section Assignment in the calling frame:
#' Unlike most R functions, \code{setTrackStyleParam} modifies its first
#' argument in place: it assigns the updated track back to the variable that was
#' passed in, so \code{setTrackStyleParam(fox2, "color", "red")} is enough and
#' \code{fox2 <- setTrackStyleParam(...)} is not needed. This works only when
#' \code{ts} is given as a name or a subsettable expression that can be assigned
#' to; calling it on a temporary value, for example
#' \code{setTrackStyleParam(importScore(f), "color", "red")}, updates nothing
#' the caller can see. Inside \code{lapply} and friends the assignment targets
#' the loop variable, which is discarded, so use the returned object there
#' instead — it is returned invisibly and carries the same change.
#' @return An object of class \code{\link{track}} with the modified style,
#' returned invisibly. The variable supplied as \code{ts} is updated as a side
#' effect.
#' @examples
#' extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
#' fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED")
#'
#' ## the track is updated in place; no re-assignment needed
#' setTrackStyleParam(fox2, "tracktype", "histogram")
#' setTrackStyleParam(fox2, "color", c("#E69F00", "#56B4E9"))
#' setTrackStyleParam(fox2, "ylim", c(0, 50))
#' setTrackStyleParam(fox2, "ylabpos", "topleft")
#' setTrackStyleParam(fox2, "ylabgp", list(cex=.8, col="gray30"))
#' fox2$style$tracktype
#'
#' ## inside lapply, use the returned value
#' trs <- lapply(list(a=fox2, b=fox2), function(.ele){
#'     setTrackStyleParam(.ele, "height", .25)
#' })
#'
#' ## invalid values are rejected straight away
#' \dontrun{
#' setTrackStyleParam(fox2, "ylablas", 5)      # must be 0:3
#' setTrackStyleParam(fox2, "xscale", list())  # use setTrackXscaleParam
#' }
#' @seealso \code{\link{setTrackXscaleParam}}, \code{\link{setTrackYaxisParam}},
#' \code{\link{trackStyle}}
setMethod("setTrackStyleParam", 
          signature(ts="track", attr="character", value="ANY"),
          function(ts, attr, value){
              if(!attr %in% c("tracktype", "color", "height", "marginTop", 
                              "marginBottom", "ylim", "ylabpos", "ylablas",
                              "ylabgp", "breaks", "NAcolor", 'ysplit'))
                  stop("attr must be a slot name (except xscale and yaxis) of trackStyle.
                       try setTrackXscaleParam for xscale slot and 
                       setTrackYaxisParam for yaxis.")
              x <- ts
              slot(x@style, attr, check = TRUE) <- value
              eval.parent(substitute(ts <- x))
              return(invisible(x))
          })
#' Method setTrackXscaleParam
#' @rdname trackStyle-class
#' @exportMethod setTrackXscaleParam
#' @aliases setTrackXscaleParam
#' @aliases setTrackXscaleParam,track,character-method
setGeneric("setTrackXscaleParam", function(ts, attr, value) 
    standardGeneric("setTrackXscaleParam"))
#' @rdname trackStyle-class
#' @aliases setTrackXscaleParam,track,character,ANY-method
#' @details 
#' \code{setTrackXscaleParam} changes the x-scale bar of a track, stored as an
#' \code{\link{xscale}} object at \code{ts@@style@@xscale}. The x-scale bar is a
#' short horizontal ruler drawn inside the track's own panel, independent of
#' the shared genomic x-axis of the whole figure, typically used to show a
#' distance such as "1K" next to a zoomed-in feature.
#'
#' \code{attr} can be any slot of \code{\link{xscale}}, or the convenience
#' name \code{"position"}:
#' \describe{
#'   \item{\code{from}, \code{to}}{Objects of class \code{\link{pos}}, the two
#'   endpoints of the scale bar. Because these are \code{pos} objects rather
#'   than bare numbers, setting them directly means constructing the object
#'   yourself, e.g. \code{new("pos", x=12345678, y=0.5, unit="native")} for a
#'   point given in native (genomic) coordinates, or a fraction with
#'   \code{unit="npc"} for a position relative to the panel. Most callers will
#'   find \code{attr="position"} (below) more convenient than setting
#'   \code{from} and \code{to} separately, since it keeps them centered and in
#'   sync with the label.}
#'   \item{\code{label}}{\code{"character"}, the text drawn on the scale bar,
#'   e.g. \code{"500 bp"}. Set directly this is used verbatim.}
#'   \item{\code{gp}}{\code{"list"} of graphical parameters for the scale bar
#'   and its label, converted to \code{\link[grid]{gpar}}, e.g.
#'   \code{list(cex=.5, col="black", lwd=2)}.}
#'   \item{\code{draw}}{\code{"logical"}, whether the scale bar is drawn at
#'   all. Default \code{FALSE}; set to \code{TRUE} even \code{from}/\code{to}
#'   (or \code{position}) have been set, the scale bar will show.}
#'   \item{\code{position}}{A shortcut that sets \code{from}, \code{to} and
#'   \code{label} together from a single scale length, instead of requiring two
#'   \code{pos} objects. \code{value} must be a \code{"list"} with named
#'   elements \code{x}, \code{y} and \code{label}:
#'     \itemize{
#'       \item \code{x}, \code{y} — the native-coordinate of the bar
#'       (genomic position and track-relative y, respectively). The bar is
#'       drawn symmetrically about this point.
#'       \item \code{label} — the width the bar should represent.
#'     }
#'   A \code{label} that cannot be coerced to numeric (\code{NA} after
#'   \code{as.numeric}) is rejected before any slot is touched.}
#' }
#'
#' As with the other slots of \code{\link{xscale}}, updated values are assigned
#' with \code{check = TRUE}, so the \code{\link{pos}} validity method (checking
#' the \code{unit} string) still applies when \code{from}/\code{to} are set
#' through \code{"position"}.
#' @section Assignment in the calling frame:
#' Like \code{\link{setTrackStyleParam}}, this method assigns the updated track
#' back to the variable passed as \code{ts} in the caller's environment, so
#' \code{setTrackXscaleParam(fox2, "draw", TRUE)} is sufficient and the result
#' need not be re-assigned. It is also returned invisibly for use in contexts
#' (such as \code{lapply}) where the calling-frame assignment has no visible
#' target.
#' @return An object of class \code{\link{track}} with the modified
#' \code{xscale}, returned invisibly. The variable supplied as \code{ts} is
#' updated as a side effect.
#' @examples
#' extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
#' fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED")
#'
#' ## simplest way: give a centre point and a width in bases
#' setTrackXscaleParam(fox2, "position",
#'                      list(x=122929675, y=4, label=500))
#' setTrackXscaleParam(fox2, "draw", TRUE)
#'
#' ## style the bar itself
#' setTrackXscaleParam(fox2, "gp", list(cex=.5, col="gray30"))
#'
#' ## equivalent, done by hand with explicit pos objects
#' setTrackXscaleParam(fox2, "from",
#'                      new("pos", x=122929675-500, y=4, unit="native"))
#' setTrackXscaleParam(fox2, "to",
#'                      new("pos", x=122929675+500, y=4, unit="native"))
#' setTrackXscaleParam(fox2, "label", "1K")
#'
#' \dontrun{
#' ## rejected: label can't be parsed as a number
#' setTrackXscaleParam(fox2, "position",
#'                      list(x=122929675, y=4, label="five hundred"))
#' }
#' @seealso \code{\link{xscale}}, \code{\link{pos}},
#' \code{\link{setTrackStyleParam}},
#' \code{\link{setTrackYaxisParam}}
setMethod("setTrackXscaleParam", 
          signature(ts="track", attr="character", value="ANY"),
          function(ts, attr, value){
              if(!attr %in% c("from", "to", "label", "gp", "draw", "position"))
                  stop("attr must be a slot name of xscale object")
              x <- ts
              if(attr=="position"){
                if(!is.list(value)){
                  stop("if attr is position, value must be a list of x, y, and label.")
                }
                if(any(!c("x", "y", "label") %in% names(value))){
                  stop("if attr is position, value must be a list of x, y, and label. eg: list(x=122929375, y=0.5, label=1000)")
                }
                label0 <- as.numeric(value$label)
                if(is.na(label0)){
                  stop("label can not be converted to number.")
                }
                slot(x@style@xscale, "from", check = TRUE)  <- new("pos", x=value$x-label0/2, y=value$y, unit="native")
                slot(x@style@xscale, "to", check = TRUE)  <- new("pos", x=value$x+label0/2, y=value$y, unit="native")
                slot(x@style@xscale, "label", check = TRUE) <- convertNum2HumanNum(value$label)
              }else{
                slot(x@style@xscale, attr, check = TRUE) <- value
              }
              eval.parent(substitute(ts <- x))
              return(invisible(x))
          })
#' Method setTrackYaxisParam
#' @rdname trackStyle-class
#' @exportMethod setTrackYaxisParam
#' @aliases setTrackYaxisParam
#' @aliases setTrackYaxisParam,track,character-method
setGeneric("setTrackYaxisParam", function(ts, attr, value) 
    standardGeneric("setTrackYaxisParam"))
#' @rdname trackStyle-class
#' @aliases setTrackYaxisParam,track,character,ANY-method
#' @details 
#' \code{setTrackYaxisParam} changes the y-axis of a track, stored as a
#' \code{\link{yaxisStyle}} object at \code{ts@@style@@yaxis}. This is the axis
#' drawn alongside the track's own data panel (score, coverage, etc.), separate
#' from the shared genomic x-axis of the whole figure.
#'
#' Accepted values of \code{attr}, all slots of \code{\link{yaxisStyle}}:
#' \describe{
#'   \item{\code{at}}{\code{"numeric"} vector of y-values at which tick marks
#'   are drawn, e.g. \code{c(0, 25, 50)}. If left at its default
#'   (\code{numeric(0)}), tick positions are chosen automatically from the
#'   plotted range.}
#'   \item{\code{label}}{\code{"logical"}, whether the numeric value of each
#'   tick in \code{at} is printed next to it. Default \code{FALSE}, i.e. ticks
#'   are drawn without their values; set to \code{TRUE} to show them.}
#'   \item{\code{gp}}{\code{"list"} of graphical parameters for the axis line,
#'   ticks and labels, converted to \code{\link[grid]{gpar}}, e.g.
#'   \code{list(cex=.6, col="blue", lwd=1.5)}.}
#'   \item{\code{draw}}{\code{"logical"}, whether the y-axis is drawn at all.
#'   Default \code{TRUE}.}
#'   \item{\code{main}}{\code{"logical"}, which side of the panel the axis is
#'   drawn on: \code{TRUE} (the default) draws it on the left of the track,
#'   \code{FALSE} on the right. Useful for telling apart the axes of two
#'   adjacent or back-to-back tracks at a glance.}
#' }
#'
#' Any other value of \code{attr}, including \code{"xscale"} or a slot name
#' belonging to \code{\link{trackStyle}} itself (like \code{"color"} or
#' \code{"height"}), is rejected, use \code{\link{setTrackStyleParam}} for the
#' plain style slots and \code{\link{setTrackXscaleParam}} for the x-scale bar.
#' Values are assigned with \code{check = TRUE}, so an ill-typed value (for
#' example a non-logical passed to \code{"draw"}) fails immediately rather than
#' at plotting time.
#' @section Assignment in the calling frame:
#' Like \code{\link{setTrackStyleParam}} and \code{\link{setTrackXscaleParam}},
#' this method assigns the updated track back to the variable passed as
#' \code{ts} in the caller's environment, so
#' \code{setTrackYaxisParam(fox2, "draw", FALSE)} is sufficient on its own. It
#' is also returned invisibly for use in contexts, such as \code{lapply}, where
#' that calling-frame assignment has no visible target.
#' @return An object of class \code{\link{track}} with the modified
#' \code{yaxis}, returned invisibly. The variable supplied as \code{ts} is
#' updated as a side effect.
#' @examples
#' extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
#' fox2 <- importScore(file.path(extdata, "fox2.bed"), format="BED")
#'
#' ## show tick values at chosen positions
#' setTrackYaxisParam(fox2, "at", c(0, 25, 50))
#' setTrackYaxisParam(fox2, "label", TRUE)
#'
#' ## style the axis and move it to the right-hand side
#' setTrackYaxisParam(fox2, "gp", list(cex=.6, col="blue"))
#' setTrackYaxisParam(fox2, "main", FALSE)
#'
#' ## hide the axis entirely
#' setTrackYaxisParam(fox2, "draw", FALSE)
#'
#' \dontrun{
#' ## rejected: "color" belongs to trackStyle, not yaxisStyle
#' setTrackYaxisParam(fox2, "color", "red")
#' }
#' @seealso \code{\link{yaxisStyle}}, \code{\link{setTrackStyleParam}},
#' \code{\link{setTrackXscaleParam}}
setMethod("setTrackYaxisParam", 
          signature(ts="track", attr="character", value="ANY"),
          function(ts, attr, value){
              if(!attr %in% c("at", "label", "gp", "draw", "main"))
                  stop("attr must be a slot name of xscale object")
              x <- ts
              slot(x@style@yaxis, attr, check = TRUE) <- value
              eval.parent(substitute(ts <- x))
              return(invisible(x))
          })

#' List of tracks
#' @description An extension of List that holds only \code{\link{track}} objects. 
#' @rdname trackList-class
#' @aliases trackList
#' @seealso \code{\link{track}}.
#' @exportClass trackList

setClass("trackList", contains="list", representation(names="vector"),
         validity=function(object){
             re <- sapply(object, class)
             if(any(re!="track"))
                 return("class of elements should be track")
             return(TRUE)
         })

#' Method seqlevelsStyle<-
#' @rdname trackList-class
#' @param x trackList object.
#' @param value values to be assigned.
#' @exportMethod seqlevelsStyle<-
#' @aliases seqlevelsStyle<-,trackList-method
setReplaceMethod("seqlevelsStyle", "trackList", 
                 function(x, value){
                     for(i in seq_along(x)){
                         seqlevelsStyle(x[[i]]) <- value
                     }
                     x
                 })


#' @rdname trackList-class
#' @param \dots Each tracks in ... becomes an element in the new 
#' trackList, in the same order. This is analogous to the list constructor, except
#' every argument in ... must be derived from \code{\link{track}}.
#' @param heightDist A vector or NA to define the height of each track.
#' @export trackList
trackList <- function(..., heightDist=NA){
    listData <- list(...)
    dots <- substitute(list(...))[-1]
    names <- as.character(sapply(dots, deparse))
    if(is.na(heightDist[1])){
        heightDist <- rep(1, length(listData))
    }else{
        if(length(heightDist)!=length(listData)){
            stop("length of heightDist should be same as length of inputs")
        }
    }
    heightDist <- heightDist * 1/sum(heightDist)
    tmpData <- list()
    recursiveList <- function(tmp, tmpData, name, heightDist){
        if(is.list(tmp)){
            if(length(tmp)==0) return(tmpData)
            for(w in 1:length(tmp)) {
                each <- tmp[[w]]
                curHeightDist <- heightDist/length(tmp)
                if(is.list(each)){
                    tmpData <- recursiveList(each, tmpData, 
                                         names(tmp)[w], curHeightDist)
                }else{
                    .name <- names(tmp)[w]
                    each@style@height <- curHeightDist
                    tmpData <- c(tmpData, each)
                    names(tmpData)[length(tmpData)] <- .name
                }
            }
        }else{
            tmp@style@height <- heightDist
            tmpData <- c(tmpData, tmp)
            names(tmpData)[length(tmpData)] <- name
        }
        tmpData
    }
    if(length(listData)>=1){
        for(i in 1:length(listData)){
            tmp <- listData[[i]]
            tmpData <- recursiveList(tmp, tmpData, names[i], heightDist[i])
            rm("tmp")
        }
    }
    listData <- tmpData
    rm("tmpData")
    
    new("trackList", listData)
}