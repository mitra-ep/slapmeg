#' @title Plot the estimated random effects from SLaPMEG seperated for the study groups
#'
#' @description This plot can provide a graphical insight into the source of effect (i.e. diffrential expression)
#' in relevent pathway.
#'
#' @param obj  An slapmeg object which is the output from \code{\link{slapmeg}} or
#' \code{\link{pairslapmeg}}. Note that the \code{fullreturn=TRUE}) must have been used.
#' @param \dots optional graphical parameters can be added with an \code{+} based on \code{\link{ggplot2}}
#' structure.
#' @return returns NULL
#'
#' @author Mitra Ebrahimpoor
#'
#' \email{m.ebrahimpoor@@lumc.nl}
#'
#' @seealso
#'
#'  \code{\link{slapmeg}}, \code{\link{multslapmeg}}, \code{\link{pairslapmeg}}
#'
#' @references
#' paper title goes here
#'
#' @examples
#'
#' \dontrun{
#' # run slapmeg for a pathway
#' #' m1 <- slapmeg()
#' #Plot the estimated random effects
#' plotslapmeg(m1)
#' }
#' @export
#' @import ggplot2

plotslapmeg<-function(obj,...){

if(class(obj)!="slapmeg") stop('The object should be the output of slapmeg function!')
if(is.null(obj$EB_pred)) stop('The full output by slapmeg function should be retrieved!')

#extract data from output
plotdata<-obj$EB_pred

data_m<-melt(plotdata[,-1], id.vars = names(plotdata)[2], value.name = "obs")
#names(data_m)[1]<-"Group"
labs_group<-unique(plotdata[,2])

#plot the EBs
ggplot(data_m, aes(x=obs)) +
  facet_wrap(.~variable, ncol = 2, scales = "free")+
  geom_density(aes(fill=data_m[,1]),alpha=0.5)+
  guides(fill=guide_legend(title="Group"))+
  labs(y= "\n Density", x = "\n Estimate")+
  scale_fill_discrete(labels = paste(labs_group))+
  theme_minimal()

}
