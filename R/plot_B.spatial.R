
#' @title  Spatial plots of factor loadings
#' @description  The functions builds maps of factor loadings. 
#'
#' @param Gibbs  Result of Gibbs sampler from DIFM function. 
#' @param areapoly  The polygon of the areas. We can obtain this through \code{readOGR} function from \code{sp} package. 
#' @param burnin  Number of burn-in. If not specified, it uses the first tenths as burn-in period.
#' @param permutation  Permutation of variables. If not specified, no permutation.
#' @param main.bool  Add title of the plots.
#' @param layout.dim  Dimension of panel layout for multiple factor loadings. If not specificed, factor loadings plots are layout in one column.
#'
#' @importFrom sp spplot
#' @importFrom gridExtra grid.arrange
#'
#' @return  Factor loadings map plots
#'
#' @export

plot_B.spatial <- function(Gibbs, areapoly, burnin = NA, permutation = NA, main.bool = TRUE, layout.dim = NA){
  
  if(is.na(burnin)){burnin <- round(dim(Gibbs$B)[1] / 10)}
  n.iter <- dim(Gibbs$B)[1]
  turns <- seq(burnin + 1, n.iter)
  r <- dim(Gibbs$B)[2]
  k <- dim(Gibbs$B)[3]
  countOrder <- c("First", "Second", "Third", "Fourth", "Fifth" ,"Sixth", "Seventh", "Eighth", "Nineth", "Tenth")
  repermute <- 1:r
  if( !is.na(permutation[1]) ){repermute <- order(permutation)}
  if( is.na(layout.dim[1]) ){
    layout.out <- cbind(1:k)
    layout.col <- 1}else{
    layout.out <- matrix(layout.dim[1]*layout.dim[2], layout.dim[1], layout.dim[2], byrow = TRUE)
    layout.col <- layout.dim[2]
  }

  B.Gibbs.mean <- apply(Gibbs$B[turns,repermute,], c(2,3), mean)
  # colblock <- GetColors(256, scheme = "RdYlGn")
  colblock <- colorRampPalette(c("red","white","blue"))
  theme.novpadding <-
    list(layout.heights =
           list(top.padding = 0,
                main.key.padding = 1,
                key.axis.padding = 0,
                axis.xlab.padding = 0,
                xlab.key.padding = 0,
                key.sub.padding = 0,
                bottom.padding = 0),
         layout.widths =
           list(left.padding = 1,
                key.ylab.padding = 0,
                ylab.axis.padding = 0,
                axis.key.padding = 0,
                right.padding = 0))
  # Cited from https://stackoverflow.com/questions/21579236/reset-margins-in-the-wireframe-plot-using-lattice-r-package
  
  DIFM.Bplots <- vector("list", length = k)
  names(DIFM.Bplots) <- paste("DIFM_B", 1:k, sep = "")
  
  for(i in 1:k){
    currentBname <- paste("DIFM_B", i, sep = "")
    areapoly[[currentBname]] <- B.Gibbs.mean[,i]
    colbreak <- c(seq(min(areapoly[[currentBname]]),-0.01, l = 129),
                  seq(0.01,max(areapoly[[currentBname]]), l = 128)) + c(-0.01, rep(0,255), 0.01)
    DIFM.Bplots[[currentBname]] <- spplot(areapoly, currentBname, at = colbreak, col.regions = colblock(256), 
                 main = "", par.settings = theme.novpadding)
  }
  
  do.call("grid.arrange", c(DIFM.Bplots, ncol = layout.col))

}
