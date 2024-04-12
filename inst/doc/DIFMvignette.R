## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE, message = FALSE, warning = FALSE--------------------
library(DIFM)
library(knitr)

## -----------------------------------------------------------------------------
data(Violent)
data(Property)
data(WestStates)
Violent <- as.matrix(Violent)
Violent <- sqrt(Violent)
Property <- as.matrix(Property)
Property <- sqrt(Property)

## ----Violent plots, fig.width = 7, fig.height = 5.5, fig.align = "center", fig.cap = "Figure 1: Number of violent crimes"----
par(mar=c(2,2,2,2))
layout(rbind(1:4, 5:8, c(9:11,0)))
for(i in 1:11){
  plot(1960:2019, Violent[,i], main = colnames(Violent)[i], type = "l", xlab = "", ylab = "")
}

## ----Property plots, fig.width = 7, fig.height = 5.5, fig.align = "center", fig.cap = "Figure 2: Number of property crimes"----
par(mar=c(2,2,2,2))
layout(rbind(1:4, 5:8, c(9:11,0)))
for(i in 1:11){
  plot(1960:2019, Property[,i], main = colnames(Property)[i], type = "l", xlab = "", ylab = "")
}

## ---- results = "hide"--------------------------------------------------------
n.iter <- 5000
n.save <- 10
G0 <- rbind(c(1,1), c(0,1))
Violent.permutation <- permutation.order(Violent, 4)
Violent <- Violent[,Violent.permutation]
Violent.Hlist <- buildH(WestStates, Violent.permutation)

set.seed(1101)

model.attributes1V <- difm.model.attributes(Violent, n.iter, n.factors = 1, G0)
hyp.parm1V <- difm.hyp.parm(model.attributes1V, Hlist = Violent.Hlist)
ViolentDIFM1 <- DIFMcpp(model.attributes1V, hyp.parm1V, Violent, every = n.save, verbose = FALSE)
ViolentAssess1 <- marginal_d_cpp(Violent, model.attributes1V, hyp.parm1V, ViolentDIFM1, verbose = FALSE)

model.attributes2V <- difm.model.attributes(Violent, n.iter, n.factors = 2, G0)
hyp.parm2V <- difm.hyp.parm(model.attributes2V, Hlist = Violent.Hlist)
ViolentDIFM2 <- DIFMcpp(model.attributes2V, hyp.parm2V, Violent, every = n.save, verbose = FALSE)
ViolentAssess2 <- marginal_d_cpp(Violent, model.attributes2V, hyp.parm2V, ViolentDIFM2, verbose = FALSE)

model.attributes3V <- difm.model.attributes(Violent, n.iter, n.factors = 3, G0)
hyp.parm3V <- difm.hyp.parm(model.attributes3V, Hlist = Violent.Hlist)
ViolentDIFM3 <- DIFMcpp(model.attributes3V, hyp.parm3V, Violent, every = n.save, verbose = FALSE)
ViolentAssess3 <- marginal_d_cpp(Violent, model.attributes3V, hyp.parm3V, ViolentDIFM3, verbose = FALSE)

model.attributes4V <- difm.model.attributes(Violent, n.iter, n.factors = 4, G0)
hyp.parm4V <- difm.hyp.parm(model.attributes4V, Hlist = Violent.Hlist)
ViolentDIFM4 <- DIFMcpp(model.attributes4V, hyp.parm4V, Violent, every = n.save, verbose = FALSE)
ViolentAssess4 <- marginal_d_cpp(Violent, model.attributes4V, hyp.parm4V, ViolentDIFM4, verbose = FALSE)

## ---- results = "hide"--------------------------------------------------------
Property.permutation <- permutation.order(Property, 4)
Property <- Property[,Property.permutation]
Property.Hlist <- buildH(WestStates, Property.permutation)

set.seed(1101)

model.attributes1P <- difm.model.attributes(Property, n.iter, n.factors = 1, G0)
hyp.parm1P <- difm.hyp.parm(model.attributes1P, Hlist = Property.Hlist)
PropertyDIFM1 <- DIFMcpp(model.attributes1P, hyp.parm1P, Property, every = n.save, verbose = FALSE)
PropertyAssess1 <- marginal_d_cpp(Property, model.attributes1P, hyp.parm1P, PropertyDIFM1, verbose = FALSE)

model.attributes2P <- difm.model.attributes(Property, n.iter, n.factors = 2, G0)
hyp.parm2P <- difm.hyp.parm(model.attributes2P, Hlist = Property.Hlist)
PropertyDIFM2 <- DIFMcpp(model.attributes2P, hyp.parm2P, Property, every = n.save, verbose = FALSE)
PropertyAssess2 <- marginal_d_cpp(Property, model.attributes2P, hyp.parm2P, PropertyDIFM2, verbose = FALSE)

model.attributes3P <- difm.model.attributes(Property, n.iter, n.factors = 3, G0)
hyp.parm3P <- difm.hyp.parm(model.attributes3P, Hlist = Property.Hlist)
PropertyDIFM3 <- DIFMcpp(model.attributes3P, hyp.parm3P, Property, every = n.save, verbose = FALSE)
PropertyAssess3 <- marginal_d_cpp(Property, model.attributes3P, hyp.parm3P, PropertyDIFM3, verbose = FALSE)

model.attributes4P <- difm.model.attributes(Property, n.iter, n.factors = 4, G0)
hyp.parm4P <- difm.hyp.parm(model.attributes4P, Hlist = Property.Hlist)
PropertyDIFM4 <- DIFMcpp(model.attributes4P, hyp.parm4P, Property, every = n.save, verbose = FALSE)
PropertyAssess4 <- marginal_d_cpp(Property, model.attributes4P, hyp.parm4P, PropertyDIFM4, verbose = FALSE)

## -----------------------------------------------------------------------------
PDtable <- matrix(NA, 2, 4)
PDtable[1,] <- c(ViolentAssess1$Maximum, ViolentAssess2$Maximum, ViolentAssess3$Maximum, ViolentAssess4$Maximum)
PDtable[2,] <- c(PropertyAssess1$Maximum, PropertyAssess2$Maximum, PropertyAssess3$Maximum, PropertyAssess4$Maximum)
PDtable <- as.data.frame(PDtable)
rownames(PDtable) <- c("Violent", "Property")
colnames(PDtable) <- paste("Factors =", 1:4)
kable(PDtable)

## ---- fig.height = 7, fig.width = 4, fig.align = "center"---------------------
oldpar <- par(mar = c(2,2,2,2))
plot_B.CI(ViolentDIFM4, permutation = Violent.permutation)
plot_X.CI(ViolentDIFM4)
par(oldpar)

## ---- fig.height = 5, fig.width = 5, fig.align = "center"---------------------
plot_sigma2.CI(ViolentDIFM4, permutation = Violent.permutation)
plot_tau.CI(ViolentDIFM4)

## ---- fig.height = 6, fig.width = 6, fig.align = "center"---------------------
plot_B.spatial(ViolentDIFM4, WestStates, layout.dim = c(2,2))

## ---- fig.height = 7, fig.width = 4, fig.align = "center"---------------------
oldpar <- par(mar = c(2,2,2,2))
plot_B.CI(PropertyDIFM4, permutation = Property.permutation)
plot_X.CI(PropertyDIFM4)
par(oldpar)

## ---- fig.height = 5, fig.width = 5, fig.align = "center"---------------------
plot_sigma2.CI(PropertyDIFM4, permutation = Property.permutation)
plot_tau.CI(PropertyDIFM4)

## ---- fig.height = 6, fig.width = 6, fig.align = "center"---------------------
plot_B.spatial(PropertyDIFM4, WestStates, layout.dim = c(2,2))

