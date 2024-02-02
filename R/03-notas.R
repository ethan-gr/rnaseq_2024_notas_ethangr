# NOTES 03

## ?model.matrix
mat <- with(
  trees, model.matrix(
    log(Volume) ~ log(Height) + log(Girth)
  )
)

colnames(mat)

head(
  model.matrix(
    log(trees$Volume) ~ log(trees$Height) +log(trees$Girth)
  )
)

summary(
  lm(
    log(Volume) ~ log(Height) + log(Girth),
    data = trees
  )
)


# Treatment table

(sampleData <- data.frame(
  genotype = rep(c("A", "B"), each = 4),
  treatment = rep(c("ctrl", "trt"), 4)
))

vd <- ExploreModelMatrix::VisualizeDesign(
  sampleData = sampleData,
  designFormula = ~ genotype + treatment,
  textSizeFitted = 4
)

cowplot::plot_grid(plotlist = vd$plotlist)


# fitted vakues, columnas
# variables dummy
# para variables categoricas
# se crea una columna booleana para cada uno de los posibles estados


# No soporta muchas comvinaciones lineales
# Si para alguna cosa fañltan datos son excentos de la tabla


###

# 2 Response resistant
# R: E(resistant;pre) - E(resistant;post)
#
# 3
# R: Para no estimar el intercepto

# ExperimentHub

