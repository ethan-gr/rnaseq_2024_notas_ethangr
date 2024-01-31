# NOTES 02
# `mailto:` for emails in postcards


## SummarizedExperiment

# Object <- colData(samples): A set of data that refers to the sample
# Object <- rowRanges(genes):
# Object <- assay(s): gene;sample: Se can have more than one table, more than one combination for gene;sample-values (e.g. counts)
# Object <- exptData: metadata
# `se <- summarizedExperiment(assays, rowData, colData, exptData)`


###### \begin{example}

## Lets build our first SummarizedExperiment object
library("SummarizedExperiment")
## ?SummarizedExperiment

## De los ejemplos en la ayuda oficial

## Creamos los datos para nuestro objeto de tipo SummarizedExperiment
## para 200 genes a lo largo de 6 muestras
nrows <- 200
ncols <- 6
## Números al azar de cuentas
set.seed(20210223)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
## Información de nuestros genes
rowRanges <- GRanges(
  rep(c("chr1", "chr2"), c(50, 150)),
  IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),   # genera valores al azar
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))
## Información de nuestras muestras
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)
## Juntamos ahora toda la información en un solo objeto de R
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

## Exploremos el objeto resultante
rse

dim(rse)
assayNames(rse)

head(rse)

rowData(rse)

rse$Treatment

###### \end{example}

# El espacio en código nos ayuda a entender de qué se trata cada elemento, porque esta en un espacio para ello
# Estoy de acuerdo

# Cuando se hace un subste se duplica en la memoria




## iSEE

# Genera una pagina web con displays que nos mustran los datos del objeto SummarizedExperiment

sce_layer <- spatialLIBD::fetch_data("sce_layer")
iSEE::iSEE(sce_layer)

## reCount

# Cuenta con datasets curados para el analisis
# Tiene diferentes propiedades de interéslibrar
library('recount3')
human_projects <- available_projects()
head(human_projects)

proj_info <- subset(human_projects, project == "SRP009615" & project_type == "data_sources")

rse_gene_SRP009615 <- create_rse(proj_info)

assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)

rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)

iSEE::iSEE(rse_gene_SRP009615)


## SHINYAPPS

## Para hacer
## https://libd.shinyapps.io/SRP009615/

## Primero necesitamos configurar RStudio con
## shinyapps.io. Para eso necesitaremos:
# install.packages("rsconnect")

## También necesitamos verificar que tengamos todos
## los paquetes en versiones nuevas. Eso
## lo podemos hacer con:
# BiocManager::valid()

## Después necesitamos copiar y pegar la información
## de nuestra cuenta (numéro y token de acceso)

## Ahora si ya podemos continuar
options(repos = BiocManager::repositories())

library("recount3")


## URL de recount3
# options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")

## ----'quick_example'--------------------------------
## Revisemos todos los proyectos con datos de humano en recount3
# human_projects <- available_projects()
#
# ## Encuentra tu proyecto de interés. Aquí usaremos
# ## SRP009615 de ejemplo
# proj_info <- subset(
#     human_projects,
#     project == "SRP009615" & project_type == "data_sources"
# )
# dput(proj_info)

## Crear la info del proyecto de forma manual
proj_info <- structure(list(
  project = "SRP009615", organism = "human", file_source = "sra",
  project_home = "data_sources/sra", project_type = "data_sources",
  n_samples = 12L
), row.names = 1838L, class = "data.frame")

## Crea un objeto de tipo RangedSummarizedExperiment (RSE)
## con la información a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

## ----"tranform_counts"------------------------------
## Convirtamos las cuentas por nucleotido a cuentas por lectura
## usando compute_read_counts().
## Para otras transformaciones como RPKM y TPM, revisa transform_counts().
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)


## ----"expand_attributes"----------------------------
## Para este estudio en específico, hagamos más fácil de usar la
## información del experimento
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)

## Crear el sitio web interactivo
iSEE::iSEE(rse_gene_SRP009615)
