###PHYLOSEQ###
################################################################################
library(phyloseq)
library(microbiome)
library(ggplot2)
library(vegan)
library(dplyr)
library(data.table)
################################################################################
#1.- Cargar los datos.
data("dietswap", package = "microbiome")
datos_ps <- dietswap
datos_ps

sample_data(datos_ps)
tax_table(datos_ps)

#################################################################################
#2.-  Curvas de rarefacción
#Genera curvas de rarefacción usando rarecurve() de vegan sobre la matriz de abundancia.

otu_rare <- otu_table(datos_ps) #crea un objeto con la table de taxones
otu_rare <- as.data.frame(t(otu_rare)) #se tiene que convertir en un data frame
sample_names <- rownames(otu_rare) #le agrega los nombres de los renglones 
#que estan presentes.

otu.rarecurve = rarecurve(otu_rare, step = 10000, col = rainbow(length(otu_rare)))
#crea la curva de rarefraccion.

  pdf("03_OUT/CUrva de rarefraccion_1.pdf", height = 8, width = 10)
  otu.rarecurve = rarecurve(otu_rare, step = 10000, col = rainbow(length(otu_rare)))
  dev.off()

#################################################################################
#3.-Diversidad alfa (α)
#Calcula y grafica  #Usa plot_richness() de phyloseq. 
  #los siguientes índices de diversidad alfa:
  #Observed (Riqueza)
  #Shannon
  #Simpson
 
