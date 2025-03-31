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
sample.variables(datos_ps)

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
  
mod_datos_ps <- prune_taxa(taxa_sums(datos_ps) > 1, datos_ps) #eliminar los que no estan.
#GRAFICA: 
plot_richness(mod_datos_ps, x = "sample",color = "group", measures= c("Observed", "Shannon", "Simpson"))

  pdf("03_OUT/PLOT/DIVERSIDAD_ALFA_1.pdf", height = 8, width = 10)
  plot_richness(mod_datos_ps, x = "sample", measures= "Shannon",color = "group")
  plot_richness(mod_datos_ps, x = "sample", measures= "Simpson",color = "group")
  plot_richness(mod_datos_ps, x = "sample", measures= "Observed",color = "group")
  dev.off()

#OBSERVADOS/RIQEUZA
estimate_richness(datos_ps, measures="Observed")
#SHANNON:
estimate_richness(datos_ps, measures= "Shannon")
#SIMPSON:
estimate_richness(datos_ps, measures= "Simpson")

#ARCHIVO CSV:
metricas_alfa <- estimate_richness(datos_ps, measures= c("Observed", "Shannon", "Simpson"))
write.csv(metricas_alfa, file = "03_OUT/DATOS/METRICAS.ALFA.DIETSWAP.csv", row.names = TRUE) 

#REFERENCIA:
#https://rdrr.io/bioc/phyloseq/man/estimate_richness.html
################################################################################
#4.- Filtrado y transformación
#Aplica un filtrado para quedarte solo con los géneros más abundantes 
  #(por ejemplo, los que tienen más del 0.1% de abundancia relativa en al menos 10% de las muestras).
#INTENTO 1:
total <- sum(otu_table(datos_ps)) #para el porcentaje.
datos_filtrados <- prune_taxa(taxa_sums(datos_ps) > (0.1 * total), datos_ps) #10%
datos_filtrados 


datos_ps_2 <- transform_sample_counts(datos_ps, function(x) x/sum(x))
(datos_ps_2  <- filter_taxa(datos_ps_2, function(x) sum(x > 0.001) > (0.1*length(x)), TRUE)) #filtrado de los datos.

#REFERENCIA:
#http://www.castrolab.org/teaching/data_analysis/intro-phyloseq.html#control-de-calidad-del-an%C3%A1lisis-de-16s

################################################################################
#5.- Diversidad beta
#Realiza una ordención PCoA utilizando distancia Bray-Curtis. Usa ordinate() y 
  #plot_ordination().

pca_datos <- ordinate(datos_ps, method = "PCoA", distance = "bray")

pdf("03_OUT/PLOT/DIVERSIDAD.BETA.GRUPOS.pdf", height = 8, width = 10)
plot_1 <- plot_ordination(datos_ps, pca_datos, color = "group")+
  geom_point(size=2) + geom_path(size=0.1) + scale_colour_hue(guide = FALSE) 
plot(plot_1)
dev.off()


pdf("03_OUT/PLOT/DIVERSIDAD.BETA.SAMPLE.pdf", height = 8, width = 10)
plot_2 <- plot_ordination(datos_ps, pca_datos, color = "sample")+ #para ver por muestra
  geom_point(size=2) + geom_path(size=0.1) + scale_colour_hue(guide = FALSE) 
plot(plot_2)
dev.off()

################################################################################
#6.- Gráfica Rank-Abundance
#Crea una gráfica de abundancia-rango para observar la dominancia de taxones.
View(otu_table(datos_ps))

# PASA LOS DATOS DEL TAXA EN PORCENTAJE:
datos_modificados <- transform_sample_counts(datos_ps, function(x) 1e+02 * x/sum(x))
#LO  MANIPULA PARA GGPLOT
clusterData <- psmelt(datos_modificados) #psmelt:Melt phyloseq data object into large data.frame
#QUITA LAS QUE NO TIENEN=0
clusterData <- filter(clusterData,Abundance > 0)

#CALCULA LA MEDIA DE LSO TAXONES (OTU) Y LOS PHYLUM
# this is where the mean is calculated and the taxa to display is chosen
tax_table(datos_ps) #Phylum:
clusterAgg <- aggregate(Abundance ~ OTU + Phylum, data=clusterData, mean)

#Ordena las abunancias, aqui solo toma las 100 más abundantes.
# filtering and picking the number to display
clusterAgg <- clusterAgg[order(-clusterAgg$Abundance),][1:100,]

#grafico ggplot
pdf("03_OUT/PLOT/RANK-ABUNDANCE.PHYLUM.pdf", height = 8, width = 10)
ggplot(clusterAgg,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point(aes(color=Phylum),size=3) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10() + ggtitle("RANK-ABUNDANCE POR PHYLUM")
dev.off()

#REFERENCIA:
#https://github.com/joey711/phyloseq/issues/631 
################################################################################
#7.- Gráficas apiladas de abundancia por taxón
#Agrupa por phylum o género y grafica la composición de cada muestra como gráfica de barras apiladas.
#Puedes hacer una versión más simple con plot_bar():
View(tax_table(datos_ps))
#NO hay taxon, solo phylum
  
plot_bar(datos_ps, x = "nationality", y = "Abundance", fill ="Phylum") + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

datos_ps_apilados <-merge_samples(datos_ps, "nationality")

sample_data(datos_ps_apilados)$nationality = factor(sample_names(datos_ps_apilados))
datos_ps_apilados_2 = transform_sample_counts(datos_ps_apilados, function(x) 100 * x/sum(x))

pdf("03_OUT/PLOT/GRAFICAS.APILADAS.PHYLUM.pdf", height = 8, width = 10)
plot_bar(datos_ps_apilados_2, x = "nationality", y = "Abundance", fill ="Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  ggtitle("GRAFICAS APILADAS PHYLUM") + theme_bw() 
dev.off()

#################################################################################
################################EJERCICIO 2######################################
#GlobalPatterns.
data("GlobalPatterns")
data_gp <- GlobalPatterns

#################################################################################
#1.- Preprocesamiento
#Filtrar taxa con menos de 5 lecturas en al menos 20% de las muestras
otu_table(datos_gp)
datos_gp <- transform_sample_counts(data_gp, function(x) x/sum(x))
(datos_gp_1  <- filter_taxa(datos_gp, function(x)  sum(x > 0.02) > (0.2*length(x)), TRUE))
#2 pq a 5 ya no cuenta datos.

#Aglomerar a nivel de Familia
View(tax_table(datos_gp_1))
  #Family
rank_names(datos_gp_1)
filtrado_gp<- tax_glom(data_gp, taxrank= "Family", NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))
#NArm elimina los NA
#valores invalidos en la taxonomia
#REFERENCIA:https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/tax_glom

#Transformar a abundancias relativas (%)
relativos_gp <- transform_sample_counts(filtrado_gp, function(x) x / sum(x))
class(relativos_gp)
#Subset para incluir solo muestras de: Soil, Feces, Skin
filtrado_solo_gp <- subset_samples(relativos_gp, SampleType %in% c("Soil", "Feces", "Skin"))
filtrado_solo_gp

#FILTRADO SIN ABUNDANCIAS RELATIVAS
filtrado_solo_gp_2 <- subset_samples(filtrado_gp, SampleType %in% c("Soil", "Feces", "Skin"))
filtrado_solo_gp_2 #NO EN ABUNDANCIAS REALTIVAS.

#REFERENCIA:https://www.rdocumentation.org/packages/phyloseq/versions/1.16.2/topics/subset_samples
################################################################################
#2.-Diversidad alfa
#Calcular 3 índices de diversidad alfa (Shannon, Simpson, Observed)
#Crear boxplots comparativos de los índices entre tipos de muestra
#Realizar prueba estadística (Kruskal-Wallis) para diferencias entre grupos
sample.variables(filtrado_solo_gp_2)

alpha_div <- estimate_richness(filtrado_solo_gp_2, measures = c("Shannon","Simpson", "Observed"))
alpha_div

plot_richness(filtrado_solo_gp_2, x = "sample",color = "SampleType", measures= "Observed") + 
  geom_boxplot()

plot_richness(filtrado_solo_gp_2, x = "sample",color = "SampleType", measures= "Shannon")+ 
  geom_boxplot()

plot_richness(filtrado_solo_gp_2, x = "sample",color = "SampleType", measures= "Simpson")+ 
  geom_boxplot()

pdf("03_OUT/PLOT/GLOBAL.PATTERNS.ALFA.pdf", height = 8, width = 10)
plot_richness(filtrado_solo_gp_2, x = "sample",color = "SampleType", measures= c("Simpson","Shannon", "Observed" ))+ 
                geom_boxplot()
dev.off()

#hacerlo un data framen para poder sacar los analisis
frame_filtrado_gp <- psmelt(filtrado_solo_gp_2)
head(frame_filtrado_gp)
frame_filtrado_gp$Abundance
frame_filtrado_gp$SampleType

prubea_estadistica<- kruskal.test(Abundance~SampleType, data = frame_filtrado_gp)
print(prubea_estadistica)
#prubea_estadistica$p.value---valores de p

#por taxon:
prubea_estadistica_2 <- frame_filtrado_gp %>% group_by(OTU) %>%
  summarise(valor_p = kruskal.test(Abundance ~ SampleType)$p.value) %>% 
  filter(valor_p < 0.05) 
prubea_estadistica_2

#REFERENCIA:https://dplyr.tidyverse.org/reference/summarise.html

################################################################################
#3.-Curvas de Rango-Abundancia
#Crear gráficas de rango-abundancia para cada tipo de muestra Usar escala log10 
  #en Y Comparar patrones entre ambientes

# Sugerencia: Usar ggplot2 + geom_line(aes(...))

View(otu_table(data_gp))

# PASA LOS DATOS DEL TAXA EN PORCENTAJE:
datos_modificados_gp <- transform_sample_counts(data_gp, function(x) 1e+02 * x/sum(x))
clusterData_gp <- psmelt(datos_modificados_gp) #psmelt:Melt phyloseq data object into large data.frame

#QUITA LAS QUE NO TIENEN=0
clusterData_gp <- filter(clusterData_gp,Abundance > 0)

#CALCULA LA MEDIA DE LSO TAXONES (OTU) Y LOS PHYLUM
# this is where the mean is calculated and the taxa to display is chosen
tax_table(data_gp) #Phylum:
clusterAgg_pg <- aggregate(Abundance ~ OTU + Phylum, data=clusterData_gp, mean)

#Ordena las abunancias, aqui solo toma las 100 más abundantes.
# filtering and picking the number to display
clusterAgg_pg <- clusterAgg_pg[order(-clusterAgg_pg$Abundance),][1:100,]

ggplot(clusterAgg_pg,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point(aes(color=Phylum),size=3) + 
  geom_line(aes(group = 1), size = 0.5, color="red")+
theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10() + ggtitle("RANK-ABUNDANCE.")

#grafico ggplot
pdf("03_OUT/PLOT/RANK-ABUNDANCE.GlobalPatterns.pdf", height = 8, width = 10)
ggplot(clusterAgg_pg,aes(x=reorder(OTU,-Abundance),y=Abundance)) +
  geom_point(aes(color=Phylum),size=3) + 
  geom_line(aes(group = 1), size = 0.5, color="red")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
  scale_y_log10() + ggtitle("RANK-ABUNDANCE.")
dev.off()

################################################################################
#4.- Perfil taxonómico
#Crear gráfico apilado de abundancia a nivel de Phylum
#Mostrar solo los 5 phyla más abundantes
#Agrupar por tipo de muestra
#Usar facet_wrap para comparar ambientes
#Incluir gráficos y comentar resultados biológicos

#clusterAgg_pg <- clusterAgg_pg[order(-clusterAgg_pg$Abundance),][1:100,]

#plot_bar(agrupados_gp_2, x = "SampleType", y = "Abundance", fill ="Phylum") + 
 # geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  #ggtitle("Agrupados") + theme_bw() +
  #facet_wrap(~ Phylum)

sample.variables(data_gp)

datos_modificados_gp <- transform_sample_counts(data_gp, function(x) 1e+02 * x/sum(x))
clusterData_gp <- psmelt(datos_modificados_gp) #psmelt:Melt phyloseq data object into large data.frame

abundancia_mayores <- clusterData_gp %>% group_by(Phylum) %>%
  summarise(abundancia_promedio = mean(Abundance)) %>%
  arrange(desc(abundancia_promedio)) 
#REFERENICA:https://dplyr.tidyverse.org/reference/summarise.html

abundantes_5 <-abundancia_mayores[1:5, ] #superan el promedio y solo son los 5 que los superan

abundantes_5_filtrado <- clusterData_gp %>% filter(Phylum %in% abundantes_5$Phylum)#para que los filtre

pdf("03_OUT/PLOT/TOP5.PHYLUM.GlobalPatterns.pdf", height = 8, width = 10)
ggplot(abundantes_5_filtrado, aes(x = SampleType, y = Abundance, fill = Phylum)) + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  ggtitle("Agrupados 5") +theme_bw() +
  facet_wrap(~ Phylum)
dev.off()

#REFERENCIA:https://evomics.org/learning/genomics/metagenomics/introductory-phyloseq-plots/

################################################################################
#5.- Diversidad Beta
#Calcular distancia Bray-Curtis
#Realizar PCoA
#Visualizar con:
  #Colores por tipo de muestra
  #Elipses de confianza del 95%
#Incluir stress plot
#Realizar PERMANOVA para diferencias entre grupos
#Interpretar resultados en contexto ecológico

sample.variables(data_gp)
pca_datos_gp <- ordinate(data_gp, method = "PCoA", distance = "bray")

pdf("03_OUT/PLOT/DIVERSIDAD.BETA.GLOBAL.PATTERNS.pdf", height = 8, width = 10)
  plot_1_gp <- plot_ordination(data_gp, pca_datos_gp, color = "SampleType")+
    geom_point(size=2) + geom_path(size=0.1) + scale_colour_hue(guide = FALSE) +
    stat_ellipse(geom = "polygon", level = 0.95,
                 fill = 4, alpha = 0.25) 
  plot(plot_1_gp)
dev.off()
#REFERENCIA:https://r-charts.com/es/correlacion/grafico-dispersion-elipses-ggplot2/
# ## S3 method for class 'metaMDS'
#goodness(object, dis, ...)
## Default S3 method:
#stressplot(object, dis, pch, p.col = "blue", l.col = "red", lwd = 2, ...) 
otu_matriz <- as.matrix(otu_table(data_gp))
otu_matriz <- otu_matriz[rowSums(otu_matriz) > 0, ]#quitar todos los 0, solo mayores
distancia_otu <- vegdist(otu_matriz, method = "bray") #matriz de distancia.
head(sample_data(data_gp))

#metaMDS(otu_matriz, distance = "bray", k = 2, try = 2, trymax = 5, 
#       engine = "monoMDS")->objeto_mds
#mi computadora no pudo cargar este archivo para realizar el stressplot :/

#stressplot(object, dis, pch, p.col = "blue", l.col = "red", 
#         lwd = 2, ...) 

#REFERENCIA:https://search.r-project.org/CRAN/refmans/vegan/html/goodness.metaMDS.html

metadata <- as(sample_data(filtrado_solo_gp_2), "data.frame")

adonis2(distance(filtrado_solo_gp_2, method="bray") ~ SampleType, data = metadata)

#REFERENCIA:https://github.com/joey711/phyloseq/issues/689



