####################
#                  #
# Copiar todo esto #
#                  #
####################
# Hecho con gusto por Luz Yolanda Rivera Álvarez (UAEH)

# Laboratorio - MAPA DE CALOR -TERMICO- with pheatmap
# DATOS GENETICOS TOMADOS DE Sahir Bhatnagar.
# PRACTICA DE CODERS

# Objetivo: Realizar un heatmap con datos geneticos 
# ------------------------------------------------------------------------------------------------------
# En este ejercicio vamos a:
# 1. Cargar nuestra matriz hipotética de datos y dataframes adicionales
# 2. Realizar varios heatmaps

#Un mapa de calor es una representaci?n gr?fica de datos que utiliza un sistema de 
#codificaci?n de colores para representar diferentes valores


#Heatmaps with pheatmap 
#Simulated data created by Sahir Bhatnagar.

#possible data pre-processing - normalization - quantile, median, etc., log transform
#not necessary here - we have log fold change data that has already been normalized

#Calculating your distance matrix (see dist objects):
#compute how similar or different you values are
#parametric - distance measures based on Pearson correlation 
#non parametric - spearman rank - replace by ranks and calculate correlation, Kendall's - relative ordering
#euclidean - shortest distance between values (has to be normalized), takes magnitude into account
#city block/Manhattan - sum of distances along each dimension
#distance 1-correlation - of all pairs of items to be clustered

#Cluster your samples (see hclust objects):
#hierarchical, organizes into a tree structure based on similarity - short branches if similar and longer branches as similarity decreases
#repeated cycles where the 2 closest remaining items (smallest distance) get joined by a branch with the length of the branch reflecting the distance between them, the distance between this item and all other remaining items are computed until only one object remains
#single linkage clustering - distance between 2 items is the minimum of all pairwise distances between items contained in x and y - fast b/c no other calculations need to be performed once you have your distance matrix
#complete linkage is the maximum of all paiwise distances between x and y 
#average linkage - mean of all pairwise distances between items contained in x and y
#k-means organize into clusters (self-chosen number) - items are randomly assigned to a cluster - the mean vector fo rall items in each hcluster is computed, items are reassigned to the cluster whose center is closest to them - random starting points so will not always get the same answer, number of trial done to deal with the randomness
#self organizing maps



install.packages("pheatmap")
library(pheatmap)

# importar datos
file.choose()

genes <- as.matrix ( 
  read.csv("/Users/luz/Documents/Doctorado/Materias/Complejidad Económica/Entregables /Semana 7/Lab 32/heatmap_data.csv",
           sep = ",",
           header= T,
           row.names = 1))

annotation_col <- read.csv("/Users/luz/Documents/Doctorado/Materias/Complejidad Económica/Entregables /Semana 7/Lab 32/annotation_col.csv",
           header= T,
           row.names = 1)
annotation_row <- read.csv("/Users/luz/Documents/Doctorado/Materias/Complejidad Económica/Entregables /Semana 7/Lab 32/annotation_row.csv",
                           header= T,
                           row.names = 1)



#Plotting with pheatmap!, graficar con pheatmap
pheatmap(genes)

#change font, cambiar tamaño de letra
pheatmap(genes, frontsize = 6)


#default is clustering rows and columns
#cluster by gene - groups of similar genes----LOS GENES ESTAN EN LOS RENGLONES
#POR DEFAULT CLUSTEA LOS RENGLONES
pheatmap(genes, frontsize = 6, cluster_rows = TRUE, cluster_cols = FALSE)

#cluster by patient - groups of similar patients 
#DEBES HACER QUE LAS COLUMNAS SE TRANFOMEN A RENGLONES
pheatmap(genes, frontsize = 6, cluster_rows = FALSE, cluster_cols = TRUE)

#usually order by both, aparecen ambos
pheatmap(genes, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE)

#seeing some patterns emerge - but what do they mean? Great time to add annotation to our plot
#anotaciones en los renglones

pheatmap(genes, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row)

#add to row first, see that genes are clustering according to the pathways they belong to
#anotaciones en las columnas también
pheatmap(genes, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, annotation_col= annotation_col)

#now have information about the drug and condition 
#GRAFICO COMPLETO G1

pheatmap(genes, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, annotation_col= annotation_col, treeheight_row = 0, treeheight_col = 0, main= "expresión genética")

#GRAFICO QUITANDO CLUSTERS (ARBOLES DE AGRUPACI?N-DENDOGRAMAS)
#take a smaller subset, sub base de datos en base a la matriz original
sub<-genes [c (1:5, 55:60), c (1:5, 20:35 , 55:60)]

#con subset 1, le añadimos título "expresión genética"
pheatmap(sub, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, annotation_col= annotation_col, treeheight_row = 0, treeheight_col = 0, main= "expresión genética")
# con subset 2 -- DESPLEGAR VALORES
#cambiamos fuentes

pheatmap(sub, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, annotation_col= annotation_col, treeheight_row = 0, treeheight_col = 0, main= "expresión genética", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE, fontsize_number = 6)
# con color
pheatmap(sub, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, annotation_col= annotation_col, treeheight_row = 0, treeheight_col = 0, main= "expresión genética", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE, fontsize_number = 6)

#llamamos a la librería viridis y tenemos diferentes fuentes como: 
#viridis, magma, plasma, cividis, inferno
install.packages("viridis")
library(viridis)
#probamos la paleta plasma
pheatmap(sub, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, 
         annotation_col= annotation_col, treeheight_row = 0, treeheight_col = 0, 
         main= "expresión genética", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE, 
         fontsize_number = 6, col= viridis_pal(option= "plasma")(6))

#probamos la paleta magma
pheatmap(sub, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, 
         annotation_col= annotation_col, treeheight_row = 0, treeheight_col = 0, 
         main= "expresión genética", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE, 
         fontsize_number = 6, col= viridis_pal(option= "magma")(6))
#probamos con viridis
pheatmap(sub, frontsize = 6, cluster_rows = TRUE, cluster_cols = TRUE, annotation_row = annotation_row, 
         annotation_col= annotation_col, treeheight_row = 0, treeheight_col = 0, 
         main= "expresión genética", fontsize = 8, annotation_legend = FALSE, display_numbers = TRUE, 
         fontsize_number = 6, col= viridis_pal(option= "vridis")(6))

# elementos adicionales 
#cálculo de la distancia entre genes
dist(sub)
#identificar el mapa de calor de la correlación de nuestros datos
#correlación entre pacientes
pheatmap(cor(sub))

#Correlación entre genes
#con la matriz traspuesta
trans<- t(sub)
pheatmap(cor(trans))


