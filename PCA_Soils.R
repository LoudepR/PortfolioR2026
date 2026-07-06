install.packages("agricolae")
library(agricolae)
library(FactoMineR)
library(factoextra)
library(cluster)
library(dplyr)

data(soil)
dat <- soil
str(dat)
summary(dat)
head(dat)

#On ne garde que les variables quantitatives.
dat_num <- dat %>% select_if(is.numeric)
dat_scaled <- scale(dat_num)
dat_scaled

#ACP
#####Explorer la structure des sols
# - identifier les gradients physico‑chimiques (acidité, fertilité, salinité…)
# - repérer les variables qui structurent le plus la variabilité
# - détecter des profils de sols (pauvres, fertiles, calcaires…)
#####Réduire la dimension
# - condenser l’information en quelques axes (2–3)
# - visualiser les relations entre variables
# - éviter les redondances (ex : Ca, Mg, CIC souvent corrélés)
#####Préparer un clustering
# - regrouper les sols en classes homogènes
# - construire une typologie agronomique


res_pca <- PCA(dat_scaled, graph = FALSE)
#eigenvalue = qtt de variance expliquée par chaque axe
res_pca$eig
	### AXE 1 : explique 33,8% de la variabilité totale
	### AXE 2 : 18,5% ainsi de suite
#Les 4 premiers axes expliquent 82.6 % de la variabilité totale
#l’ACP capture efficacement la structure multivariée du jeu de données
	#graph des dimensions ?
fviz_eig(res_pca,
         addlabels = TRUE,
         ylim = c(0, 40),
         barfill = "#377EB8",
         barcolor = "black")


#axes principaux
res_pca$var$coord
#Chaque ligne = une variable
#Chaque colonne = contribution à un axe
	#exemple : le pH est fortement corrélé à l'axe 1 (0,89)
	# MO est quasi entièrement portée par l'axe 3 (0,93)

###AXE 1 : pH, Ca, Mg, texture (sand, clay) : pH = 0.89, Ca = 0.80, clay = 0.83, sand = –0.87
#Axe 1 = gradient chimico-textural  = Sols calcaires, argileux, pH élevé ↔ sols sableux, acides.

###AXE 2 : (salinité, sodium, potassium) : EC = 0.69, Na = 0.73, K = 0.57
#Axe 2 = gradient salin / alcalinisation

###AXE 3 : MO = 0.93, CIC = 0.72, Fe = 0.72, P = 0.65
#Axe 3 = gradient de fertilité organique et minérale

###AXE 4 : rapport Ca/Mg et structure : Ca_Mg = 0.88, Mg = –0.64, Zn = 0.67
#Axe 4 = équilibre cationique Ca/Mg


############ EIGENVALUES : 
######Combien d’axes retenir ?
#4 axes : 82,6% structure claire
######L’ACP est-elle pertinente ?
#Oui : les premiers axes capturent une grande partie de la variabilité.

############ coordonnées des variables :
######Quels gradients agronomiques structurent les sols ?
#Axe 1 : chimie + texture
#Axe 2 : salinité
#Axe 3 : fertilité organique
#Axe 4 : équilibre Ca/Mg
###la fertilité, la salinité et la texture sont les trois grands gradients structurants.


#visualisation des variables
fviz_pca_var(res_pca,
             col.var = "cos2",
             gradient.cols = c("#E41A1C", "#377EB8"),
             repel = TRUE)

fviz_pca_var(res_pca,
             axes = c(1, 3),
             col.var = "cos2",
             gradient.cols = c("#E41A1C", "#377EB8"),
             repel = TRUE)


#visualisation des sols 
fviz_pca_ind(res_pca,
             geom = c("point", "text"),
             col.ind = res_pca$ind$cos2[,1],
             gradient.cols = c("#E41A1C", "#377EB8"),
             repel = TRUE,
             legend.title = "Qualité (cos²)")
#Chaque point représente un sol, et sa couleur indique la qualité de 
#représentation (cos²) dans ce plan :
	#bleu = très bien représenté
	#rouge = mal représenté
####Les sols sont dispersés dans les quatre quadrants :
# - ils ne se regroupent pas tous sur un même profil,
# - il existe une variabilité multivariée réelle,
# - les deux premiers axes capturent une part importante de cette variabilité (≈ 52 % cumulés).
##Le nuage n’est pas compact : les sols ont des signatures physico‑chimiques différentes.
#Bleu foncé : le sol est très bien expliqué par ces deux axes → interprétation fiable.
#Rouge / orange : le sol est mal représenté → son profil dépend d’axes supérieurs (Dim3, Dim4…).

#regarder d'autres axes, pour rappel avec eig :
#Dim1 = 33.8 %
#Dim2 = 18.5 %
#Dim3 = 18.1 %
#Dim4 = 12.2 %

#gradient fait : Dim1–Dim2 → gradient chimico‑textural + salinité
#a regarder : 
	#Dim1–Dim3 → gradient chimico‑textural + fertilité organique
fviz_pca_ind(res_pca,
             axes = c(1, 3),
             geom = c("point", "text"),
             col.ind = res_pca$ind$cos2[,3],   # qualité sur Dim3
             gradient.cols = c("#E41A1C", "#377EB8"),
             repel = TRUE,
             legend.title = "Qualité (cos² Dim3)")
#Dim1–Dim3 révèle les sols atypiques en matière organique.

	#Dim2–Dim3 → gradient salinité + fertilité organique
fviz_pca_ind(res_pca,
             axes = c(2, 3),
             geom = c("point", "text"),
             col.ind = res_pca$ind$cos2[,3],   # qualité de représentation sur Dim3
             gradient.cols = c("#E41A1C", "#377EB8"),
             repel = TRUE,
             legend.title = "Qualité (cos²)")




###############clustering

#Clustering hiérarchique (classification des sols)
#Distance + clustering agglomératif
dist_mat <- dist(dat_scaled, method = "euclidean")
hc <- agnes(dist_mat, method = "ward")
#dendro
fviz_dend(hc, k = NULL)

fviz_dend(hc, k = 3,
          k_colors = c("#E41A1C", "#4DAF4A", "#377EB8"))

clusters <- cutree(hc, k = 3)
fviz_cluster(list(data = dat_scaled, cluster = clusters),
             palette = c("#E41A1C", "#4DAF4A", "#377EB8"),
             ellipse.type = "convex",
             geom = "point")


#combiner pca et cluster
clusters <- cutree(hc, k = 3)
clusters <- factor(clusters)
fviz_pca_ind(res_pca,
             geom = c("point", "text"),
             col.ind = clusters,
             palette = c("#E41A1C", "#4DAF4A", "#377EB8"),
             addEllipses = TRUE,
             pointsize = 2,
             repel = TRUE,
             legend.title = "Groupes de sols")

#identifier des groupes de sols présentant des propriétés physico-chimiques 
#similaires afin de construire une typologie agronomique exploitable.

dat_scaled_df <- dat_scaled %>%
  as.data.frame() %>%
  mutate(sol = rownames(dat_scaled),
         cluster = clusters)

centroids <- aggregate(dat_scaled, by = list(cluster = clusters), FUN = mean)
parangons <- dat_scaled_df %>%
  group_by(cluster) %>%
  mutate(dist_to_centroid = sqrt(rowSums((across(where(is.numeric)) - centroids[cluster, -1])^2))) %>%
  slice_min(dist_to_centroid, n = 1) %>%   # ← le point clé
  ungroup() %>%
  select(sol, cluster, dist_to_centroid)






