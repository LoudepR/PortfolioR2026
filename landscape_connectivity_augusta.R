#############################################################################
# ANALYSE DE PAYSAGE ECOLOGIQUE - CONNECTIVITE FONCTIONNELLE

# ---------------------------------------------------------------------------
# 1. PACKAGES ET DONNEES
# ---------------------------------------------------------------------------

# Decommenter si besoin d'installer :
# install.packages(c("landscapemetrics", "terra", "raster", "gdistance",
#                     "ggplot2", "dplyr"))

library(terra)             # gestion des rasters modernes
library(raster)            # requis pour la compatibilite avec gdistance
library(landscapemetrics)  # metriques de fragmentation
library(gdistance)         # cartes de cout et connectivite
library(ggplot2)           # visualisation
library(dplyr)             # manipulation de donnees

set.seed(42)

# Parametres d'affichage reutilises dans tout le script, pour un rendu
# graphique homogene (taille de legende, zoom sur la zone d'interet)
LEGEND_CEX <- 1.5
ZOOM_XLIM  <- c(1253000, 1264000)
ZOOM_YLIM  <- c(1253000, 1260000)

# Jeu de donnees integre : occupation du sol NLCD autour d'Augusta (USA)
# Classes typiques : foret, zones agricoles, zones urbaines, eau, etc.
# NB : les jeux de donnees de {landscapemetrics} sont stockes au format
# PackedSpatRaster (serialisable) : il faut les "depaqueter" avec
# terra::unwrap() avant toute utilisation.
landscape <- terra::unwrap(landscapemetrics::augusta_nlcd)
print(landscape)

plot(landscape, main = "Occupation du sol brute (NLCD) - Augusta",
     plg = list(cex = LEGEND_CEX),
     mar = c(3, 3, 3, 12))

# Table de correspondance des codes NLCD (simplifiee)
# 41,42,43 = foret ; 21-24 = urbain/bati ; 81,82 = agricole ; 11 = eau
freq_table <- terra::freq(landscape)
print(freq_table)

# ---------------------------------------------------------------------------
# 2. EXPLORATION DU PAYSAGE BRUT
# ---------------------------------------------------------------------------

# Verification de la coherence du raster (resolution, projection, etendue)
check_landscape(landscape)

# ---------------------------------------------------------------------------
# 3. RECLASSIFICATION : HABITAT vs NON-HABITAT
# ---------------------------------------------------------------------------
# Notre espece fictive depend de la foret continue. On cree une carte
# binaire : 1 = habitat favorable (foret), 0 = non-habitat.

codes_foret <- c(41, 42, 43)  # foret decidue, conifere, mixte (NLCD)

habitat <- terra::classify(
  landscape,
  rcl = matrix(c(codes_foret, rep(1, length(codes_foret))), ncol = 2),
  others = 0
)
names(habitat) <- "habitat"

levels(habitat) <- data.frame(
  value = c(1, 0),
  label = c("Habitat", "Matrice")
)

plot(habitat, main = "Habitat vs Matrice",
     col = c("grey85", "forestgreen"),
     plg = list(cex = LEGEND_CEX))

# ---------------------------------------------------------------------------
# 4. METRIQUES DE FRAGMENTATION
# ---------------------------------------------------------------------------

# --- Niveau paysage : vue d'ensemble de la fragmentation globale
metrics_paysage <- calculate_lsm(
  habitat,
  level = "landscape",
  metric = c("np",    # nombre total de patchs
             "pd",    # densite de patchs
             "shdi",  # diversite de Shannon
             "ed")    # densite de lisiere (edge density)
)
print(metrics_paysage)

# --- Niveau classe : on se concentre sur la classe "habitat" (valeur 1)
metrics_classe <- calculate_lsm(
  habitat,
  level = "class",
  metric = c("np", "area_mn", "area_cv", "clumpy", "cohesion")
) %>%
  filter(class == 1)
print(metrics_classe)

# --- Niveau patch : liste de tous les patchs individuels d'habitat
metrics_patch <- calculate_lsm(
  habitat,
  level = "patch",
  metric = c("area", "core", "shape")
) %>%
  filter(class == 1) %>%
  arrange(desc(value))
head(metrics_patch, 10)

# Visualisation des patchs d'habitat individuels
show_patches(habitat, class = 1, labels = FALSE)

# ---------------------------------------------------------------------------
# 5. CARTE DE COUT (RESISTANCE AU DEPLACEMENT)
# ---------------------------------------------------------------------------
# Hypothese ecologique : la foret = deplacement facile (cout faible),
# les milieux ouverts/agricoles = deplacement moderement couteux,
# l'urbain = quasi-infranchissable (cout tres eleve).

rcl_cout <- matrix(c(
  41, 1,   42, 1,   43, 1,              # foret : cout faible
  81, 5,   82, 5,                       # agricole : cout moyen
  21, 50,  22, 80,  23, 100, 24, 150,   # urbain croissant : cout fort
  11, 1000                              # eau : quasi infranchissable
), ncol = 2, byrow = TRUE)

cout <- terra::classify(landscape, rcl = rcl_cout, others = 10)
names(cout) <- "cout_resistance"

plot(cout, main = "Carte de cout de deplacement (resistance)",
     plg = list(cex = LEGEND_CEX))

# ---------------------------------------------------------------------------
# 6. CONNECTIVITE FONCTIONNELLE (gdistance)
# ---------------------------------------------------------------------------
cout_raster <- raster::raster(cout)

# La transition layer encode la conductance (inverse du cout) entre
# cellules voisines ; fonction de conductance moyenne harmonique,
# classique en ecologie du paysage (directions = 8).
tr <- transition(cout_raster, transitionFunction = function(x) 1 / mean(x),
                  directions = 8)

# Deux corrections geometriques distinctes selon l'usage en aval :
tr_c <- geoCorrection(tr, type = "c")  # pour costDistance() / shortestPath()
tr_r <- geoCorrection(tr, type = "r")  # pour passage()

# --- Labellisation directe des patchs d'habitat (8-connectivite) -----------
forest_mask <- terra::classify(habitat, cbind(0, NA))       # ne garde que habitat = 1
patch_id_r  <- terra::patches(forest_mask, directions = 8)  # 1 ID unique par patch connecte

# --- Selection des 3 plus grands patchs -------------------------------------
patch_sizes <- terra::freq(patch_id_r)
patch_sizes <- patch_sizes[order(-patch_sizes$count), ]
top3_ids <- patch_sizes$value[1:3]

# --- Centroide = moyenne des coordonnees des cellules de chaque patch ------
vals <- terra::values(patch_id_r)[, 1]
xy   <- terra::xyFromCell(patch_id_r, which(!is.na(vals)))
df   <- data.frame(id = vals[!is.na(vals)], x = xy[, 1], y = xy[, 2])

pts_df <- df %>%
  filter(id %in% top3_ids) %>%
  group_by(id) %>%
  summarise(x = mean(x), y = mean(y), .groups = "drop")

pts <- as.matrix(pts_df[, c("x", "y")])
print(pts)  # verification : doit contenir exactement 3 lignes distinctes

# --- Distances entre les 3 plus grands massifs forestiers -------------------
dist_cout <- costDistance(tr_c, pts)
print(dist_cout)

dist_euclid <- dist(pts)
print(dist_euclid)

# ---------------------------------------------------------------------------
# 6.1 CARTE DE PASSAGE ET SENSIBILITE AU PARAMETRE THETA
# ---------------------------------------------------------------------------
# theta pondere l'arbitrage entre deux comportements de deplacement :
#   - theta faible  -> dispersion large, proche d'une marche aleatoire
#   - theta eleve   -> trajets optimises vers la route la moins couteuse
# Calibrage sur l'ordre de grandeur reel de dist_cout (~10 000-15 000) :
# theta doit rester dans une plage ou theta x cout ~ 1-10 pour eviter les
# valeurs numeriquement nulles (theta trop eleve) ou une carte non
# discriminante (theta trop faible).
summary(as.vector(dist_cout))

passage_map <- passage(tr_r, origin = pts[1, ], goal = pts[2, ], theta = 1e-4)
plot(passage_map, main = "Corridors de connectivite potentielle (theta = 1e-4)")

# Test de robustesse sur une plage de theta
thetas <- c(1e-5, 1e-4, 1e-3)
par(mfrow = c(1, 3))
for (th in thetas) {
  pm <- passage(tr_r, origin = pts[1, ], goal = pts[2, ], theta = th)
  plot(pm, main = paste("theta =", th))
}
par(mfrow = c(1, 1))

# Rehaussement visuel du corridor pour la carte de synthese finale
# (log1p + mise a l'echelle pour ameliorer le contraste des faibles valeurs)
passage_boost <- log1p(passage_map * 10000)

# ---------------------------------------------------------------------------
# 6.2 CHEMIN DE MOINDRE COUT (LCP) ET CARTE DE SYNTHESE
# ---------------------------------------------------------------------------
lcp <- shortestPath(tr_c, pts[1, ], pts[2, ], output = "SpatialLines")

# Carte de synthese (vue d'ensemble) : fond de corridor sous le LCP
# pour garantir que le trajet reste toujours visible au premier plan
plot(habitat, col = c("grey90", "forestgreen"), main = "Corridor et LCP - vue d'ensemble")
plot(passage_boost, add = TRUE, alpha = 0.5)
lines(lcp, col = "red", lwd = 3)

# Carte de synthese (zoom sur la zone d'interet)
# NB : les rasters sont recadres en amont avec crop() plutot que par
# xlim/ylim au moment du plot(), pour eviter tout desalignement entre
# les couches empilees (raster categoriel + raster continu + vecteur).
zoom_ext <- terra::ext(ZOOM_XLIM[1], ZOOM_XLIM[2], ZOOM_YLIM[1], ZOOM_YLIM[2])
habitat_zoom <- terra::crop(habitat, zoom_ext)
passage_zoom <- terra::crop(rast(passage_boost), zoom_ext)

plot(habitat_zoom, col = c("grey90", "forestgreen"),
     main = "Corridor et LCP - zoom zone d'interet",
     plg = list(cex = LEGEND_CEX))
plot(passage_zoom, add = TRUE, alpha = 0.5)
lines(lcp, col = "red", lwd = 3)

# ---------------------------------------------------------------------------
# 7. COMPARAISON DISTANCE EUCLIDIENNE vs DISTANCE DE COUT
# ---------------------------------------------------------------------------
# Un ratio >> 1 indique que la fragmentation / resistance du paysage
# allonge fortement le trajet reel par rapport a la ligne droite :
# signal d'une connectivite fonctionnelle degradee.

ratio_connectivite <- as.numeric(dist_cout) / as.numeric(dist_euclid)
cat("Ratio distance de cout / distance euclidienne :\n")
print(round(ratio_connectivite, 2))

