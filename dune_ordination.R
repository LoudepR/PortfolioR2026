## 1. Chargement des packages et des donnees
library(vegan)

data(dune)
data(dune.env)

## 2. NMDS (ordination non contrainte)
set.seed(123)
nmds <- metaMDS(dune, distance = "bray", k = 2, trymax = 100)
plot(nmds, type = "t")

## 3. DCA - choix entre reponse unimodale (CCA) ou lineaire (RDA)
dca_mod <- decorana(dune)
dca_mod
	#Axis lengths         3.70 (deux méthodes se valent)

## 4. CCA (reponse unimodale)
cca_mod <- cca(dune ~ Management + Use + Moisture, data = dune.env)
summary(cca_mod)
plot(cca_mod, display = c("sites", "species", "bp"),
     main = "CCA - dune ~ Management + Use + Moisture")
anova(cca_mod)
anova(cca_mod, by = "term")
anova(cca_mod, by = "axis")
scores(cca_mod, display = "species")

## 5. RDA (reponse lineaire, test de sensibilite)
rda_mod <- rda(dune ~ Management + Use + Moisture, data = dune.env)
summary(rda_mod)
plot(rda_mod, display = c("sites", "species", "bp"),
     main = "RDA - dune ~ Management + Use + Moisture")
anova(rda_mod)
anova(rda_mod, by = "term")
anova(rda_mod, by = "axis")

## 6. PERMANOVA (test global, independant de l'ordination)
perm_mod <- adonis2(dune ~ Management + Use + Moisture,
                     data = dune.env,
                     method = "bray", permutations = 999)
perm_mod

## 7. Visualisation par gradient (Moisture)
ordiplot(nmds, type = "n")
points(nmds, display = "sites",
       col = as.numeric(dune.env$Moisture),
       pch = 19)
legend("topright", legend = levels(dune.env$Moisture),
       col = 1:length(levels(dune.env$Moisture)), pch = 19,
       title = "Moisture")