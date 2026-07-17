# Chargement des packages nécessaires
library(agridat)
library(lme4)
library(ggplot2)
library(dplyr)
library(viridis)
library(DHARMa)
# install.packages("SpATS")  # décommenter si le package n'est pas déjà installé
library(SpATS)

# --- Chargement et préparation des données -----------------------------
data(stroup.nin)
dat <- stroup.nin
str(dat)
summary(dat)

# Conversion du rendement en t/ha (1 bu/ac = 0.06725 t/ha pour le blé)
dat <- dat %>%
  mutate(yield_t_ha = yield * 0.06725)

# --- Visualisation spatiale globale -------------------------------------
ggplot(dat, aes(x = col, y = row, fill = yield_t_ha)) +
  geom_tile(color = "grey60", linewidth = 0.3) +
  scale_fill_viridis(
    option = "C", direction = 1, na.value = "white",
    name = "Rendement (t/ha)",
    limits = c(min(dat$yield_t_ha, na.rm = TRUE), max(dat$yield_t_ha, na.rm = TRUE)),
    guide = guide_colorbar(title.position = "top", title.hjust = 0.5,
                           barwidth = 15, barheight = 0.8)
  ) +
  coord_equal() +
  labs(title = "Répartition spatiale des rendements sur le champ expérimental",
       subtitle = "Jeu de données Stroup.nin – 56 génotypes de blé, essai de rendement",
       x = "Colonnes du champ", y = "Rangées du champ") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey30"),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold")
  )

# --- Visualisation par bloc ----------------------------------------------
ggplot(dat, aes(x = col, y = row, fill = yield_t_ha)) +
  geom_tile(color = "grey30") +
  facet_wrap(~ rep) +
  scale_fill_viridis(option = "D", na.value = "white") +
  coord_equal() +
  labs(title = "Rendement par bloc (effets spatiaux visibles)",
       x = "Colonne", y = "Rangée", fill = "Rendement (t/ha)") +
  theme_minimal()

# --- Modèle RCB de référence + diagnostic spatial (DHARMa) ---------------
d_clean <- na.omit(dat[, c("yield", "gen", "rep", "row", "col")])
mod_rcb_clean <- lmer(yield ~ gen + (1 | rep), data = d_clean)

simres <- simulateResiduals(mod_rcb_clean)
plotResiduals(simres, d_clean$row)
plotResiduals(simres, d_clean$col)
testSpatialAutocorrelation(simres, x = d_clean$col, y = d_clean$row)

# --- Modèle spatial SpATS --------------------------------------------------
# NB : le terme de bloc (rep) n'est pas inclus ici — la composante spatiale
# continue (PSANOVA) est supposée capturer la variabilité auparavant
# attribuée au bloc. Si tu veux conserver un effet de bloc distinct,
# ajoute fixed = ~ rep (les résultats changeront alors).
mod_spatial <- SpATS(response = "yield",
                     spatial = ~ PSANOVA(row, col, nseg = c(10, 10)),
                     genotype = "gen",
                     genotype.as.random = FALSE,  # explicite : on veut des BLUEs
                     data = dat)

summary(mod_spatial)

# --- Extraction des BLUEs et tableau final --------------------------------
blues <- predict(mod_spatial, which = "gen")

brut <- dat %>%
  group_by(gen) %>%
  summarise(mean_yield = mean(yield, na.rm = TRUE))

tableau_final <- brut %>%
  left_join(blues, by = "gen") %>%
  mutate(
    mean_yield_t_ha = mean_yield * 0.06725,
    predicted_t_ha = predicted.values * 0.06725,
    standard_error_t_ha = standard.errors * 0.06725
  ) %>%
  select(gen, mean_yield_t_ha, predicted_t_ha, standard_error_t_ha) %>%
  arrange(desc(predicted_t_ha))

print(tableau_final)


# ---------------------------------------------------------------------------
# COMPARAISON : modèle spatial avec vs sans effet de bloc (rep)
# ---------------------------------------------------------------------------

# --- Modèle spatial SANS effet de bloc (celui déjà utilisé) ---------------
mod_spatial_sans_rep <- SpATS(response = "yield",
                               spatial = ~ PSANOVA(row, col, nseg = c(10, 10)),
                               genotype = "gen",
                               genotype.as.random = FALSE,
                               data = dat)

# --- Modèle spatial AVEC effet de bloc fixe --------------------------------
mod_spatial_avec_rep <- SpATS(response = "yield",
                               fixed = ~ rep,
                               spatial = ~ PSANOVA(row, col, nseg = c(10, 10)),
                               genotype = "gen",
                               genotype.as.random = FALSE,
                               data = dat)

# --- Extraction des BLUEs pour les deux modèles ----------------------------
blues_sans_rep <- predict(mod_spatial_sans_rep, which = "gen") %>%
  select(gen, predicted_sans_rep = predicted.values)

blues_avec_rep <- predict(mod_spatial_avec_rep, which = "gen") %>%
  select(gen, predicted_avec_rep = predicted.values)

# --- Fusion et comparaison --------------------------------------------------
comparaison <- blues_sans_rep %>%
  left_join(blues_avec_rep, by = "gen") %>%
  mutate(
    predicted_sans_rep_t_ha = predicted_sans_rep * 0.06725,
    predicted_avec_rep_t_ha = predicted_avec_rep * 0.06725,
    ecart_t_ha = predicted_avec_rep_t_ha - predicted_sans_rep_t_ha,
    rang_sans_rep = rank(-predicted_sans_rep),
    rang_avec_rep = rank(-predicted_avec_rep),
    changement_rang = rang_avec_rep - rang_sans_rep
  ) %>%
  arrange(desc(abs(changement_rang)))

print(comparaison)

# --- Indicateurs synthétiques de comparaison --------------------------------

# Corrélation entre les deux jeux de BLUEs (proche de 1 = modèles quasi équivalents)
cor_blues <- cor(comparaison$predicted_sans_rep, comparaison$predicted_avec_rep, use = "complete.obs")
cat("Corrélation entre les deux jeux de BLUEs :", round(cor_blues, 4), "\n")

# Nombre de génotypes ayant changé de rang, et amplitude du plus grand changement
cat("Génotypes ayant changé de rang :", sum(comparaison$changement_rang != 0, na.rm = TRUE),
    "sur", nrow(comparaison), "\n")
cat("Plus grand changement de rang observé :", max(abs(comparaison$changement_rang), na.rm = TRUE), "\n")

# Comparaison des critères d'information des deux modèles (AIC/BIC équivalents pour SpATS)
cat("\n--- Comparaison des ajustements ---\n")
cat("Deviance sans rep :", mod_spatial_sans_rep$deviance, "\n")
cat("Deviance avec rep :", mod_spatial_avec_rep$deviance, "\n")

# --- Visualisation des écarts de classement ---------------------------------
ggplot(comparaison, aes(x = predicted_sans_rep_t_ha, y = predicted_avec_rep_t_ha)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Comparaison des BLUEs : modèle spatial avec vs sans effet de bloc",
       x = "Rendement ajusté sans rep (t/ha)",
       y = "Rendement ajusté avec rep (t/ha)") +
  theme_minimal(base_size = 12)