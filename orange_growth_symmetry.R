#############################################################################
# CROISSANCE DU TRONC D'ORANGERS : DYNAMIQUE SYMETRIQUE OU ASYMETRIQUE ?

# ---- Packages ----
required_packages <- c("nlme", "ggplot2", "dplyr", "brms")
to_install <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(nlme)
library(ggplot2)
library(dplyr)
library(brms)

data(Orange, package = "datasets")
str(Orange)

# ---------------------------------------------------------------------------
# 1. STRUCTURE REPETEE : MODELE NAIF VS MODELE MIXTE

ggplot(Orange, aes(x = age, y = circumference, color = Tree, group = Tree)) +
  geom_point(size = 2) +
  geom_line(alpha = 0.5) +
  labs(title = "Croissance du tronc de 5 orangers",
       x = "Age (jours depuis le 31/12/1968)",
       y = "Circonference du tronc (mm)") +
  theme_minimal()

mod_log_naif  <- nls(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)
mod_gomp_naif <- nls(circumference ~ SSgompertz(age, Asym, b2, b3), data = Orange)

mod_log_mixte <- nlme(
  circumference ~ SSlogis(age, Asym, xmid, scal),
  data = Orange,
  fixed = Asym + xmid + scal ~ 1,
  random = Asym ~ 1 | Tree,
  start = coef(mod_log_naif)
)

mod_gomp_mixte <- nlme(
  circumference ~ SSgompertz(age, Asym, b2, b3),
  data = Orange,
  fixed = Asym + b2 + b3 ~ 1,
  random = Asym ~ 1 | Tree,
  start = coef(mod_gomp_naif)
)

cat("--- Naif vs mixte : impact de la structure repetee ---\n")
print(AIC(mod_log_naif, mod_gomp_naif, mod_log_mixte, mod_gomp_mixte))

# ---------------------------------------------------------------------------
# 2. MODELE DE RICHARDS : TEST FORMEL SYMETRIQUE VS ASYMETRIQUE

#   C(t) = Asym / (1 + nu * exp(-k * (t - t0)))^(1/nu)

richards_fn <- deriv(
  ~ Asym / (1 + nu * exp(-k * (age - t0)))^(1 / nu),
  namevec = c("Asym", "t0", "k", "nu"),
  function.arg = c("age", "Asym", "t0", "k", "nu")
)

start_richards <- list(
  Asym = coef(mod_log_naif)["Asym"],
  t0   = coef(mod_log_naif)["xmid"],
  k    = 1 / coef(mod_log_naif)["scal"],
  nu   = 1
)

mod_richards_naif <- nls(
  circumference ~ richards_fn(age, Asym, t0, k, nu),
  data = Orange,
  start = start_richards
)

mod_richards_mixte <- nlme(
  circumference ~ richards_fn(age, Asym, t0, k, nu),
  data = Orange,
  fixed = Asym + t0 + k + nu ~ 1,
  random = Asym ~ 1 | Tree,
  start = unlist(start_richards)
)
summary(mod_richards_mixte)

# LRT : necessite un ajustement par maximum de vraisemblance (ML), pas
# REML, pour comparer des modeles a effets fixes differents.
mod_log_mixte_ML      <- update(mod_log_mixte,      method = "ML")
mod_richards_mixte_ML <- update(mod_richards_mixte, method = "ML")

cat("\n--- Test formel : symetrique (logistique) vs asymetrique (Richards) ---\n")
print(anova(mod_log_mixte_ML, mod_richards_mixte_ML))

cat("\nEstimation du parametre de forme nu :\n")
print(intervals(mod_richards_mixte_ML, which = "fixed"))

# ---------------------------------------------------------------------------
# 3. DIAGNOSTICS DE STRUCTURE (effets aleatoires, correlation, variance)

mod_log_mixte_full <- nlme(
  circumference ~ SSlogis(age, Asym, xmid, scal),
  data = Orange,
  fixed = Asym + xmid + scal ~ 1,
  random = pdDiag(Asym + scal ~ 1),
  start = fixef(mod_log_mixte)
)
mod_log_corAR1   <- update(mod_log_mixte, correlation = corAR1(form = ~ age | Tree))
mod_log_varPower <- update(mod_log_mixte, weights = varPower())

cat("\n--- Diagnostics complementaires ---\n")
cat("Effet aleatoire sur le rythme (scal), en plus de l'asymptote :\n")
print(anova(mod_log_mixte, mod_log_mixte_full))
cat("Correlation temporelle des residus (corAR1) :\n")
print(anova(mod_log_mixte, mod_log_corAR1))
cat("Heteroscedasticite (varPower) :\n")
print(anova(mod_log_mixte, mod_log_varPower))

# ---------------------------------------------------------------------------
# 4. TRADUCTION BIOLOGIQUE : DATE D'INFLEXION DE LA CROISSANCE

modele_retenu <- mod_log_mixte  # remplacer par mod_richards_mixte si retenu

age_grid <- seq(min(Orange$age), max(Orange$age), length.out = 500)
circ_pop <- predict(modele_retenu, newdata = data.frame(age = age_grid), level = 0)
pas <- diff(age_grid)[1]
vitesse <- (c(circ_pop[-1], NA) - c(NA, circ_pop[-length(circ_pop)])) / (2 * pas)

df_vitesse <- data.frame(age = age_grid, circonference = circ_pop, vitesse = vitesse)
age_inflexion <- df_vitesse$age[which.max(df_vitesse$vitesse)]
prop_asym <- 100 * df_vitesse$circonference[which.max(df_vitesse$vitesse)] /
  fixef(modele_retenu)["Asym"]

cat("\n--- Date d'inflexion (vitesse de croissance maximale) ---\n")
cat("Age au point d'inflexion :", round(age_inflexion, 0), "jours\n")
cat("Proportion de l'asymptote atteinte :", round(prop_asym, 1),
    "% (50% attendu si symetrique, ~37% si Gompertz)\n")

ggplot(df_vitesse, aes(x = age, y = vitesse)) +
  geom_line(linewidth = 1, color = "darkgreen") +
  geom_vline(xintercept = age_inflexion, linetype = "dashed", color = "red") +
  labs(title = "Vitesse de croissance du tronc au cours du temps",
       subtitle = paste("Point d'inflexion estime a", round(age_inflexion, 0), "jours"),
       x = "Age (jours)", y = "Vitesse de croissance (mm/jour)") +
  theme_minimal()

# ---------------------------------------------------------------------------
# 5. LES 5 ORANGERS DIFFERENT-ILS PAR LEUR TAILLE FINALE ?

start_vals <- list(
  Asym = rep(coef(mod_log_naif)["Asym"], nlevels(Orange$Tree)),
  xmid = coef(mod_log_naif)["xmid"],
  scal = coef(mod_log_naif)["scal"]
)

mod_log_asym_par_arbre <- nls(
  circumference ~ SSlogis(age, Asym[Tree], xmid, scal),
  data = Orange,
  start = start_vals
)

cat("\n--- Les arbres different-ils par leur taille finale (asymptote) ? ---\n")
print(anova(mod_log_naif, mod_log_asym_par_arbre))

cat("\nEcarts individuels a l'asymptote moyenne (effets aleatoires) :\n")
print(ranef(modele_retenu))

# ---------------------------------------------------------------------------
# 6. LIMITE : DIVERGENCE DES MODELES HORS DE LA PLAGE OBSERVEE

age_extrapole <- seq(min(Orange$age), max(Orange$age) * 1.8, length.out = 300)
extra_grid <- data.frame(age = age_extrapole)
extra_grid$log_pred  <- predict(mod_log_mixte,  newdata = extra_grid, level = 0)
extra_grid$gomp_pred <- predict(mod_gomp_mixte, newdata = extra_grid, level = 0)

ggplot(extra_grid, aes(x = age)) +
  geom_line(aes(y = log_pred),  color = "blue", linetype = "dashed", linewidth = 1) +
  geom_line(aes(y = gomp_pred), color = "red",  linewidth = 1) +
  geom_vline(xintercept = max(Orange$age), linetype = "dotted", color = "grey40") +
  geom_point(data = Orange, aes(x = age, y = circumference), inherit.aes = FALSE, alpha = 0.4) +
  labs(title = "Divergence des modeles hors de la plage de donnees observee",
       subtitle = "Ligne pointillee verticale = derniere date observee",
       x = "Age (jours)", y = "Circonference predite (mm)") +
  theme_minimal()

# ---------------------------------------------------------------------------
# 7. VALIDATION CROISEE PAR UNE APPROCHE BAYESIENNE (brms)

prior_log <- c(
  prior(normal(170, 50),  nlpar = "Asym"),
  prior(normal(700, 200), nlpar = "xmid"),
  prior(normal(300, 100), nlpar = "scal", lb = 0)
)

mod_bayes <- brm(
  bf(circumference ~ Asym / (1 + exp((xmid - age) / scal)),
     Asym ~ 1 + (1 | Tree),
     xmid ~ 1,
     scal ~ 1,
     nl = TRUE),
  data = Orange,
  prior = prior_log,
  chains = 4, iter = 4000, cores = 4,
  control = list(adapt_delta = 0.95)
)

summary(mod_bayes)

cat("\n--- Comparaison frequentiste (nlme) vs bayesien (brms) ---\n")
cat("Asymptote moyenne - nlme :", round(fixef(mod_log_mixte)["Asym"], 1), "mm\n")
cat("Asymptote moyenne - brms :", round(fixef(mod_bayes)["Asym_Intercept", "Estimate"], 1), "mm\n")