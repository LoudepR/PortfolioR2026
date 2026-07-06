###############################################
# GLM avancé – Effet des traitements insecticides
# Jeu de données : beall.webworms (agridat)
###############################################

library(agridat)
library(dplyr)
library(ggplot2)
library(car)
library(MASS)
library(emmeans)
library(viridis)

#Chargement des données
data("beall.webworms")
dat <- beall.webworms
str(dat)
summary(dat)
head(dat)

# Variables :
# y = nombre de webworms (comptage)
# spray = traitement insecticide (oui/non)
# lead = traitement arsenical (oui/non)
# trt = combinaison des traitements (4 niveaux)
# block = bloc expérimental
# row, col = structure spatiale

#C’est un comptage, on ne peut utiliser un modèle linéaire classique (lm) → GLM avec une distribution adaptée (Poisson, NB, etc).

#Visualiser les données
unique(dat[, c("trt", "spray", "lead")])
#T1 aucun traitement (témoin : aucune protection)
#T2 spray seul : action rapide, effet de contact
#T3 lead arsenate seul : effet plus lent, résiduel
#T4 combinaison : effet synergique, efficacité maximale

#Visualisation spatiale de l'infestation sur le champ
ggplot(dat, aes(x = row, y = col, fill = y)) +
  geom_tile() +
  scale_fill_viridis(option = "magma", name = "Webworms") +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(title = "Cartographie de l'infestation de webworms",
       x = "Ligne",
       y = "Colonne")

#Visualisation par bloc (RCB)
ggplot(dat, aes(x = col, y = row, color = y)) +
  geom_point(size = 3) +
  scale_color_viridis(option = "plasma", name = "Webworms") +
  coord_equal() +
  facet_wrap(~block) +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution spatiale des webworms par bloc",
       x = "Colonne",
       y = "Ligne")
#### PAS LISIBLE

# Filtrer uniquement par bloc
ggplot(dat |> filter(block == "B1"), 
  aes(x = col, y = row, color = y)) +
  geom_point(size = 3) +
  scale_color_viridis(option = "plasma", name = "Webworms") +
  coord_equal() +
  theme_minimal(base_size = 14) +
  labs(title = "Distribution spatiale des webworms – Bloc 1",
       x = "Colonne",
       y = "Ligne")

—---------------------------------
#BOXPLOT - Abondance des webworms selon les traitements
ggplot(dat, aes(x = trt, y = y, fill = trt)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(title = "Abondance des webworms selon les traitements",
       x = "Traitement",
       y = "Nombre de webworms")

#Boxplot final avec résultats GLMnb
ggplot(dat, aes(x = trt, y = y, fill = trt)) +
  geom_boxplot(alpha = 0.7) +
  geom_text(aes(x = "T1", y = 8, label = "a"), size = 6) +
  geom_text(aes(x = "T2", y = 8, label = "c"), size = 6) +
  geom_text(aes(x = "T3", y = 8, label = "b"), size = 6) +
  geom_text(aes(x = "T4", y = 8, label = "c"), size = 6) +
  geom_text(aes(x = "T1", y = 7.5, label = "***"), size = 5) +
  geom_text(aes(x = "T2", y = 7.5, label = "***"), size = 5) +
  geom_text(aes(x = "T3", y = 7.5, label = "***"), size = 5) +
  geom_text(aes(x = "T4", y = 7.5, label = "ns"), size = 5) +
  scale_fill_manual(
    values = c(
      "T1" = "#E41A1C",  # rouge
      "T2" = "#4DAF4A",  # vert
      "T3" = "#377EB8",  # bleu
      "T4" = "#984EA3"   # violet
    )
  ) +
  scale_x_discrete(labels = c(
    "T1" = "Témoin",
    "T2" = "Spray",
    "T3" = "Lead",
    "T4" = "Spray+Lead"
  )) +
  labs(
    x = "Traitement",
    y = "Nombre de webworms"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 16, face = "bold")
  )


####### analyse boxplot ########

#T1 (témoin) :
#→ médiane la plus élevée, grande dispersion, plusieurs valeurs extrêmes.
#→ forte infestation, absence de traitement.
#T2 (spray seul) :
#→ médiane plus basse, distribution resserrée.
#→ effet insecticide visible, réduction de l’abondance.
#T3 (lead arsenate seul) :
#→ médiane encore plus basse que T2.
#→ effet résiduel du plomb arsenical, efficace contre les larves.
#T4 (spray + lead) :
#→ médiane la plus faible, très peu de dispersion.
#→ combinaison synergique : efficacité maximale, infestation quasi nulle.

---------------

#Commencer par modele GLM poisson parce que c’est le modèle de référence 
#pour des comptages indépendants.

#GLM Poisson avec interaction spray × lead
mod_pois <- glm(y ~ spray * lead + block,
                data = dat,
                family = poisson())
summary(mod_pois)
#En R, spray * lead = spray + lead + spray:lead  
#→ le terme spray:lead représente l’effet spécifique de la combinaison.

####### analyse GLM POISSON ########

#INTERCEPT = le log du nombre moyen de webworms dans la condition reference 
	#reference : B bloc1 - aucun traitement 
#0.664 -> exp(0.66459)≈1.94 webworms

#EFFETS PRINCIPAUX
#sprayY = -1.02043 -  Très significatif (p < 2e-16).
#exp(−1.02043)≈0.36 
#Le spray réduit l’abondance à 36 % du témoin → −64 % d’infestation.
#LeadY = -0.49628  -  Très significatif (7.41e-11).
#exp(-0.49628)≈0.36 
#Le spray réduit l’abondance à 61 % du témoin → −39 % d’infestation.
#sprayY:leadY = 0.29425 -  Très significatif (p = 0.034).
#Ce coefficient est positif, ce qui peut surprendre.
#Mais attention : il s’ajoute aux effets principaux.
#AINSI - (combinaison → (spray + lead + interaction))
#exp(−1.02043−0.49628+0.29425)= exp(−1.22246)≈0.295
#Le spray réduit l’abondance à 29 % du témoin → −71 % % d’infestation.
#L'interaction des deux traitements reste meilleure que chaque traitement seul

#Effet des blocs (variabilité spaciale);
#Plusieurs blocs significativement négatifs (B11, B12, B13, B2, B3, B4, B5, B9)
#→ infestation plus faible que dans le bloc de référence (B1).
#Quelques blocs non significatifs (B10, B6, B7, B8)
#→ infestation comparable au bloc de référence.
#Il existe une forte hétérogénéité spatiale dans le champ.
#C’est normal pour un ravageur qui forme des foyers d’infestation.

#Null deviance: 1955.9
#Residual deviance: 1598.4
#Le modèle explique une partie de la variabilité, mais la deviance résiduelle 
#est encore élevée.
	##########Indice de surdispersion à vérifier 

—-------

#Vérification de la surdispersion
dispersion <- sum(residuals(mod_pois, type = "pearson")^2) / mod_pois$df.residual
dispersion
#≈1	OK	Poisson adapté
#>1.5	Surdispersion	Passer à quasi-Poisson ou NB
#<1	Sous-dispersion	Rare, vérifier les données
# Si surdispersion > 1.5 → quasi-Poisson ou NB
			#1.25 -> poisson possible mais privilégier nb

#ajuster le modèle à la variabilité réelle.
#CHOIX NB
mod_nb <- glm.nb(y ~ spray * lead + block, data = dat)
summary(mod_nb)
AIC(mod_pois, mod_nb)
#L’AIC du modèle NB est plus bas (2990 < 3027) → il s’ajuste mieux aux données.

	####### analyse GLM NB ########	

#EFFETS PRINCIPAUX
#sprayY = −1.01067 -  Très significatif (p < 2e-16).
#exp(−1.01067)≈0.36 
#Le spray réduit l’abondance à 36 % du témoin → −64 % d’infestation.
#LeadY = −0.50071  -  Très significatif (2.48e‑08).
#exp(−0.50071)≈0.61 
#Le spray réduit l’abondance à 61 % du témoin → −39 % d’infestation.
#sprayY:leadY = +0.29088 -  peu significatif (p = 0.0596).
#Ce coefficient est positif, ce qui peut surprendre.
#Mais attention : il s’ajoute aux effets principaux.
#AINSI - (combinaison → (spray + lead + interaction))
#exp(−1.01067−0.50071+0.29088)= exp(−1.2205)≈0.295
#Le spray réduit l’abondance à 29 % du témoin → −71 % % d’infestation.
#L'interaction des deux traitements reste meilleure que chaque traitement seul
#L’interaction spray × lead présente une tendance à la significativité 
#(p = 0.0596), suggérant un effet combiné potentiellement complémentaire entre 
#les deux traitements.
#Bien que marginale, cette tendance reste biologiquement cohérente avec 
#l’efficacité observée sur le terrain.

—---------

#FAIRE UN POST-HOC
#OBJ d’identifier précisément quelles modalités diffèrent entre elles après 
#le GLM, car le modèle indique qu’il existe un effet global, mais pas quelles 
#comparaisons sont significatives. 
Anova(mod_nb, type = 2)

posthoc <- emmeans(mod_nb, pairwise ~ spray * lead)
posthoc
#Lecture des moyennes (emmeans) :
#Y Y < Y N < N Y < N N → la combinaison est la plus efficace.

#Lecture des contrastes (comparaisons Tukey)
#Les contrastes testent les différences entre traitements
#| Comparaison | z.ratio | p.value  | Interprétation |
#| N N − Y N   | 9.865   | < 0.0001 | spray seul - significatif |
#| N N − N Y   | 5.574   | < 0.0001 | plomb seul - significatif |
#| N N − Y Y   | 11.183  | < 0.0001 | combinaison - très significative |
#| Y N − N Y   | −4.668  | < 0.0001 | spray - plus fort que plomb |
#| Y N − Y Y   | 1.671   | 0.3392   | combinaison ≈ spray seul (différence non significative) |
#| N Y − Y Y   | 6.230   | < 0.0001 | combinaison - plus forte que plomb seul |
#Le test post‑hoc (méthode de Tukey) montre que tous les traitements diffèrent 
#significativement du témoin. Le spray seul et la combinaison spray + lead 
#entraînent les plus fortes diminutions du nombre de webworms. La combinaison 
#n’est pas significativement plus efficace que le spray seul, mais elle reste 
#la plus performante globalement, confirmant un effet complémentaire entre les 
#deux produits.

—------

#Predictions valeurs glm nb
#Convertir les coefficients du modèle (sur l’échelle du log) en nombre moyen 
#de webworms.
#spray seul → ~36 % du témoin
#lead seul → ~61 % du témoin
#combinaison → ~29 % du témoin
# les valeurs exactes, en unités biologiques (nombre de larves).
#Fixé sur bloc1 (reference)

#Cette étape sert à transformer les coefficients du modèle en valeurs 
#biologiquement interprétables (nombre moyen de webworms) et à visualiser 
#clairement l’effet des traitements. Les prédictions permettent de comparer les 
#traitements sur l’échelle naturelle de la réponse, tandis que le graphique 
#synthétise les résultats du modèle négatif binomial de manière intuitive et 
#pédagogique.

newdata <- expand.grid(
  spray = levels(dat$spray),
  lead = levels(dat$lead),
  block = levels(dat$block)[1])
newdata$pred <- predict(mod_nb, newdata = newdata, type = "response")
newdata$pred

#(pour retrouver les données à la main)
exp(0.65857)
exp(0.65857 - 1.01067)
exp(0.65857 - 0.50071)
exp(0.65857 - 1.01067 - 0.50071 + 0.29088)

#Les prédictions du modèle négatif binomial sont obtenues en additionnant les 
#coefficients correspondant à chaque combinaison de traitements (spray, lead, 
#interaction), puis en appliquant l’exponentielle pour revenir sur l’échelle 
#naturelle du nombre de webworms. Les valeurs obtenues (1.93, 0.70, 1.17, 0.57) 
#représentent les abondances moyennes attendues pour les quatre combinaisons de 
#traitements.

#Visualisation des prédictions
# Créer une colonne traitement simple
newdata$traitement <- factor(
  paste(newdata$spray, newdata$lead),
  levels = c("N N", "Y N", "N Y", "Y Y"),
  labels = c("Témoin", "Spray", "Lead", "Combinaison")
)

# Barplot
ggplot(newdata, aes(x = traitement, y = pred, fill = traitement)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(pred, 2)), vjust = -0.7, size = 5) +
  scale_fill_manual(
    values = c(
      "Témoin"      = "#E41A1C",  # rouge
      "Spray"       = "#4DAF4A",  # vert
      "Lead"        = "#377EB8",  # bleu
      "Combinaison" = "#984EA3"   # violet
    )
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  ) +
  labs(
    title = "Nombre de webworms prédit par traitement (GLM NB)",
    y = "Webworms prédits"
  ) +
  ylim(0, max(newdata$pred) * 1.25)



—---------------------
install.packages("AER")
library(AER)
disp_test <- dispersiontest(mod_pois)
disp_test

install.packages("COMPoissonReg")
library(COMPoissonReg)
mod_com <- glm.cmp(y ~ spray * lead + block, data = dat)
summary(mod_com)

dat %>%
  group_by(spray, lead) %>%
  summarise(
    Moy = mean(y),
    Var = var(y),
    ID = Var / Moy
  )
#Si l’indice de dispersion (ID) diminue dans les groupes traités,
#cela prouve que les traitements réduisent l’agrégation.

#L’indice de dispersion (ID) diminue nettement entre le témoin (1.66) et les 
#traitements (≈ 1.2–1.3).
#Or, ID>1 indique une agrégation (surdispersion), et plus ID est proche de 1, 
#plus la population est homogène.
#Les insecticides tuent les individus dans les zones denses → les survivants 
#sont plus dispersés.
#La population devient moins agrégée, donc plus homogène spatialement.







