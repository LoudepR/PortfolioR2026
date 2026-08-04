# ============================================================================
# CMR Jolly-Seber -- bunting (package Rcapture)
# Question de recherche : une fois les individus transitoires exclus, la
# survie hivernale du bruant azure est-elle constante dans le temps, et
# quelle valeur prend-elle ?
# ============================================================================

# 1. Setup et chargement des donnees
library(Rcapture)

data(bunting)
str(bunting)
head(bunting)

n_occasions <- 8
n_histoires_possibles <- 2^n_occasions - 1
cat("Nombre d'historiques de capture possibles (hors 00000000):",
    n_histoires_possibles, "\n")
cat("Nombre d'individus observes au total:", sum(bunting[, "freq"]), "\n")

# 2. Exploration descriptive
desc <- descriptive(bunting, dfreq = TRUE)
desc

n_col <- 15
n_row <- 7

waffle_df <- data.frame(
  x = rep(1:n_col, times = n_row),
  y = rep(n_row:1, each = n_col)
)
waffle_df <- waffle_df[1:100, ]

waffle_df$categorie <- factor(
  c(rep("1 fois (transitoire)", 85),
    rep("2-3 fois", 13),
    rep("4+ fois (noyau resident)", 2)),
  levels = c("1 fois (transitoire)", "2-3 fois", "4+ fois (noyau resident)")
)

p1 <- ggplot(waffle_df, aes(x = x, y = y, fill = categorie)) +
  geom_tile(color = "white", linewidth = 1.5, width = 0.9, height = 0.9) +
  scale_fill_manual(values = c(
    "1 fois (transitoire)" = "#d9e8df",
    "2-3 fois" = "#7fb69e",
    "4+ fois (noyau resident)" = "#1b4332"
  )) +
  coord_equal() +
  labs(title = "Sur 100 bruants observes,\nseuls 2 sont de vrais residents fideles",
       fill = NULL) +
  theme_void(base_size = 11) +
  theme(legend.position = "bottom",
      legend.text = element_text(size = 12),
      legend.key.size = unit(1, "cm"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12))

p1


# 3. Comparaison naive : Lincoln-Petersen (Chapman) entre occasions consecutives
chapman_lp <- function(n1, n2, m2) {
  N_hat <- ((n1 + 1) * (n2 + 1) / (m2 + 1)) - 1
  var_N <- ((n1 + 1) * (n2 + 1) * (n1 - m2) * (n2 - m2)) / (((m2 + 1)^2) * (m2 + 2))
  data.frame(N_hat = N_hat, se = sqrt(var_N))
}

lp_results <- data.frame(
  periode = character(),
  n1 = integer(), n2 = integer(), m2 = integer(),
  N_hat = numeric(), se = numeric(),
  stringsAsFactors = FALSE
)

hist_mat <- as.matrix(bunting[, paste0("p", 1:n_occasions)])
freq_vec <- bunting[, "freq"]

for (i in 1:(n_occasions - 1)) {
  n1 <- sum(freq_vec[hist_mat[, i] == 1])
  n2 <- sum(freq_vec[hist_mat[, i + 1] == 1])
  m2 <- sum(freq_vec[hist_mat[, i] == 1 & hist_mat[, i + 1] == 1])
  est <- chapman_lp(n1, n2, m2)
  lp_results <- rbind(lp_results, data.frame(
    periode = paste0(1972 + i, "-", 1972 + i + 1),
    n1 = n1, n2 = n2, m2 = m2,
    N_hat = est$N_hat, se = est$se
  ))
}

lp_results

# 4. Modele de Jolly-Seber -- version de base (m1)
op_m1 <- openp(bunting, dfreq = TRUE)
op_m1$model.fit
op_m1$survivals
op_m1$N
op_m1$birth

plot(op_m1)


library(ggplot2)

phi_df <- data.frame(
  transition = factor(c("1->2","2->3","3->4","4->5","5->6","6->7"),
                       levels = c("1->2","2->3","3->4","4->5","5->6","6->7")),
  estimate = c(0.4882479, 0.2227044, 0.5511586, 0.3237007, 0.4082043, 0.3329193),
  se = c(0.14234679, 0.04924800, 0.12811812, 0.05480251, 0.05892365, 0.04426334)
)

ggplot(phi_df, aes(x = transition, y = estimate)) +
  geom_pointrange(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se),
                   color = "#2e8b57", linewidth = 0.8) +
  labs(title = "Survie estimee tres instable selon les periodes",
       x = "Transition", y = "Survie estimee (phi)") +
  theme_minimal(base_size = 12)



# 5. Hypothese de transience : exclusion des individus captures une seule fois
keep_no_transient <- rowSums(hist_mat) > 1
sum(freq_vec[!keep_no_transient])
sum(freq_vec[!keep_no_transient]) / sum(freq_vec)

op_m2 <- openp(bunting, dfreq = TRUE, keep = keep_no_transient)
op_m2$model.fit

dev_diff <- op_m1$model.fit["deviance"] - op_m2$model.fit["deviance"]
df_diff  <- op_m1$model.fit["df"] - op_m2$model.fit["df"]

plot(op_m2)


comparaison <- data.frame(
  modele = factor(c("m1 (avec transitoires)", "m2 (sans transitoires)"),
                   levels = c("m1 (avec transitoires)", "m2 (sans transitoires)")),
  deviance = c(219.41, 125.18)
)

ggplot(comparaison, aes(x = modele, y = deviance, fill = modele)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = deviance), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("#a8c9b8", "#1b4332")) +
  labs(title = "Retirer les transitoires ameliore nettement l'ajustement",
       x = NULL, y = "Deviance residuelle") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")


# 6. Retrait des residus de Pearson extremes
pearson_resid <- residuals(op_m2$glm, type = "pearson")
keep_pearson_ok <- pearson_resid < 4

idx_no_transient <- which(keep_no_transient)
idx_final <- idx_no_transient[keep_pearson_ok]

keep_final <- rep(FALSE, nrow(bunting))
keep_final[idx_final] <- TRUE

op_m3 <- openp(bunting, dfreq = TRUE, keep = keep_final)
op_m3$model.fit

survie_comparaison <- cbind(
  m2_avec_residus_extremes = op_m2$survivals[, 1],
  m3_sans_residus_extremes = op_m3$survivals[, 1]
)
survie_comparaison


survie_df <- data.frame(
  transition = rep(c("2->3","3->4","4->5","5->6","6->7"), 2),
  modele = rep(c("m2", "m3"), each = 5),
  estimate = c(0.4851117, 0.6742944, 0.7287239, 0.5176471, 0.5559809,
               0.4815109, 0.6188964, 0.7013263, 0.5012495, 0.5512693)
)

ggplot(survie_df, aes(x = transition, y = estimate, color = modele, group = modele)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("m2" = "#7fb69e", "m3" = "#1b4332")) +
  labs(title = "Survie stable malgre le retrait des residus extremes",
       x = "Transition", y = "Survie estimee (phi)", color = NULL) +
  theme_minimal(base_size = 12)



# 7. Test d'homogeneite de la survie (moindres carres generalises + chi2)
phi_vec <- op_m3$survivals[2:6, 1]
cov_phi <- op_m3$cov[paste0("phi ", 2:6), paste0("phi ", 2:6)]

k <- length(phi_vec)
ones <- rep(1, k)
siginv <- solve(cov_phi)

phi_commun <- as.numeric((t(ones) %*% siginv %*% phi_vec) /
                          (t(ones) %*% siginv %*% ones))
se_phi_commun <- as.numeric(1 / sqrt(t(ones) %*% siginv %*% ones))

chi2_stat <- as.numeric(t(phi_vec - phi_commun) %*% siginv %*% (phi_vec - phi_commun))
chi2_df <- k - 1
chi2_pval <- 1 - pchisq(chi2_stat, df = chi2_df)

data.frame(
  phi_commun = phi_commun, se = se_phi_commun,
  chi2 = chi2_stat, df = chi2_df, pvalue = chi2_pval
)


se_vec <- sqrt(diag(cov_phi))

forest_df <- data.frame(
  transition = factor(c("2->3","3->4","4->5","5->6","6->7"),
                       levels = c("2->3","3->4","4->5","5->6","6->7")),
  estimate = phi_vec,
  se = se_vec
)

ggplot(forest_df, aes(x = estimate, y = transition)) +
  geom_vline(xintercept = phi_commun, color = "#1b4332", linewidth = 0.8) +
  annotate("rect",
           xmin = phi_commun - 1.96*se_phi_commun,
           xmax = phi_commun + 1.96*se_phi_commun,
           ymin = -Inf, ymax = Inf, fill = "#1b4332", alpha = 0.1) +
  geom_pointrange(aes(xmin = estimate - 1.96*se, xmax = estimate + 1.96*se),
                   color = "#2e8b57", size = 0.6) +
  labs(title = "Survie constante non rejetee (p = 0,77)",
       x = "Survie estimee (phi)", y = "Transition") +
  theme_minimal(base_size = 12)


# 8. Question complementaire : taux de renouvellement de la population

data(bunting)

op_m1 <- openp(bunting, dfreq = TRUE)

B_vec <- op_m1$birth[2:6, 1]
N_vec <- op_m1$N[3:7, 1]

turnover_df <- data.frame(
  transition = factor(c("2->3", "3->4", "4->5", "5->6", "6->7"),
                       levels = c("2->3", "3->4", "4->5", "5->6", "6->7")),
  B = B_vec,
  N = N_vec
)
turnover_df$taux_renouvellement <- turnover_df$B / turnover_df$N

turnover_df
mean(turnover_df$taux_renouvellement)