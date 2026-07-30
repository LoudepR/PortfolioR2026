## 1. Installation et chargement (Keras necessite un backend TensorFlow,
##    a installer une seule fois via install_keras() puis redemarrer la session)
library(keras)
library(tensorflow)

## 2. Chemin du dataset (structure : un sous-dossier par classe)
base_dir <- "chemin/vers/dataset/"
classes <- list.files(base_dir)
cat("Nombre de classes detectees :", length(classes), "\n")

## 3. Generateurs d'images (augmentation + split train/validation)
train_datagen <- image_data_generator(
  rescale = 1/255,
  rotation_range = 20,
  width_shift_range = 0.1,
  height_shift_range = 0.1,
  zoom_range = 0.1,
  horizontal_flip = TRUE,
  validation_split = 0.2
)

train_generator <- flow_images_from_directory(
  base_dir, train_datagen,
  target_size = c(124, 124), batch_size = 32,
  class_mode = "categorical", color_mode = "rgb",
  subset = "training"
)

val_generator <- flow_images_from_directory(
  base_dir, train_datagen,
  target_size = c(124, 124), batch_size = 32,
  class_mode = "categorical", color_mode = "rgb",
  subset = "validation"
)

## 4. Architecture du modele (CNN a 3 couches de convolution)
model <- keras_model_sequential() %>%
  layer_conv_2d(filters = 16, kernel_size = 3, activation = "relu",
                input_shape = c(124, 124, 3)) %>%
  layer_max_pooling_2d(pool_size = 2) %>%
  layer_conv_2d(filters = 32, kernel_size = 3, activation = "relu") %>%
  layer_max_pooling_2d(pool_size = 2) %>%
  layer_conv_2d(filters = 64, kernel_size = 3, activation = "relu") %>%
  layer_max_pooling_2d(pool_size = 2) %>%
  layer_flatten() %>%
  layer_dense(units = 128, activation = "relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = train_generator$num_classes, activation = "softmax")

model %>% compile(
  optimizer = "adam",
  loss = "categorical_crossentropy",
  metrics = "accuracy"
)

## 5. Entrainement (2 sessions de 8 epochs, sous-echantillonnage a 500
##    steps/epoch pour rester compatible avec un entrainement CPU)
history_1 <- model %>% fit(
  train_generator,
  steps_per_epoch = 500,
  epochs = 8,
  validation_data = val_generator,
  validation_steps = 100,
  callbacks = list(
    callback_early_stopping(monitor = "val_loss", patience = 3, restore_best_weights = TRUE)
  )
)
save_model_hdf5(model, "modele_plantvillage_v1.h5")

history_2 <- model %>% fit(
  train_generator,
  steps_per_epoch = 500,
  epochs = 8,
  validation_data = val_generator,
  validation_steps = 100,
  callbacks = list(
    callback_early_stopping(monitor = "val_loss", patience = 3, restore_best_weights = TRUE)
  )
)
plot(history_2)
save_model_hdf5(model, "modele_plantvillage_v2.h5")

## 6. Evaluation detaillee
library(caret)
library(dplyr)
library(ggplot2)
library(viridis)

val_datagen_eval <- image_data_generator(rescale = 1/255, validation_split = 0.2)

val_generator_eval <- flow_images_from_directory(
  base_dir, val_datagen_eval,
  target_size = c(124, 124), batch_size = 32,
  class_mode = "categorical", color_mode = "rgb",
  subset = "validation",
  shuffle = FALSE   # necessaire pour aligner predictions et verites
)

n_steps <- ceiling(val_generator_eval$n / val_generator_eval$batch_size)
predictions_proba <- model %>% predict(val_generator_eval, steps = n_steps)
predictions_classe <- apply(predictions_proba, 1, which.max) - 1
vraies_classes <- val_generator_eval$classes
noms_classes <- names(val_generator_eval$class_indices)[
  order(unlist(val_generator_eval$class_indices))
]

n_reel <- length(vraies_classes)
predictions_classe <- predictions_classe[1:n_reel]

## 7. Matrice de confusion
conf_mat <- confusionMatrix(
  factor(predictions_classe, levels = 0:(length(noms_classes) - 1), labels = noms_classes),
  factor(vraies_classes, levels = 0:(length(noms_classes) - 1), labels = noms_classes)
)
print(conf_mat$overall)

conf_df <- as.data.frame(conf_mat$table)
names(conf_df) <- c("Predit", "Vrai", "Freq")

ggplot(conf_df, aes(x = Vrai, y = Predit, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis(option = "magma", trans = "log1p", name = "Nb images") +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Matrice de confusion - classification des maladies foliaires",
       subtitle = paste("Precision globale :", round(conf_mat$overall["Accuracy"] * 100, 1), "%"),
       x = "Classe reelle", y = "Classe predite")

## 8. Precision par classe
precision_par_classe <- data.frame(
  classe = rownames(conf_mat$byClass),
  sensibilite = conf_mat$byClass[, "Sensitivity"],
  precision = conf_mat$byClass[, "Pos Pred Value"],
  n_images = conf_mat$byClass[, "Detection Rate"] * n_reel
) %>%
  arrange(sensibilite)

ggplot(head(precision_par_classe, 15), aes(x = reorder(classe, sensibilite), y = sensibilite)) +
  geom_col(fill = "#E41A1C", alpha = 0.8) +
  coord_flip() +
  labs(title = "15 classes les moins bien reconnues",
       x = NULL, y = "Sensibilite (rappel)") +
  theme_minimal(base_size = 11)

## 9. Paires de classes les plus confondues
confusions <- conf_df %>%
  filter(Vrai != Predit, Freq > 0) %>%
  arrange(desc(Freq))

ggplot(head(confusions, 15),
       aes(x = reorder(paste(Vrai, "->", Predit), Freq), y = Freq)) +
  geom_col(fill = "#377EB8", alpha = 0.8) +
  coord_flip() +
  labs(title = "15 confusions les plus frequentes entre classes",
       x = NULL, y = "Nombre d'images mal classees") +
  theme_minimal(base_size = 10)

## 10. Exemples concrets d'images mal classees
indices_erreurs <- which(predictions_classe != vraies_classes)
chemins_fichiers <- val_generator_eval$filepaths

set.seed(42)
indices_tires <- sample(indices_erreurs, 10)

exemples_erreurs <- data.frame(
  fichier = basename(chemins_fichiers[indices_tires]),
  vraie_classe = noms_classes[vraies_classes[indices_tires] + 1],
  classe_predite = noms_classes[predictions_classe[indices_tires] + 1]
)
print(exemples_erreurs)

## 11. Prediction sur une nouvelle image (fonction reutilisable)
predire_maladie <- function(chemin_image, modele, noms_classes, taille = c(124, 124)) {

  img <- image_load(chemin_image, target_size = taille)
  img_array <- image_to_array(img)

  # Meme normalisation que celle appliquee pendant l'entrainement
  # (rescale = 1/255 dans image_data_generator)
  img_array <- img_array / 255
  img_array <- array_reshape(img_array, c(1, dim(img_array)))

  proba <- modele %>% predict(img_array)

  classe_predite_id <- which.max(proba) - 1  # -1 : Keras est 0-based
  classe_predite_nom <- noms_classes[classe_predite_id + 1]
  confiance <- max(proba) * 100

  top3_idx <- order(proba, decreasing = TRUE)[1:3]
  top3 <- data.frame(
    classe = noms_classes[top3_idx],
    probabilite_pct = round(proba[top3_idx] * 100, 1)
  )

  cat("Image analysee :", basename(chemin_image), "\n")
  cat("Prediction principale :", classe_predite_nom,
      "(confiance :", round(confiance, 1), "%)\n\n")
  cat("Top 3 des predictions :\n")
  print(top3)

  return(invisible(list(
    classe_predite = classe_predite_nom,
    confiance = confiance,
    top3 = top3
  )))
}

## 12. Utilisation sur une image externe au jeu d'entrainement
resultat <- predire_maladie(
  chemin_image = "chemin/vers/image_test.jpg",
  modele = model,
  noms_classes = noms_classes
)