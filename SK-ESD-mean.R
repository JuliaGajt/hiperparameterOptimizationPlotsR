# Instalacja pakietów
install.packages("ScottKnottESD")

# Ładowanie pakietu
library(ScottKnottESD)
required_packages <- c(
  "ScottKnottESD", "readr", "ggplot2", "gridExtra", "tidyverse",
  "psych", "FSA", "lattice", "coin", "PMCMRplus", "rcompanion", "DescTools"
)

if (!require(DescTools)) {
  install.packages("DescTools")
}
library(DescTools)

if (!require(tidyverse)) {
  install.packages("tidyverse")
}
library(tidyverse)

# Instalacja i ładowanie potrzebnych pakietów
if (!require("ScottKnottESD")) {
  install.packages("ScottKnottESD", dependencies = TRUE)
}
library(ScottKnottESD)

# Instalacja i załadowanie pakietu dplyr, jeśli jeszcze nie został załadowany
if (!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)

# Import libraries
library(ScottKnottESD)
library(readr)
library(ggplot2)

# Ustawienie ścieżki do folderu z plikami CSV
folder_path <- "C:\\Users\\julia\\PycharmProjects\\plotsAnalysis\\ml-tuning-08-16-2024-splits-10" # Zmień na odpowiednią ścieżkę
folder_path1 <- "C:\\Users\\julia\\PycharmProjects\\plotsAnalysis\\ml-tuning-08-09-2024-splits-5" # Zmień na odpowiednią ścieżkę
folder_path2 <- "C:\\Users\\julia\\PycharmProjects\\plotsAnalysis\\ml-tuning-08-07-2024-splits-3" # Zmień na odpowiednią ścieżkę

# Pobieranie listy plików CSV w folderze
file_list <- list.files(path = folder_path, pattern = "*.csv", full.names = TRUE)
file_list1 <- list.files(path = folder_path1, pattern = "*.csv", full.names = TRUE)
file_list2 <- list.files(path = folder_path2, pattern = "*.csv", full.names = TRUE)

# Wczytywanie i łączenie wszystkich plików CSV w jeden zbiór danych
data_list <- lapply(file_list, read.csv)
data_list1 <- lapply(file_list1, read.csv)
data_list2 <- lapply(file_list2, read.csv)
data <- do.call(rbind, c(data_list, data_list1, data_list2))
data <- data[data$Tuned == 1, ]
data <- data %>%
  filter(!(Validation_method %in% c("Bootstrap", "HoldOut") & n_splits %in% c(5, 10)))
combined_data <- data.frame()


# Filtrowanie danych dla metody optymalizacji 'grid search'

dataset <- unique(data$Dataset)
plot_list <- list()
kendall_results <- data.frame()
plot_list_dataset <- list()

# Mapowanie metod
validation_mapping <- c(
  "RepeatedKFold" = "RKF",
  "KFold" = "KF",
  "RepeatedHoldOut" = "MC",
  "HoldOut" = "HO",
  "Bootstrap" = "B"
)

# Mapowanie metod
optimization_mapping <- c(
  "grid search" = "GS",
  "bayes search" = "BS",
  "random search" = "RS",
  "SA" = "SA",
  "OPTUNA" = "TPE"
)


for (dat in dataset) {
  kendall_data <- data.frame()
  regressors <- unique(data$Regressor)
  pom_data <- data[data$Dataset == dat, ]

  for (reg in regressors) {
    combined_data_1 <- pom_data[pom_data$Regressor == reg, ]

    # Zastosowanie mapowania
    combined_data_1 <- combined_data_1 %>%
      mutate(Validation_method = factor(Validation_method, levels = names(validation_mapping), labels = validation_mapping))

    combined_data_1 <- combined_data_1 %>%
      mutate(Hyper_param_method = factor(Hyper_param_method, levels = names(optimization_mapping), labels = optimization_mapping))

    # Tworzenie unikalnego identyfikatora dla każdej kombinacji
    combined_data_1$ComboID <- with(combined_data_1, gsub(" ", ".", paste(Hyper_param_method, Validation_method, n_splits, sep = ".")))

    filtered_data <- combined_data_1

    # Tworzenie kolumny 'Group' na podstawie kolumny 'Dataset'
    filtered_data$Group <- filtered_data$ComboID

    # Usunięcie pustych grup
    filtered_data <- filtered_data[filtered_data$Group != "", ]

    # Upewnienie się, że mamy wystarczającą liczbę wierszy dla każdej grupy
    num_rows_per_group <- nrow(filtered_data) / length(unique(filtered_data$Group))

    # Przekształcenie danych do formatu szerokiego
    data_wide <- data.frame(matrix(ncol = length(unique(filtered_data$Group)), nrow = num_rows_per_group))
    names(data_wide) <- unique(filtered_data$Group)

    # Wypełnianie tabeli przekształconymi danymi
    for (group in unique(filtered_data$Group)) {
      data_wide[[group]] <- filtered_data$MAE[filtered_data$Group == group]
    }

    # Usunięcie wierszy zawierających NA w całym zestawie danych
    data_wide <- na.omit(data_wide)
    names(data_wide) <- gsub("\\.(HO|B)\\.3$", ".\\1", names(data_wide))

    # Przeprowadzenie testu Scott-Knott ESD
    sk_results <- sk_esd(as.matrix(data_wide), v = "np")

    max_rank <- max(sk_results$groups)

    sk_ranks <- data.frame(
      model = names(sk_results$groups),
      rank = paste0(max_rank + 1 - sk_results$groups)
    )

    # Inicjalizacja wartości min i max
    global_min <- Inf
    global_max <- -Inf

    # Aktualizacja wartości min i max
    local_min <- min(as.matrix(data_wide))
    local_max <- max(as.matrix(data_wide))

    if (local_min < global_min) {
      global_min <- local_min
    }

    if (local_max > global_max) {
      global_max <- local_max
    }

    plot_data <- melt(as.matrix(data_wide))
    plot_data <- merge(plot_data, sk_ranks, by.x = "Var2", by.y = "model")

    methods_colors <- hsv(h = seq(0.55, 0.75, length.out = 75), s = 1, v = 1)

    names(methods_colors) <- levels(plot_data$Var2)

    vals <- data.frame(iteration_id = 1:10, MAE = as.matrix(data_wide)[1:10, ])
    vals$method <- reg

    kendall_data <- rbind(kendall_data, vals)

    plot_data$dataset <- dat
    plot_data$regressor <- reg

    combined_data <- rbind(combined_data, plot_data)
  }

  combined_data$rank <- as.numeric(combined_data$rank)

  dataset_data <- combined_data %>%
    filter(dataset == dat) %>%
    arrange(Var2, rank) %>%
    group_by(Var2) %>%
    mutate(mean_rank = mean(rank)) # Obliczanie średniej rangi dla każdej kombinacji

  dataset_plot <- ggplot(
    dataset_data,
    aes(x = reorder(Var2, mean_rank), y = rank, fill = Var2)
  ) +
    ggtitle(dat) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, notch = FALSE, alpha = 0.7, fatten = 1) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      aes(group = 1),
      color = "blue",
      size = 0.4,
      linetype = "solid"
    ) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      position = position_dodge(width = 0.8),
      width = 0.7,
      fatten = 1,
      col = "black",
      linetype = "dashed"
    ) +
    scale_fill_manual(values = methods_colors) +
    ylab("Ranga") +
    xlab(NULL) +
    theme_bw() +
    scale_y_continuous(breaks = seq(min(dataset_data$rank), max(dataset_data$rank), by = 1)) + # Zakładając, że 'rank' to nazwa kolumny
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15), # Zwiększony rozmiar czcionki dla nazw metod
      axis.text.y = element_text(size = 15), # Zwiększony rozmiar czcionki dla wartości na osi Y
      axis.title.y = element_text(size = 15, vjust = 0.5), # Zwiększony rozmiar czcionki dla etykiety "Ranga"
      text = element_text(size = 10),
      legend.position = "none",
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(colour = "grey50", fill = NA, size = 1),
      axis.title.x = element_text(size = 10)
    )
  plot_list_dataset[[paste("Plot", dat)]] <- dataset_plot
}

if (length(plot_list_dataset) > 0) {
  library(gridExtra)
  combined_plot <- do.call(grid.arrange, c(plot_list_dataset, ncol = 1))
  ggsave("combined_plots_datasets.png", combined_plot, width = 20, height = length(plot_list_dataset) * 3)
}

plot_list_regressor <- list()

for (reg in regressors) {
  regressor_data <- combined_data %>%
    filter(regressor == reg) %>%
    arrange(Var2, rank) %>%
    group_by(Var2) %>%
    mutate(mean_rank = mean(rank)) # Obliczanie średniej rangi dla każdej kombinacji

  regressor_plot <- ggplot(
    regressor_data,
    aes(x = reorder(Var2, mean_rank), y = rank, fill = Var2)
  ) +
    ggtitle(reg) +
    geom_boxplot(outlier.shape = 16, outlier.size = 1, notch = FALSE, alpha = 0.7, fatten = 1) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      aes(group = 1),
      color = "blue", # Ustawienie koloru linii mediany na czerwony
      size = 0.4, # Grubość linii mediany
      linetype = "solid" # Styl linii (przerywana)
    ) +
    stat_summary(
      fun = mean,
      geom = "crossbar",
      position = position_dodge(width = 0.8),
      width = 0.7,
      fatten = 1,
      col = "black",
      linetype = "dashed"
    ) +
    scale_fill_manual(values = methods_colors) +
    ylab("Ranga") +
    xlab(NULL) +
    theme_bw() +
    scale_y_continuous(breaks = seq(min(dataset_data$rank), max(dataset_data$rank), by = 1)) + # Zakładając, że 'rank' to nazwa kolumny
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15), # Zwiększony rozmiar czcionki dla nazw metod
      axis.text.y = element_text(size = 15), # Zwiększony rozmiar czcionki dla wartości na osi Y
      axis.title.y = element_text(size = 15, vjust = 0.5), # Zwiększony rozmiar czcionki dla etykiety "Ranga"
      text = element_text(size = 10),
      legend.position = "none",
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(colour = "grey50", fill = NA, size = 1),
      axis.title.x = element_text(size = 10)
    )

  plot_list_regressor[[paste("Plot", reg)]] <- regressor_plot
}

if (length(plot_list_regressor) > 0) {
  library(gridExtra)
  combined_plot <- do.call(grid.arrange, c(plot_list_regressor, ncol = 1))
  ggsave("combined_plots_ml.png", combined_plot, width = 20, height = length(plot_list_regressor) * 3)
}
