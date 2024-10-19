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

if (!require("ScottKnottESD")) {
  install.packages("ScottKnottESD", dependencies = TRUE)
}
library(ScottKnottESD)

if (!require(dplyr)) {
  install.packages("dplyr")
}
library(dplyr)

# Import libraries
library(ScottKnottESD)
library(readr)
library(ggplot2)
library(caret)

# Ustawienie ścieżki do folderu z plikami CSV
folder_path <- "C:\\Users\\julia\\PycharmProjects\\plotsAnalysis\\ml-tuning-08-16-2024-splits-10"
folder_path1 <- "C:\\Users\\julia\\PycharmProjects\\plotsAnalysis\\ml-tuning-08-09-2024-splits-5"
folder_path2 <- "C:\\Users\\julia\\PycharmProjects\\plotsAnalysis\\ml-tuning-08-07-2024-splits-3"

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

dataset <- unique(data$Dataset)
plot_list <- list()
kendall_results <- data.frame()
# plot_list_dataset <- list()

# Mapowanie metod walidacji
validation_mapping <- c(
  "RepeatedKFold" = "RKF",
  "KFold" = "KF",
  "RepeatedHoldOut" = "MC",
  "HoldOut" = "HO",
  "Bootstrap" = "B"
)

# Mapowanie metod optymalizacji
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
    filtered_data <- pom_data[pom_data$Regressor == reg, ]

    # Zastosowanie mapowania
    filtered_data <- filtered_data %>%
      mutate(Validation_method = factor(Validation_method, levels = names(validation_mapping), labels = validation_mapping))
    filtered_data <- filtered_data %>%
      mutate(Hyper_param_method = factor(Hyper_param_method, levels = names(optimization_mapping), labels = optimization_mapping))

    # Tworzenie unikalnego identyfikatora dla każdej kombinacji
    filtered_data$ComboID <- with(filtered_data, gsub(" ", ".", paste(Hyper_param_method, Validation_method, n_splits, sep = ".")))

    # Tworzenie kolumny 'Group' na podstawie kolumny kombinacji
    filtered_data$Group <- filtered_data$ComboID

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

    # Przygotowanie ramki danych do generowania wizualizacji
    plot_data <- melt(as.matrix(data_wide))
    plot_data <- merge(plot_data, sk_ranks, by.x = "Var2", by.y = "model")

    # Generowanie wizualizacji z niestandardowym zakresem osi Y i swobodnym rozstawem osi X
    methods_colors <- hsv(h = seq(0.55, 0.75, length.out = 75), s = 1, v = 1)

    names(methods_colors) <- levels(plot_data$Var2)

    vals <- data.frame(iteration_id = 1:10, MAE = as.matrix(data_wide)[1:10, ])
    vals$method <- reg

    kendall_data <- rbind(kendall_data, vals)

    plot_data$dataset <- dat
    plot_data$regressor <- reg

    combined_data <- rbind(combined_data, plot_data)
  }

  kendall_res <- data.frame()

  # Analiza Kendall W dla różnych klasyfikatorów
  for (regressor_name in unique(kendall_data$method)) {
    regressor_data <- kendall_data %>%
      filter(method == regressor_name) %>%
      select(-iteration_id, -method)

    # Sprawdzenie, czy dane dla regresora istnieją
    if (nrow(regressor_data) > 0) {
      # Obliczenie statystyki Kendall W
      KendallW.MAE <- KendallW(regressor_data, correct = TRUE, test = TRUE)

      # Utworzenie wiersza z wynikami
      result_row <- data.frame(
        Folder = dat,
        Regressor = regressor_name,
        KendallW_Statistic = KendallW.MAE$estimate,
        KendallW_P_Value = KendallW.MAE$p.value
      )

      kendall_res <- rbind(kendall_res, result_row)
    }
  }

  kendall_results <- rbind(kendall_results, kendall_res)

  write.csv(kendall_results, "kendallW_.csv")
}

# Konwersja kolumny rank na typ numeryczny
combined_data$rank <- as.numeric(combined_data$rank)

scaled_ranks <- combined_data %>%
  group_by(regressor, dataset) %>%
  mutate(scaled_rank = caret::preProcess(as.data.frame(rank), method = "range", rangeBounds = c(1, 100)) %>%
    predict(as.data.frame(rank)) %>%
    unlist()) %>% # Konwersja wyników na prosty wektor
  ungroup()

scaled_ranks <- scaled_ranks %>%
  arrange(Var2, scaled_rank) %>%
  group_by(Var2) %>%
  mutate(mean_rank = mean(scaled_rank)) # Obliczanie średniej rangi dla każdej kombinacji


# Tworzenie wykresu pudełkowego
boxplot_overall <- ggplot(
  data = scaled_ranks,
  aes(x = reorder(Var2, mean_rank), y = scaled_rank, fill = Var2)
) +
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
  scale_y_continuous(breaks = c(1, 25, 50, 72, 100), limits = c(1, 100)) + # Skala od 1-100 i wybrane wartości na osi Y
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 15, vjust = 0.5),
    text = element_text(size = 10),
    legend.position = "none",
    panel.spacing = unit(1, "lines"),
    panel.border = element_rect(colour = "grey50", fill = NA, size = 1),
    axis.title.x = element_text(size = 10)
  )

print(boxplot_overall)

ggsave("boxplot.png", plot = boxplot_overall, width = 18, height = 5, units = "in")


scaled_ranks <- combined_data %>%
  group_by(regressor, dataset) %>%
  mutate(scaled_rank = caret::preProcess(as.data.frame(rank), method = "range", rangeBounds = c(1, 100)) %>%
    predict(as.data.frame(rank)) %>%
    unlist()) %>% # Konwersja wyników na prosty wektor
  ungroup()

pivot_table_regressor <- scaled_ranks %>%
  select(regressor, Var2, scaled_rank) %>%
  mutate(scaled_rank = as.numeric(scaled_rank)) %>%
  group_by(regressor, Var2) %>%
  summarize(mean_rank = mean(scaled_rank), .groups = "keep") %>%
  spread(key = Var2, value = mean_rank)

pivot_table_regressor <- pivot_table_regressor %>%
  column_to_rownames(var = "regressor")

average_row_regressor <- pivot_table_regressor %>%
  summarize_all(mean) %>%
  mutate(regressor = "Średnia") %>%
  column_to_rownames(var = "regressor")

pivot_table_regressor <- bind_rows(pivot_table_regressor, average_row_regressor)

pivot_table_dataset <- scaled_ranks %>%
  select(dataset, Var2, scaled_rank) %>%
  mutate(scaled_rank = as.numeric(scaled_rank)) %>%
  group_by(dataset, Var2) %>%
  summarize(mean_rank = mean(scaled_rank), .groups = "keep") %>%
  spread(key = Var2, value = mean_rank)

pivot_table_dataset <- pivot_table_dataset %>%
  column_to_rownames(var = "dataset")

average_row_dataset <- pivot_table_dataset %>%
  summarize_all(mean) %>%
  mutate(dataset = "Średnia") %>%
  column_to_rownames(var = "dataset")

pivot_table_dataset <- bind_rows(pivot_table_dataset, average_row_dataset)

pivot_table_regressor[] <- lapply(pivot_table_regressor, function(x) round(x, 2))
pivot_table_dataset[] <- lapply(pivot_table_dataset, function(x) round(x, 2))

transposed_df <- as.data.frame(t(pivot_table_regressor))
transposed_df$Średnia <- as.numeric(transposed_df$Średnia)
sorted_df <- transposed_df[order(transposed_df$Średnia), ]
pivot_table_regressor_sorted <- as.data.frame(t(sorted_df))

transposed_df2 <- as.data.frame(t(pivot_table_dataset))
transposed_df2$Średnia <- as.numeric(transposed_df$Średnia)
sorted_df2 <- transposed_df2[order(transposed_df$Średnia), ]
pivot_table_dataset_sorted <- as.data.frame(t(sorted_df2))

pivot_table_dataset_sorted <- pivot_table_dataset_sorted[, order(as.numeric(pivot_table_dataset_sorted["Średnia", ]), decreasing = FALSE)]
pivot_table_regressor_sorted <- pivot_table_regressor_sorted[, order(as.numeric(pivot_table_regressor_sorted["Średnia", ]), decreasing = FALSE)]

write.csv(pivot_table_regressor_sorted, "pivot_table_regressor_sorted.csv")
write.csv(pivot_table_dataset_sorted, "pivot_table_dataset_sorted.csv")
