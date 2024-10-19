# Instalacja pakietów
install.packages("ScottKnottESD")

# Ładowanie pakietu
library(ScottKnottESD)
required_packages <- c(
    "ScottKnottESD", "readr", "ggplot2", "gridExtra", "tidyverse",
    "psych", "FSA", "lattice", "coin", "PMCMRplus", "rcompanion", "DescTools"
)

library(caret)
library(dplyr)
library(ggplot2)
library(gridExtra) # Do łączenia wielu wykresów

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

if (!require(caret)) {
    install.packages("caret")
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

regressors <- unique(data$Regressor)
datasets <- unique(data$Dataset)

plot_list <- list()
plot_list_overall <- list()
kendall_data <- data.frame()
plot_list_preliminary <- list()

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

for (dat in datasets) {
    combined_data_dataset <- data[data$Dataset == dat, ]

    # Zastosowanie mapowania
    combined_data_dataset <- combined_data_dataset %>%
        mutate(Validation_method = factor(Validation_method, levels = names(validation_mapping), labels = validation_mapping))

    combined_data_dataset <- combined_data_dataset %>%
        mutate(Hyper_param_method = factor(Hyper_param_method, levels = names(optimization_mapping), labels = optimization_mapping))

    for (reg in regressors) {
        combined_data <- combined_data_dataset[combined_data_dataset$Regressor == reg, ]

        # Tworzenie unikalnego identyfikatora dla każdej kombinacji
        combined_data$ComboID <- with(combined_data, gsub(" ", ".", paste(Hyper_param_method, Validation_method, n_splits, sep = ".")))

        filtered_data <- combined_data

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

        # Usunięcie wybranych kolumn z tabeli data_wide
        names(data_wide) <- gsub("\\.(HO|B)\\.3$", ".\\1", names(data_wide))

        # Przeprowadzenie testu Scott-Knott ESD
        sk_results <- sk_esd(as.matrix(data_wide), v = "np")

        max_rank <- max(sk_results$groups)

        sk_ranks <- data.frame(
            model = names(sk_results$groups),
            rank = paste0(max_rank + 1 - sk_results$groups)
        )

        # Ekstrakcja nazw metod walidacji i optymalizacji z nazwy modelu
        sk_ranks$Validation_method <- sapply(strsplit(sk_ranks$model, "\\."), function(x) x[2])
        sk_ranks$Hyper_param_method <- sapply(strsplit(sk_ranks$model, "\\."), function(x) x[1])
        sk_ranks$n_splits <- sapply(strsplit(sk_ranks$model, "\\."), function(x) x[3])

        # Łączenie metody walidacji z n_splits
        sk_ranks$Validation_method_with_num_splits <- paste(sk_ranks$Validation_method, sk_ranks$n_splits, sep = ".")

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
        # Obliczanie średniego MAE dla każdego Var2
        mean_mae <- aggregate(value ~ Var2, data = plot_data, FUN = mean)

        # Sortowanie plot_data według średniego MAE
        plot_data$Var2 <- factor(plot_data$Var2, levels = mean_mae$Var2[order(mean_mae$value)])

        rank_plot <- ggplot(data = plot_data, aes(x = Var2, y = value, fill = Var2)) +
            geom_boxplot(outlier.shape = 16, outlier.size = 1, notch = FALSE, alpha = 0.7, fatten = 1) +
            stat_summary(
                fun = mean, geom = "errorbar", aes(ymin = ..y.., ymax = ..y..),
                width = 0.7, color = "pink", size = 0.25
            ) + # Dodanie średniej jako różowe linie
            ylim(global_min, global_max) +
            facet_grid(~rank, scales = "free_x", space = "free_x") +
            scale_fill_manual(values = methods_colors) +
            ggtitle(reg) +
            xlab("") + # Usunięcie nazwy osi X
            ylab("MAE") + # Zmiana nazwy osi Y na MAE
            theme_bw() +
            theme(
                plot.title = element_text(size = 14, hjust = 0.5), # Zwiększono rozmiar tytułu
                axis.text.x = element_text(angle = 45, hjust = 1, size = 10.5), # Zwiększono rozmiar tekstu na osi X
                axis.text.y = element_text(size = 10.5), # Zwiększono rozmiar tekstu na osi Y
                text = element_text(size = 10), # Zwiększono ogólny rozmiar tekstu
                legend.position = "none",
                panel.spacing = unit(1, "lines"),
                panel.border = element_rect(colour = "grey50", fill = NA, size = 1),
                axis.title.x = element_text(size = 12), # Zwiększono rozmiar tytułu osi X
                axis.title.y = element_text(size = 12) # Zwiększono rozmiar tytułu osi Y
            )
        plot_list[[reg]] <- rank_plot

        # Konwersja kolumny rank na typ numeryczny
        sk_ranks$rank <- as.numeric(sk_ranks$rank)

        # Inicjalizacja wartości min i max
        global_min <- Inf
        global_max <- -Inf

        # Aktualizacja wartości min i max
        local_min <- min((sk_ranks$rank))
        local_max <- max((sk_ranks$rank))

        if (local_min < global_min) {
            global_min <- local_min
        }

        if (local_max > global_max) {
            global_max <- local_max
        }
        print(combined_data)

        # Tworzenie wykresu pudełkowego
        boxplot_overall_rank <- ggplot(
            data = sk_ranks %>%
                arrange(Validation_method, rank) %>%
                group_by(Validation_method) %>%
                mutate(mean_rank = mean(rank)),
            aes(x = reorder(Validation_method, mean_rank), y = rank, fill = Validation_method)
        ) +
            ggtitle(paste(dat, " - ", reg, sep = "")) +
            ylim(global_min, global_max) +
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
            scale_y_continuous(breaks = seq(min(sk_ranks$rank), max(sk_ranks$rank), by = 1), limits = c(global_min, global_max)) + # Zakładając, że 'rank' to nazwa kolumny
            theme(
                plot.title = element_text(size = 22, hjust = 0.5),
                axis.text.x = element_text(angle = 45, hjust = 1, size = 22), # Zwiększony rozmiar czcionki dla nazw metod
                axis.text.y = element_text(size = 16), # Zwiększony rozmiar czcionki dla wartości na osi Y
                axis.title.y = element_text(size = 22, vjust = 0.5), # Zwiększony rozmiar czcionki dla etykiety "Ranga"
                text = element_text(size = 18),
                legend.position = "none",
                panel.spacing = unit(1, "lines"),
                panel.border = element_rect(colour = "grey50", fill = NA, size = 1),
                axis.title.x = element_text(size = 22)
            )

        plot_list_overall[[paste("Boxplot", reg, dat)]] <- boxplot_overall_rank
    }

    # Aktualizacja wartości min i max
    min_val <- min(combined_data_dataset$MAE, na.rm = TRUE)
    max_val <- max(combined_data_dataset$MAE, na.rm = TRUE)
    breaks <- pretty(range(min_val, max_val), 1000)

    methods_colors <- hsv(h = seq(0.55, 0.75, length.out = length(unique(combined_data_dataset$Hyper_param_method))), s = 1, v = 1)
    names(methods_colors) <- levels(combined_data_dataset$Hyper_param_method)

    # Obliczanie średniej wartości MAE dla każdej metody walidacji
    mean_mae <- aggregate(MAE ~ Validation_method, data = combined_data_dataset, FUN = mean)

    # Sortowanie metod walidacji według średniego MAE
    combined_data_dataset$Validation_method <- factor(combined_data_dataset$Validation_method, levels = mean_mae$Validation_method[order(mean_mae$MAE)])

    # Tworzenie wykresu z posortowanymi metodami walidacji
    plot_preliminary <- ggplot(data = combined_data_dataset, aes(x = Validation_method, y = MAE, fill = Hyper_param_method)) +
        geom_boxplot(outlier.shape = 16, outlier.size = 1, notch = FALSE, alpha = 0.7, fatten = 1) +
        scale_y_continuous() +
        facet_grid(~Hyper_param_method, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = methods_colors) +
        ggtitle(dat) +
        xlab("Metody") +
        ylab("MAE") +
        theme_bw() +
        theme(
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            text = element_text(size = 11),
            legend.position = "none",
            panel.spacing = unit(1, "lines"),
            panel.border = element_rect(colour = "grey50", fill = NA, size = 1),
            axis.title.x = element_text(size = 11),
            axis.title.y = element_text(size = 11)
        )

    plot_list_preliminary[[paste("Boxplot-", dat)]] <- plot_preliminary

    if (length(plot_list) > 0) {
        combined_plot <- do.call(grid.arrange, c(plot_list, ncol = 1))
        image_name <- paste("combined_plots_MAE_", dat, ".png", sep = "")
        ggsave(image_name, combined_plot, width = 12, height = length(plot_list) * 3)
    }
}

if (length(plot_list_overall) > 0) {
    combined_plot_1 <- do.call(grid.arrange, c(plot_list_overall, ncol = 6))
    image_name_1 <- paste("combined_plots_val_without_splits_mean_rank.png", sep = "")
    ggsave(image_name_1, combined_plot_1, width = 22, height = length(plot_list) * 2.3, limitsize = FALSE)
}

if (length(plot_list_preliminary) > 0) {
    combined_plot_1 <- do.call(grid.arrange, c(plot_list_preliminary, ncol = 2))
    image_name_1 <- "preliminary_boxplots.png"
    ggsave(image_name_1, combined_plot_1, width = 8, height = length(plot_list_preliminary) * 1.25)
}
