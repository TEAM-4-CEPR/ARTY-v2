# Mettre cette ligne au tout début de ton script R, avant même les appels aux bibliothèques.
# Cela augmente la taille maximale des fichiers uploadés à 4000 Mo (4 GB).
options(shiny.maxRequestSize = 4000 * 1024^2)
options(shiny.fullStacktraces = TRUE) # Pour un débogage plus facile

library(shiny)
library(ggplot2) # Pour les boxplots améliorés
library(ggpubr)  # Pour stat_compare_means et les tests statistiques
library(plotly)
library(dplyr)
library(DESeq2) # Assure-toi que ce package est installé manuellement (BiocManager::install("DESeq2"))
library(data.table) # Pour fread, lecture rapide des CSV
library(DT) # Pour les tableaux interactifs
library(rlang) # Pour !!sym()
library(rstatix) # Pour shapiro_test, et les fonctions de test 2 à 2 comme wilcox_test, t_test etc.
library(tibble) # Pour column_to_rownames

ui <- shinyUI(fluidPage(
  
  titlePanel("ARTY_V2 - Analyse de Données Multi-Fonction"),
  tabsetPanel(
    
    tabPanel("Boxplot & Barrres",
             titlePanel("Exploration des Données et Visualisation"),
             sidebarLayout(
               sidebarPanel(
                 fileInput('file1', 'Choisissez un fichier CSV pour l\'analyse (Boxplot/Barre)',
                           accept=c('text/csv',
                                    'text/comma-separated-values,text/plain',
                                    '.csv')),
                 tags$br(),
                 checkboxInput('header', 'En-tête', TRUE),
                 radioButtons('sep', 'Séparateur',
                              c(Semicolon=';', Comma=',', Tab='\t'),
                              ';'),
                 radioButtons('quote', 'Guillemet',
                              c(None='', 'Double Quote'='"', 'Single Quote'="'"),
                              '"'),
                 uiOutput('variable_numerique_selector'), # Remplacement de xcol_ui
                 uiOutput('variable_facteur_selector'),   # Remplacement de ycol_ui
                 
                 # Nouvelle option pour le mode de comparaison
                 radioButtons("comparison_mode", "Mode de comparaison:",
                              choices = c("Comparaison globale" = "global",
                                          "Comparaisons 2 à 2 (personnalisé)" = "pairwise"),
                              selected = "global"),
                 
                 # UI pour sélectionner les groupes pour les comparaisons 2 à 2
                 uiOutput("pairwise_groups_selector"),
                 
                 checkboxInput("add_comparisons_boxplot", "Ajouter les tests statistiques sur le plot", FALSE),
                 
                 hr(),
                 h4("Information sur le test statistique choisi :"),
                 verbatimTextOutput("chosen_test_info"), # Affiche le test choisi
                 hr(),
                 h4("Prévisualisation des Données"),
                 DTOutput("contents") # Ajout pour la prévisualisation interactive
                 # Anciens contrôles d'exportation supprimés car Plotly n'est plus utilisé ici
               ),
               mainPanel(
                 plotOutput('Myplot_boxplot', height = "600px") # Changé de plotlyOutput à plotOutput
                 # Anciens résultats Shapiro-Wilk supprimés car intégrés dans le flux de décision du test
               )
             )
    ),
    tabPanel("PCA", titlePanel("Analyse en Composantes Principales (PCA)"),
             sidebarPanel(
               fileInput('file_pca', 'Fichier des données de comptage (CSV)',
                         accept=c('text/csv',
                                  'text/comma-separated-values,text/plain',
                                  '.csv')),
               fileInput('file_group_pca', 'Fichier des groupes pour annotation (CSV)',
                         accept=c('text/csv',
                                  'text/comma-separated-values,text/plain',
                                  '.csv')),
               uiOutput('group_col_pca_ui'),
               radioButtons(inputId = "plot_type" , label = "Sélectionner le type de graphique", choices = c("2D", "3D"))
             ),
             mainPanel(
               plotlyOutput('Myplot_pca')
             )),
    tabPanel("Heatmap", titlePanel("Heatmap"),
             sidebarPanel(
               fileInput('file_heatmap', 'Fichier des données de comptage (CSV)',
                         accept=c('text/csv',
                                  'text/comma-separated-values,text/plain',
                                  '.csv'))
             ),
             mainPanel(
               plotlyOutput('MyPlot_heatmap')
             )),
    tabPanel("Design & DE", titlePanel("Analyse d'Expression Différentielle (DESeq2)"),
             sidebarLayout(
               sidebarPanel(
                 fileInput('file_counts_de', 'Fichier des données de comptage (CSV, avec ID en première colonne)',
                           accept=c('text/csv',
                                    'text/comma-separated-values,text/plain',
                                    '.csv')),
                 fileInput('file_coldata_de', 'Fichier de données d\'échantillon (CSV, avec les groupes)',
                           accept=c('text/csv',
                                    'text/comma-separated-values,text/plain',
                                    '.csv')),
                 uiOutput('group_col_de_ui'),
                 downloadButton("downloadDesign", "Télécharger le Design (colData)"),
                 hr(),
                 actionButton("runDESeq2", "Lancer l'analyse DESeq2"),
                 downloadButton("downloadDEresults", "Télécharger les Résultats DE")
               ),
               mainPanel(
                 h4("Table de design (colData)"),
                 tableOutput('design_table'),
                 h4("Résultats de l'analyse d'expression différentielle"),
                 tableOutput("DE_results_table")
               )
             )),
    tabPanel("VolcanoPlot", titlePanel("Volcano Plot"),
             sidebarPanel(
               fileInput('file_de_results', 'Fichier de résultats DE (CSV)',
                         accept=c('text/csv',
                                  'text/comma-separated-values,text/plain',
                                  '.csv')),
               sliderInput("log2fc_threshold", "Seuil log2FoldChange:", min = 0, max = 3, value = 1, step = 0.1),
               sliderInput("padj_threshold", "Seuil p-value ajustée (-log10):", min = 0, max = 10, value = 1.3, step = 0.1)
             ),
             mainPanel(
               plotlyOutput('Myplot_volcano')
             ))
  )
))

server <- shinyServer(function(input, output, session) {
  
  # --- Reactive expressions pour la lecture des fichiers (centralisées) ---
  
  # Fonction utilitaire pour lire les CSV avec fread
  read_csv_robust <- reactive({
    function(file_input) {
      req(file_input)
      # fread gère le paramètre quote différemment, on le retire pour la simplicité
      df <- fread(file_input$datapath, header = input$header, sep = input$sep, data.table = FALSE)
      names(df) <- make.names(names(df), unique = TRUE) # Nettoyage des noms de colonnes
      return(df)
    }
  })
  
  # Données pour l'onglet "Boxplot & Barres"
  data_main <- reactive({
    read_csv_robust()(input$file1)
  })
  
  # Données pour l'onglet PCA
  data_pca <- reactive({
    req(input$file_pca)
    df <- read_csv_robust()(input$file_pca)
    # Assurer que la première colonne est bien utilisée comme noms de lignes pour la transposition
    if ("id" %in% names(df)) {
      row.names(df) <- df$id
      df$id <- NULL
    } else {
      warning("La colonne 'id' n'a pas été trouvée pour la PCA. Utilisation de la première colonne.")
      row.names(df) <- df[, 1]
      df <- df[, -1]
    }
    return(df)
  })
  
  data_group_pca <- reactive({
    req(input$file_group_pca)
    df <- read_csv_robust()(input$file_group_pca)
    # Assurer que la première colonne est bien utilisée comme noms de lignes
    if (!is.null(df) && ncol(df) > 0 && !is.null(names(df)[1])) {
      # Assurez-vous que la première colonne est un identifiant unique
      if (any(duplicated(df[,1]))) {
        validate("La première colonne du fichier de groupes PCA contient des identifiants d'échantillons dupliqués. Veuillez vous assurer qu'elle est unique.")
      }
      row.names(df) <- df[, 1]
      df <- df[, -1, drop = FALSE]
    }
    return(df)
  })
  
  
  # Données pour l'onglet Heatmap
  data_heatmap <- reactive({
    req(input$file_heatmap)
    df <- read_csv_robust()(input$file_heatmap)
    if ("id" %in% names(df)) {
      row.names(df) <- df$id
      df$id <- NULL
    } else {
      warning("La colonne 'id' n'a pas été trouvée pour la heatmap. Utilisation de la première colonne.")
      row.names(df) <- df[, 1]
      df <- df[, -1]
    }
    return(df)
  })
  
  # Données pour l'onglet DE (counts et colData)
  data_counts_de <- reactive({
    req(input$file_counts_de)
    df <- read_csv_robust()(input$file_counts_de)
    if ("id" %in% names(df)) {
      row.names(df) <- df$id
      df$id <- NULL
    } else {
      warning("La colonne 'id' n'a pas été trouvée pour les données de comptage DE. Utilisation de la première colonne.")
      row.names(df) <- df[, 1]
      df <- df[, -1]
    }
    df <- as.matrix(df)
    mode(df) <- "integer" # DESeq2 exige des entiers
    return(df)
  })
  
  data_coldata_de <- reactive({
    req(input$file_coldata_de)
    df <- read_csv_robust()(input$file_coldata_de)
    if (!is.null(df) && ncol(df) > 0) {
      # Assurez-vous que la première colonne est un identifiant unique pour les noms de lignes
      if (any(duplicated(df[,1]))) {
        validate("La première colonne du fichier d'échantillons DE contient des identifiants d'échantillons dupliqués. Veuillez vous assurer qu'elle est unique.")
      }
      row.names(df) <- df[, 1]
      df <- df[, -1, drop = FALSE]
    }
    return(df)
  })
  
  # Données pour l'onglet VolcanoPlot
  data_de_results <- reactive({
    req(input$file_de_results)
    df <- read_csv_robust()(input$file_de_results)
    if (!"log2FoldChange" %in% names(df) || !"padj" %in% names(df)) {
      validate(
        "Le fichier de résultats DE doit contenir les colonnes 'log2FoldChange' et 'padj'."
      )
    }
    df$log2FoldChange <- as.numeric(df$log2FoldChange)
    df$padj <- as.numeric(df$padj)
    return(df)
  })
  
  
  # --- UI dynamiques pour les sélecteurs de colonnes (Boxplot & Barres) ---
  
  # UI pour le sélecteur de variable numérique
  output$variable_numerique_selector <- renderUI({
    df <- data_main()
    if (!is.null(df)) {
      numeric_cols <- names(df)[sapply(df, is.numeric)]
      selectInput("variable_numerique", "Sélectionnez une variable numérique:", choices = numeric_cols)
    }
  })
  
  # UI pour le sélecteur de variable facteur
  output$variable_facteur_selector <- renderUI({
    df <- data_main()
    if (!is.null(df)) {
      # On peut vouloir inclure les colonnes character/factor comme facteurs
      factor_cols <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
      selectInput("variable_facteur", "Sélectionnez une variable facteur (groupes):", choices = factor_cols)
    }
  })
  
  # UI pour le sélecteur de groupes pour comparaisons 2 à 2
  output$pairwise_groups_selector <- renderUI({
    req(input$variable_facteur, input$comparison_mode)
    df <- data_main()
    
    if (!is.null(df) && input$comparison_mode == "pairwise") {
      # S'assurer que la variable facteur est bien un facteur pour obtenir les niveaux
      factor_col <- as.factor(df[[input$variable_facteur]])
      
      # Obtenir les niveaux uniques et les compter
      levels_factor <- levels(factor_col)
      
      # Créer toutes les paires possibles pour la sélection
      all_possible_pairs <- list()
      if (length(levels_factor) >= 2) {
        for (i in 1:(length(levels_factor) - 1)) {
          for (j in (i + 1):length(levels_factor)) {
            pair_name <- paste(levels_factor[i], "vs", levels_factor[j])
            all_possible_pairs[[pair_name]] <- c(levels_factor[i], levels_factor[j])
          }
        }
      }
      
      if (length(all_possible_pairs) == 0) {
        return(p("Pas assez de groupes (au moins 2) ou de données pour les comparaisons 2 à 2."))
      }
      
      selectInput("selected_pairwise_groups", 
                  "Sélectionnez les comparaisons paires :", 
                  choices = names(all_possible_pairs), 
                  multiple = TRUE,
                  selectize = TRUE) # selectize = TRUE pour une meilleure UX
    }
  })
  
  output$group_col_pca_ui <- renderUI({
    req(data_group_pca())
    selectInput('group_pca_colname', 'Nom de la colonne de groupe pour PCA:', choices = names(data_group_pca()))
  })
  
  output$group_col_de_ui <- renderUI({
    req(data_coldata_de())
    selectInput('group_de_colname', 'Nom de la colonne de groupe pour DE:', choices = names(data_coldata_de()))
  })
  
  
  # --- Logique de l'onglet "Boxplot & Barres" ---
  
  # Prévisualisation des données du Boxplot/Barre
  output$contents <- renderDT({
    df <- data_main()
    req(df)
    datatable(head(df, 10))
  })
  
  # Fonction pour exécuter le test de Shapiro-Wilk par groupe
  run_shapiro_test <- function(data, variable_numerique, variable_facteur) {
    groups <- unique(data[[variable_facteur]])
    shapiro_p_values <- list()
    
    for (group in groups) {
      subset_data <- data %>% filter(!!sym(variable_facteur) == group)
      numeric_values <- subset_data[[variable_numerique]]
      
      # Shapiro-Wilk nécessite au moins 3 échantillons
      if (length(numeric_values) >= 3) {
        shapiro_test_result <- tryCatch({
          shapiro.test(numeric_values)
        }, error = function(e) {
          warning(paste("Erreur Shapiro pour le groupe", group, ":", e$message))
          list(p.value = NA)
        })
        shapiro_p_values[[as.character(group)]] <- shapiro_test_result$p.value
      } else {
        shapiro_p_values[[as.character(group)]] <- NA # Pas assez de données pour le test
      }
    }
    return(shapiro_p_values)
  }
  
  # Reactive pour stocker l'information sur le test choisi
  test_info <- reactiveVal(NULL)
  
  # Boxplot avec comparaisons (p-values affichées sur le plot)
  output$Myplot_boxplot <- renderPlot({ # Changé en renderPlot
    data_boxplot <- data_main()
    req(data_boxplot, input$variable_numerique, input$variable_facteur)
    
    # Convertir en facteur
    data_boxplot[[input$variable_facteur]] <- as.factor(data_boxplot[[input$variable_facteur]])
    
    # Exécuter les tests de normalité de Shapiro-Wilk par groupe
    shapiro_results <- run_shapiro_test(data_boxplot, input$variable_numerique, input$variable_facteur)
    
    # Déterminer si les données sont globalement normales (tous les groupes avec assez d'observations)
    # Ignorer les NA si certains groupes ont moins de 3 observations pour le test de Shapiro
    normal_p_values <- unlist(shapiro_results[!is.na(shapiro_results)])
    is_normal <- all(normal_p_values > 0.05)
    
    # Déterminer la méthode de test et l'information à afficher
    n_groups <- length(unique(data_boxplot[[input$variable_facteur]]))
    
    global_method_to_use <- ""
    pairwise_method_to_use <- ""
    info_message <- ""
    
    if (n_groups < 2) {
      info_message <- "Moins de 2 groupes pour la variable facteur. Pas de test statistique possible."
      test_info(info_message)
      # Renvoyer le plot sans annotations
      p <- ggplot(data_boxplot, aes_string(x = input$variable_facteur, y = input$variable_numerique)) +
        geom_boxplot(aes_string(fill = input$variable_facteur)) +
        labs(title = paste("Boxplot de", input$variable_numerique, "par", input$variable_facteur),
             x = input$variable_facteur,
             y = input$variable_numerique) + # Correction ici
        theme_minimal()
      return(p)
    }
    
    if (is_normal) {
      if (n_groups == 2) {
        global_method_to_use <- "t.test" # Pour le message d'info global
        pairwise_method_to_use <- "t.test" # Pour les comparaisons 2 à 2
        info_message <- paste0("Données considérées normales (Shapiro p > 0.05 pour les groupes avec >= 3 obs). Test utilisé : ", global_method_to_use, ".")
      } else { # n_groups > 2
        global_method_to_use <- "anova"
        pairwise_method_to_use <- "t.test" # ANOVA global, mais t-test pour 2 à 2
        info_message <- paste0("Données considérées normales (Shapiro p > 0.05 pour les groupes avec >= 3 obs). Test global utilisé : ", global_method_to_use, ". Pour les comparaisons 2 à 2 : ", pairwise_method_to_use, ".")
      }
    } else {
      if (n_groups == 2) {
        global_method_to_use <- "wilcox.test"
        pairwise_method_to_use <- "wilcox.test"
        info_message <- paste0("Données considérées non normales (au moins un Shapiro p <= 0.05 ou < 3 obs). Test utilisé : ", global_method_to_use, ".")
      } else { # n_groups > 2
        global_method_to_use <- "kruskal.test"
        pairwise_method_to_use <- "wilcox.test" # Kruskal global, mais Wilcoxon pour 2 à 2
        info_message <- paste0("Données considérées non normales (au moins un Shapiro p <= 0.05 ou < 3 obs). Test global utilisé : ", global_method_to_use, ". Pour les comparaisons 2 à 2 : ", pairwise_method_to_use, ".")
      }
    }
    
    test_info(info_message)
    
    # Création du plot de base
    p <- ggplot(data_boxplot, aes_string(x = input$variable_facteur, y = input$variable_numerique)) +
      geom_boxplot(aes_string(fill = input$variable_facteur)) +
      labs(title = paste("Boxplot de", input$variable_numerique, "par", input$variable_facteur),
           x = input$variable_facteur,
           y = input$variable_numerique) +
      theme_minimal()
    
    if (input$add_comparisons_boxplot) {
      if (input$comparison_mode == "pairwise") {
        req(input$selected_pairwise_groups) # S'assurer que des groupes sont sélectionnés
        
        # Reconstruire la liste de comparaisons à partir des sélections de l'utilisateur
        comparisons_list_selected <- lapply(input$selected_pairwise_groups, function(pair_string) {
          strsplit(pair_string, " vs ")[[1]]
        })
        
        if (length(comparisons_list_selected) > 0) {
          # Ajustement dynamique de y.position pour éviter le chevauchement
          max_y <- max(data_boxplot[[input$variable_numerique]], na.rm = TRUE)
          # Utilise une proportion de l'étendue des données pour l'incrément
          y_range <- max_y - min(data_boxplot[[input$variable_numerique]], na.rm = TRUE)
          y_pos_increment <- y_range * 0.08 # Ajuste cet incrément si nécessaire
          
          # La première comparaison commence juste au-dessus du max
          current_y_pos <- max_y + y_pos_increment 
          
          for (k in seq_along(comparisons_list_selected)) {
            p <- p + stat_compare_means(comparisons = comparisons_list_selected[k], # Passer une seule comparaison à la fois
                                        label = "p.format",
                                        method = pairwise_method_to_use, # Utiliser la méthode 2 à 2
                                        hide.ns = TRUE,
                                        tip.length = 0.01, # Longueur des lignes de connexion plus courte
                                        bracket.size = 0.2, # Taille des crochets plus petite
                                        label.y = current_y_pos, # Positionnement dynamique
                                        size = 3.5) # Taille du texte plus petite
            
            current_y_pos <- current_y_pos + y_pos_increment # Incrémenter pour la prochaine comparaison
          }
          # Ajuster les limites de l'axe Y pour inclure toutes les annotations
          p <- p + expand_limits(y = current_y_pos)
          
        } else {
          showNotification("Aucune paire de groupes sélectionnée pour les comparaisons 2 à 2.", type = "warning")
        }
        
      } else if (input$comparison_mode == "global") {
        # Ajouter la comparaison globale avec p-value sur le plot
        p <- p + stat_compare_means(label.x.npc = "left", # Aligner le label à gauche
                                    label.y.npc = "top", # Placer le label en haut du plot
                                    method = global_method_to_use, # Utiliser la méthode globale
                                    label = "p.format",
                                    aes(label = paste0("Méthode : ", global_method_to_use, "\nP-value : ", ..p.format..)), # Reconstruire le label manuellement
                                    size = 4) # Taille du texte de la p-value
      }
    }
    
    p
  })
  
  # Affichage de l'information sur le test choisi
  output$chosen_test_info <- renderPrint({
    test_info()
  })
  
  
  # --- Logique de l'onglet PCA ---
  pca_results <- reactive({
    df_pca <- data_pca() # data_pca() est déjà la matrice transposée avec row.names = échantillon
    df_group_pca <- data_group_pca()
    req(input$group_pca_colname)
    
    mat_pca <- as.matrix(df_pca) # df_pca devrait déjà être la matrice avec échantillons en lignes
    mode(mat_pca) <- "numeric"
    mat_pca <- na.omit(mat_pca)
    
    if (ncol(mat_pca) < 2) {
      validate("Pas assez de colonnes numériques pour la PCA après le nettoyage des NA.")
    }
    
    # Vérification des dimensions après na.omit, au cas où toutes les lignes seraient vides
    if (nrow(mat_pca) < 2) {
      validate("Pas assez d'échantillons (lignes) valides pour la PCA après le nettoyage des NA.")
    }
    
    prin_comp <- prcomp(mat_pca, rank. = 3, scale. = TRUE)
    
    components <- as.data.frame(prin_comp[["x"]])
    
    # Correction: s'assurer que les noms de lignes de components et df_group_pca correspondent
    # Les rownames de components sont les noms des échantillons.
    # df_group_pca doit avoir les mêmes noms de lignes.
    if (!all(rownames(components) %in% rownames(df_group_pca))) {
      validate("Les identifiants des échantillons dans le fichier de données PCA ne correspondent pas aux identifiants dans le fichier de groupes pour PCA. Veuillez vérifier les noms des échantillons.")
    }
    
    # Filtrer df_group_pca pour qu'il ait les mêmes échantillons et le même ordre que components
    groups_col <- df_group_pca[rownames(components), input$group_pca_colname, drop = FALSE]
    components$groups <- as.factor(groups_col[[input$group_pca_colname]]) # Extraire la colonne en tant que vecteur
    
    # Ces inversions sont parfois nécessaires selon l'implémentation, à garder si l'orientation te convient.
    components$PC2 <- -components$PC2
    components$PC3 <- -components$PC3
    
    variance_explained <- summary(prin_comp)[["importance"]]['Proportion of Variance', ]
    components$variance_explained <- variance_explained
    
    return(components)
  })
  
  output$Myplot_pca <- renderPlotly({
    pca_data <- pca_results()
    req(pca_data)
    
    tot_explained_variance_ratio <- 100 * sum(pca_data$variance_explained[1:min(3, length(pca_data$variance_explained))])
    
    if (input$plot_type == "2D") {
      fig <- plot_ly(pca_data, x = ~PC1, y = ~PC2, type = 'scatter', color = ~groups,
                     colors = c('#636EFA','#EF553B', '#00CC96', '#AB63FA', '#FF6692', '#B6E880', '#FF97FF', '#FECB52'),
                     mode = 'markers', text = rownames(pca_data), hoverinfo = 'text') %>%
        layout(
          title = paste0("PCA 2D (Variance Expliquée Totale: ", round(tot_explained_variance_ratio, 2), "%)"),
          legend=list(title=list(text=input$group_pca_colname)),
          plot_bgcolor='#e5ecf6',
          xaxis = list(title = paste0("PC1 (", round(pca_data$variance_explained['PC1'] * 100, 2), "%)"), zerolinecolor = "#ffff", zerolinewidth = 2, gridcolor='#ffff'),
          yaxis = list(title = paste0("PC2 (", round(pca_data$variance_explained['PC2'] * 100, 2), "%)"), zerolinecolor = "#ffff", zerolinewidth = 2, gridcolor='#ffff')
        )
      fig
    } else if (input$plot_type == "3D") {
      fig <- plot_ly(pca_data, x = ~PC1, y = ~PC2, z = ~PC3, color = ~groups,
                     colors = c('#636EFA','#EF553B', '#00CC96', '#AB63FA', '#FF6692', '#B6E880', '#FF97FF', '#FECB52')) %>% add_markers(size = 8) %>%
        layout(
          title = paste0("PCA 3D (Variance Expliquée Totale: ", round(tot_explained_variance_ratio, 2), "%)"),
          scene = list(bgcolor = "#e5ecf6",
                       xaxis = list(title = paste0("PC1 (", round(pca_data$variance_explained['PC1'] * 100, 2), "%)")),
                       yaxis = list(title = paste0("PC2 (", round(pca_data$variance_explained['PC2'] * 100, 2), "%)")),
                       zaxis = list(title = paste0("PC3 (", round(pca_data$variance_explained['PC3'] * 100, 2), "%)")))
        )
      fig
    }
  })
  
  # --- Logique de l'onglet Heatmap ---
  heatmap_data_matrix <- reactive({
    df_heatmap <- data_heatmap() # data_heatmap() est déjà la matrice avec noms de lignes = gènes
    mat <- as.matrix(df_heatmap)
    mode(mat) <- "numeric"
    mat[is.na(mat)] <- 0 # Gérer les NA si tu veux qu'ils soient 0 dans la heatmap
    return(mat)
  })
  
  output$MyPlot_heatmap <- renderPlotly({
    mat <- heatmap_data_matrix()
    req(mat)
    
    if (nrow(mat) < 2 || ncol(mat) < 2) {
      validate("Pas assez de données pour générer une Heatmap. Nécessite au moins 2 lignes et 2 colonnes.")
    }
    
    p <- plot_ly(x = colnames(mat), y = rownames(mat), z = mat, type = "heatmap") %>%
      layout(margin = list(l = 120),
             xaxis = list(side = "top"),
             yaxis = list(autorange = "reversed")
      )
    p
  })
  
  # --- Logique de l'onglet Design & DE ---
  coldata_for_deseq <- reactive({
    df_coldata <- data_coldata_de()
    req(df_coldata, input$group_de_colname)
    
    df_coldata[[input$group_de_colname]] <- as.factor(df_coldata[[input$group_de_colname]])
    return(df_coldata)
  })
  
  output$design_table <- renderTable({
    coldata_for_deseq()
  }, rownames = TRUE)
  
  deseq_results <- eventReactive(input$runDESeq2, {
    counts_data <- data_counts_de()
    coldata_data <- coldata_for_deseq()
    
    if (is.null(counts_data) || ncol(counts_data) < 2) {
      validate("Veuillez charger un fichier de données de comptage valide (au moins 2 colonnes).")
    }
    if (is.null(coldata_data) || nrow(coldata_data) < 2) {
      validate("Veuillez charger un fichier de données d'échantillon valide (au moins 2 lignes).")
    }
    
    # Correction: S'assurer que les échantillons sont alignés et que les noms de colonnes des comptes
    # correspondent aux noms de lignes du colData.
    # Filtrer colData pour ne garder que les échantillons présents dans counts_data et les ordonner.
    if (!all(colnames(counts_data) %in% rownames(coldata_data))) {
      stop("Les noms de colonnes des données de comptage ne correspondent pas aux noms de lignes du fichier d'échantillons (colData).")
    }
    coldata_data <- coldata_data[colnames(counts_data), , drop = FALSE]
    
    
    if (nlevels(coldata_data[[input$group_de_colname]]) < 2) {
      validate("La colonne de groupe doit avoir au moins deux niveaux distincts pour l'analyse d'expression différentielle.")
    }
    
    design_formula <- as.formula(paste("~", input$group_de_colname))
    
    withProgress(message = 'Lancement de DESeq2', value = 0, {
      dds <- DESeqDataSetFromMatrix(countData = counts_data,
                                    colData = coldata_data,
                                    design = design_formula)
      incProgress(0.3, detail = "Exécution de DESeq...")
      dds <- DESeq(dds)
      incProgress(0.7, detail = "Récupération des résultats...")
      res <- results(dds)
      incProgress(1, detail = "Terminé !")
    })
    return(as.data.frame(res))
  })
  
  output$DE_results_table <- renderTable({
    deseq_results()
  }, rownames = TRUE)
  
  output$downloadDesign <- downloadHandler(
    filename = function() {
      paste("Design_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(coldata_for_deseq(), file, row.names = TRUE)
    }
  )
  
  output$downloadDEresults <- downloadHandler(
    filename = function() {
      paste("DESeq2_Results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(deseq_results(), file, row.names = TRUE)
    }
  )
  
  # --- Logique de l'onglet VolcanoPlot ---
  volcano_data_processed <- reactive({
    df_volcano <- data_de_results()
    req(df_volcano, input$log2fc_threshold, input$padj_threshold)
    
    df_volcano <- df_volcano[!is.na(df_volcano$log2FoldChange) & !is.na(df_volcano$padj), ]
    
    # Ajouter une petite constante aux padj = 0 pour éviter -Inf
    df_volcano$padj_neglog10 <- -log10(df_volcano$padj + .Machine$double.eps)
    
    df_volcano <- mutate(df_volcano, condition = case_when(
      log2FoldChange > input$log2fc_threshold & padj_neglog10 >= input$padj_threshold ~ "Up-régulé (significatif)",
      log2FoldChange < -input$log2fc_threshold & padj_neglog10 >= input$padj_threshold ~ "Down-régulé (significatif)",
      TRUE ~ "Non significatif"
    ))
    return(df_volcano)
  })
  
  output$Myplot_volcano <- renderPlotly({
    volcano_data <- volcano_data_processed()
    req(volcano_data)
    
    p <- plot_ly(data = volcano_data, x = ~log2FoldChange, y = ~padj_neglog10,
                 text = ~rownames(volcano_data),
                 mode = "markers", color = ~condition,
                 colors = c("Up-régulé (significatif)" = "red", "Down-régulé (significatif)" = "blue", "Non significatif" = "gray")) %>%
      layout(title = paste0("Volcano Plot (Seuils: |log2FC| > ", input$log2fc_threshold, ", -log10(padj) > ", input$padj_threshold, ")"),
             xaxis = list(title = "log2(FoldChange)"),
             yaxis = list(title = "-log10(Adjusted p-value)"))
    p
  })
  
})

shinyApp(ui, server)