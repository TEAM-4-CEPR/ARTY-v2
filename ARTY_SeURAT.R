# Mettre cette ligne au tout début de ton script R pour augmenter la limite de taille des fichiers
options(shiny.maxRequestSize = 4000 * 1024^2)

# IMPORTANT: Active les stack traces complètes pour le débogage
options(shiny.fullStacktraces = TRUE)

# Charger les librairies nécessaires
library(shiny)
library(Seurat)      # Le package principal pour l'analyse unicellulaire
library(ggplot2)     # Pour les visualisations statiques
library(plotly)      # Pour les visualisations interactives
library(dplyr)       # Pour la manipulation de données (comme mutate, filter)
library(data.table)  # For fread, for fast metadata reading if needed
library(Matrix)      # Essentiel pour lire les fichiers .mtx (matrices creuses)
library(cowplot)     # Utile pour combiner des plots
library(patchwork)   # Pour combiner des plots ggplot2
library(scDblFinder) # Pour la détection des doublets (Bioconductor)
library(RcppAnnoy)   # Dépendance de scDblFinder (peut être nécessaire de la charger explicitement)
library(SingleCellExperiment) # Pour la classe SingleCellExperiment, nécessaire pour scDblFinder
library(scater)      # Souvent utile avec SingleCellExperiment, pour manipulations, etc.
library(scCustomize) # Nouvelle librairie pour les visualisations d'expression génique
library(scales)      # Pour la fonction hue_pal() afin de générer des couleurs cohérentes
library(DT)          # Pour les tableaux interactifs
library(tibble)      # Pour rownames_to_column, si besoin
library(rlang)       # Pour !!sym() dans dplyr::arrange
library(tidyverse)   # Required by baranal and get_pathway (for piping, etc.)
library(ggalluvial)  # Required by baranal

# =========================================================================
# Fonctions Utilitaires (issues de script_utils.R)
# =========================================================================

# Note: La fonction prep_GO dépend d'OrgDb et d'un package d'annotation d'organisme (ex: org.Hs.eg.db)
# et n'est pas intégrée directement ici car elle nécessite des configurations supplémentaires
# (choix de l'organisme, installation du package OrgDb, etc.).

# Fonction baranal - Visualisation des proportions de clusters par échantillon
# Nécessite une colonne 'orig.ident' ou similaire dans les métadonnées pour identifier les échantillons.
baranal <- function(seurat_object) {
  # Assurez-vous que les packages nécessaires sont chargés, ou utilisez ::
  # library(tidyverse) # Déjà chargé globalement
  # library(scCustomize) # Déjà chargé globalement
  # library(ggalluvial) # Déjà chargé globalement
  
  # Cette fonction suppose l'existence de `orig.ident` dans les métadonnées de l'objet Seurat
  # pour catégoriser par échantillon. Si ce n'est pas le cas, vous devrez l'ajouter
  # ou modifier cette fonction pour utiliser une autre colonne d'échantillon.
  
  cluster_stats <- as.data.frame(scCustomize::Cluster_Stats_All_Samples(seurat_object = seurat_object))
  
  # Correction pour s'assurer que la ligne de total n'est pas incluse si elle est présente
  # La condition `row_number() <= n()-1` est basée sur le script original `Utils.R`
  # Elle est censée enlever la dernière ligne qui est souvent un total ou une ligne récapitulative
  if (nrow(cluster_stats) > 1 && any(grepl("%", names(cluster_stats)))) {
    # Si la dernière ligne contient des totaux (souvent en pourcentage ou des sommes)
    # ou si elle est clairement une ligne de résumé comme dans l'exemple original,
    # nous la filtrons. Sinon, gardons toutes les lignes de cluster.
    # Une approche plus robuste serait de ne pas compter les lignes de 'total' ou 'sum' si elles sont nommées.
    # Pour l'instant, je garde la logique du script Utils.R.
    cluster_stats <- cluster_stats %>% filter(row_number() <= n()-1)
  }
  
  cluster_stats <- cluster_stats %>% select(Cluster, ends_with("%"))
  tab <- cluster_stats %>% tidyr::gather(key = "keys", value = "values", -Cluster)
  
  # Créer une palette de couleurs basée sur le nombre de clusters
  n_clusters <- length(unique(tab$Cluster))
  # Utiliser la palette de couleurs stockée dans l'objet Seurat si elle existe
  # (Cette palette est générée et stockée par la réactive `cluster_colors` dans le serveur)
  if (!is.null(seurat_object@misc$cluster_colors)) {
    cluster_pal <- seurat_object@misc$cluster_colors
  } else {
    # Fallback si aucune palette n'est stockée (moins probable avec le serveur actuel)
    set.seed(42) # Pour la reproductibilité
    cluster_pal <- scales::hue_pal()(n_clusters)
  }
  
  p <- tab %>% ggplot(aes(y = values, x = keys, fill = Cluster)) +
    geom_flow(aes(alluvium = Cluster)) +
    geom_col(width = 0.5, alpha = 0.8, color = "black") + # geom_col pour les barres de base
    scale_fill_manual(values = cluster_pal) + # Utiliser la palette générée
    theme_minimal() +
    labs(x = "Échantillon", y = "Proportion de cellules", title = "Proportion de Clusters par Échantillon") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}


# Fonction get_pathway - FeaturePlot pour les gènes d'une voie
# Utilise une liste de gènes prédéfinie ou une liste fournie par l'utilisateur
get_pathway <- function(seurat_object, pathway_name, reduction = "umap") {
  # Définitions des listes de gènes de la fonction `pathway_gene_list` du script Utils
  # Note: Ces listes devraient idéalement être chargées ou définies globalement
  # ou dans une fonction qui les encapsule pour éviter la répétition.
  # Pour cet exemple, je les inclus directement ici, avec des noms capitalisés
  # pour correspondre aux symboles géniques humains courants.
  
  TCR_signaling <- c("CD3D", "CD3E", "CD3G", "CD247", "LCK", "ZAP70", "LAT", "LCP2",
                     "PTPRC", "CD28", "CTLA4", "ICOS", "PDCD1", "BTN3A1", "TRAC",
                     "TRBC1", "TRBC2", "TRGC1", "TRGC2", "TRDC")
  
  B_cell_pathway <- c("CD19", "CD79A", "CD79B", "MS4A1", "BLK", "BTK", "SYK", "PIK3CA",
                      "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3",
                      "AKT1", "AKT2", "AKT3", "MAPK1", "MAPK3", "MAPK8", "MAPK9",
                      "NFKBIA", "IKBKB", "CHUK", "JUN", "FOS", "MYC", "CR2", "CD22", "FCGR2B")
  
  Fas_pathway <- c("CASP3", "CASP6", "CASP7", "CASP8", "FADD", "FAS", "FASLG", "CFLAR",
                   "MAP2K4", "MAP3K1", "MAPK8", "PARP1", "TP53", "LMNA", "RB1", "DAXX")
  
  Il1R_signaling <- c("IL1A", "IL1B", "IL1R1", "IL1RAP", "IL1RN", "IRAK1", "IRAK2", "IRAK3",
                      "MYD88", "TRAF6", "CHUK", "IKBKB", "NFKBIA", "RELA", "RELB", "JUN",
                      "FOS", "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAP2K3", "MAP2K6")
  
  
  # Sélectionner la liste de gènes en fonction du nom du pathway
  gene_list <- switch(pathway_name,
                      "TCR_signaling" = TCR_signaling,
                      "B_cell_pathway" = B_cell_pathway,
                      "Fas_pathway" = Fas_pathway,
                      "Il1R_signaling" = Il1R_signaling,
                      stop("Pathway non reconnu ou non défini.")
  )
  
  # Filtrer les gènes qui sont réellement présents dans l'objet Seurat
  genes_present <- intersect(gene_list, rownames(seurat_object))
  
  if (length(genes_present) == 0) {
    stop(paste("Aucun gène du pathway '", pathway_name, "' n'a été trouvé dans l'objet Seurat."))
  }
  
  # Générer un plot pour chaque gène du pathway
  plots_list <- list()
  for (gene in genes_present) {
    # scCustomize::FeaturePlot_scCustom crée un ggplot2
    p <- scCustomize::FeaturePlot_scCustom(seurat_object = seurat_object, features = gene, reduction = reduction) +
      ggtitle(paste0(gene, " Expression (", pathway_name, ")"))
    plots_list[[gene]] <- p
  }
  
  # Combiner les plots si Patchwork est utilisé, ou juste le premier si un seul
  if (length(plots_list) == 1) {
    return(plots_list[[1]])
  } else if (length(plots_list) > 1) {
    # Utiliser patchwork pour combiner les plots
    combined_plot <- patchwork::wrap_plots(plots_list) +
      plot_annotation(title = paste("Expression des Gènes du Pathway:", pathway_name))
    return(combined_plot)
  } else {
    return(NULL) # Aucun plot si aucune gène présent
  }
}


# =========================================================================
# Définition de l'interface utilisateur (UI)
# =========================================================================
ui <- fluidPage(
  titlePanel("Analyse de Données scRNA-seq avec Seurat"),
  
  navbarPage(
    "Application scRNA-seq", # Titre de la barre de navigation
    id = "main_navbar", # ID pour gérer les onglets côté serveur
    
    # ---------------------------------------------------------------------
    # Onglet "Preprocessing"
    # ---------------------------------------------------------------------
    tabPanel("Preprocessing",
             icon = icon("cogs"), # Icône d'engrenage
             sidebarLayout(
               sidebarPanel(
                 h4("1. Charger les données"),
                 tabsetPanel(
                   id = "load_data_tabs",
                   tabPanel("Cell Ranger",
                            value = "cell_ranger_load_tab",
                            helpText("Veuillez uploader les trois fichiers générés par Cell Ranger (souvent dans le dossier 'filtered_feature_bc_matrix')."),
                            fileInput('mtx_file', 'Fichier matrix.mtx.gz ou .mtx',
                                      accept = c('.mtx.gz', '.mtx')),
                            fileInput('features_file', 'Fichier features.tsv.gz ou .tsv',
                                      accept = c('.tsv.gz', '.tsv', '.txt')),
                            fileInput('barcodes_file', 'Fichier barcodes.tsv.gz ou .tsv',
                                      accept = c('.tsv.gz', '.tsv', '.txt')),
                            actionButton("load_cellranger_data", "Charger les données Cell Ranger")
                   ),
                   tabPanel("Fichier Seurat (.rds)",
                            value = "rds_load_tab",
                            helpText("Charger un objet Seurat pré-existants au format .rds."),
                            fileInput('rds_file', 'Fichier Seurat (.rds)',
                                      accept = c('.rds')),
                            actionButton("load_rds_data", "Charger l'objet Seurat (.rds)")
                   )
                 ), # Fin tabsetPanel load_data_tabs
                 hr(),
                 
                 h4("2. Contrôle Qualité (QC)"),
                 tabsetPanel(
                   id = "qc_tabs", # Un ID pour gérer les onglets du QC
                   tabPanel("Filtrage Manuel",
                            value = "manual_qc_tab",
                            h5("Définir les seuils de filtrage:"),
                            sliderInput("min_features", "Min. Gènes par Cellule (nFeature_RNA):",
                                        min = 0, max = 5000, value = 200, step = 50),
                            sliderInput("max_features", "Max. Gènes par Cellule (nFeature_RNA):",
                                        min = 1000, max = 10000, value = 2500, step = 100),
                            sliderInput("max_counts", "Max. UMI par Cellule (nCount_RNA):",
                                        min = 1000, max = 100000, value = 10000, step = 1000),
                            sliderInput("max_mito", "Max. % Mitochondries (percent.mt):",
                                        min = 0, max = 100, value = 5, step = 0.5),
                            actionButton("apply_manual_qc", "Appliquer le filtrage manuel"),
                            helpText("Note: Les seuils optimaux varient selon l'expérience et la qualité des données.")
                   ),
                   tabPanel("Filtrage Automatisé (ddqcR)",
                            value = "ddqcr_qc_tab",
                            helpText("L'approche ddqcR utilise des méthodes statistiques pour identifier les cellules de mauvaise qualité."),
                            actionButton("run_ddqcr", "Lancer ddqcR pour le filtrage automatique"),
                            conditionalPanel(
                              condition = "input.run_ddqcr > 0",
                              helpText("Le processus ddqcR peut prendre quelques instants.")
                            )
                   )
                 ), # Fin tabsetPanel QC
                 hr(),
                 
                 h4("3. Détection des Doublets (scDblFinder)"),
                 h5("Paramètres scDblFinder:"),
                 sliderInput("scdbl_dbr", "dbr (Taux de doublets):",
                             min = 0.01, max = 0.2, value = 0.06, step = 0.005),
                 helpText("dbr: Taux de doublets attendu (ex: 0.06 pour 6%)."),
                 actionButton("run_doublet_detection", "Lancer la détection des doublets"),
                 helpText("Nécessite normalisation, scaling, PCA, UMAP et clustering."),
                 hr(),
                 
                 h4("4. Appliquer le Filtrage des Doublets"),
                 checkboxInput("remove_doublets_checkbox", "Retirer les cellules doublets", value = TRUE),
                 actionButton("apply_doublet_removal", "Appliquer le filtrage des doublets"),
                 helpText("Ceci créera l'objet final pour les analyses suivantes."),
                 hr()
               ), # Fin sidebarPanel Preprocessing
               
               mainPanel(
                 h3("Statut des données Seurat"),
                 textOutput("seurat_status"),
                 hr(),
                 h4("Résumé des données initiales"),
                 verbatimTextOutput("seurat_initial_summary"),
                 hr(),
                 h4("Visualisation des métriques QC (avant filtrage)"),
                 plotlyOutput("qc_metrics_plot_initial"),
                 hr(),
                 h4("Visualisation des métriques QC (après filtrage)"),
                 plotlyOutput("qc_metrics_plot_filtered"),
                 hr(),
                 conditionalPanel(
                   condition = "input.apply_manual_qc > 0 || input.run_ddqcr > 0",
                   h4("Résumé des données après filtrage"),
                   verbatimTextOutput("seurat_filtered_summary")
                 ),
                 conditionalPanel(
                   condition = "input.run_ddqcr > 0",
                   h4("Résultats ddqcR"),
                   verbatimTextOutput("ddqcr_results_summary")
                 ),
                 hr(),
                 
                 conditionalPanel(
                   condition = "input.run_doublet_detection > 0",
                   h4("Résultats de la Détection des Doublets"),
                   verbatimTextOutput("doublet_summary"),
                   plotlyOutput("doublet_plot"),
                   helpText("Les cellules en rouge sont identifiées comme des doublets, les bleues comme des singlets.")
                 ),
                 hr(),
                 
                 conditionalPanel(
                   condition = "input.apply_doublet_removal > 0",
                   h4("Résumé de l'objet Seurat final (après retrait des doublets)"),
                   verbatimTextOutput("seurat_final_summary")
                 )
               ) # Fin mainPanel Preprocessing
             ) # Fin sidebarLayout Preprocessing
    ), # Fin tabPanel "Preprocessing"
    
    # ---------------------------------------------------------------------
    # Onglet "Visualisation"
    # ---------------------------------------------------------------------
    tabPanel("Visualisation",
             icon = icon("chart-bar"), # Icône de graphique
             sidebarLayout(
               sidebarPanel(
                 h4("Visualisation UMAP"),
                 sliderInput("umap_resolution", "Résolution de Clustering (pour UMAP):",
                             min = 0.1, max = 2, value = 0.5, step = 0.1),
                 actionButton("update_umap_clustering", "Mettre à jour UMAP & Clustering"),
                 hr(),
                 
                 h4("Calcul des Module Scores"),
                 textAreaInput("module_genes_input", "Entrez les gènes pour chaque module (une ligne par module, gènes séparés par des virgules, premier mot = nom du module):",
                               rows = 5,
                               placeholder = "Exemple:\nModule1,GENE1,GENE2,GENE3\nModule2,GENE4,GENE5"),
                 numericInput("ctrl_module_score", "Paramètre 'ctrl' pour AddModuleScore:", value = 100, min = 1),
                 actionButton("calculate_module_scores", "Calculer les Module Scores"),
                 hr(),
                 
                 h4("Visualisation de l'Expression Génique & Module Scores"),
                 selectizeInput("gene_select", "Sélectionner un ou plusieurs gènes/module scores:",
                                choices = NULL, # Sera rempli dynamiquement
                                multiple = TRUE,
                                options = list(placeholder = 'Commencez à taper un nom de gène ou de module...')),
                 selectInput("plot_type", "Type de Plot d'Expression:",
                             choices = c("UMAP (FeaturePlot)" = "feature",
                                         "Violon (VlnPlot)" = "violin",
                                         "Dot Plot" = "dot",
                                         "Plot de Pathway" = "pathway_plot", # Nouvelle option
                                         "Plot de Proportion Clusters" = "baranal_plot"), # Nouvelle option
                             selected = "feature"),
                 conditionalPanel(
                   condition = "input.plot_type == 'violin' || input.plot_type == 'dot'",
                   selectInput("group_by_var", "Grouper par (pour VlnPlot/DotPlot):",
                               choices = NULL) # Sera rempli dynamiquement
                 ),
                 conditionalPanel(
                   condition = "input.plot_type == 'pathway_plot'",
                   selectInput("pathway_select", "Sélectionner un Pathway:",
                               choices = c("TCR_signaling", "B_cell_pathway", "Fas_pathway", "Il1R_signaling"),
                               selected = "TCR_signaling")
                 ),
                 actionButton("plot_gene_expression", "Afficher le Plot"),
                 hr()
               ), # Fin sidebarPanel Visualisation
               
               mainPanel(
                 h3("Visualisations"),
                 conditionalPanel(
                   condition = "input.update_umap_clustering > 0",
                   h4("UMAP par Cluster (Résolution ajustée)"),
                   plotlyOutput("clustered_umap_plot")
                 ),
                 conditionalPanel(
                   condition = "input.plot_gene_expression > 0",
                   h4("Expression Génique, Module Score, Pathway ou Proportion de Clusters"),
                   uiOutput("gene_expression_or_pathway_output") # Utiliser uiOutput pour gérer plotOutput ou plotlyOutput
                 )
               ) # Fin mainPanel Visualisation
             ) # Fin sidebarLayout Visualisation
    ), # Fin tabPanel "Visualisation"
    
    # ---------------------------------------------------------------------
    # Onglet "Analyse d'Expression Différentielle"
    # ---------------------------------------------------------------------
    tabPanel("Analyse d'Expression Différentielle",
             icon = icon("flask"), # Icône de flasque
             sidebarLayout(
               sidebarPanel(
                 h4("Paramètres de l'analyse DE"),
                 selectInput("de_method", "Méthode de test:",
                             choices = c("Wilcoxon" = "wilcox",
                                         "t-test" = "t",
                                         "Bimod" = "bimod",
                                         "MAST" = "MAST"),
                             selected = "wilcox"),
                 numericInput("de_logfc_threshold", "Seuil de log2FC (min):", value = 0.25, min = 0),
                 numericInput("de_min_pct", "Pourcentage min. de cellules (min.pct):", value = 0.1, min = 0, max = 1, step = 0.05),
                 checkboxInput("de_only_pos", "Marqueurs positifs seulement (pour FindAllMarkers):", value = FALSE),
                 hr(),
                 
                 h4("Comparer tous les clusters"),
                 actionButton("run_all_markers", "Lancer l'analyse DE (Tous les clusters)"),
                 helpText("Compare chaque cluster vs. tous les autres."),
                 hr(),
                 
                 h4("Comparer un cluster vs un autre"),
                 selectInput("de_cluster1", "Cluster 1:", choices = NULL),
                 selectInput("de_cluster2", "Cluster 2:", choices = NULL),
                 actionButton("run_cluster_vs_cluster_markers", "Lancer l'analyse DE (Cluster vs Cluster)"),
                 helpText("Compare le Cluster 1 au Cluster 2."),
                 hr(),
                 
                 h4("Affichage des résultats"),
                 sliderInput("de_pvalue_cutoff", "Seuil de p-value ajustée (max):", min = 0, max = 0.1, value = 0.05, step = 0.001),
                 sliderInput("de_display_logfc_cutoff", "Seuil de log2FC pour l'affichage (min):", min = 0, max = 2, value = 0.25, step = 0.05),
                 selectInput("de_sort_by", "Trier par:",
                             choices = c("p_val_adj", "avg_log2FC", "pct_diff"),
                             selected = "p_val_adj"),
                 checkboxInput("de_sort_desc", "Trier par ordre décroissant:", value = FALSE),
                 downloadButton("download_de_results", "Télécharger les résultats DE"),
                 hr(),
                 
                 h4("Visualisation des Top Marqueurs"),
                 sliderInput("de_num_top_markers", "Nombre de top marqueurs à visualiser:",
                             min = 1, max = 50, value = 10, step = 1),
                 actionButton("plot_top_markers", "Afficher les Top Marqueurs")
               ), # Fin sidebarPanel DE
               
               mainPanel(
                 h3("Résultats de l'Analyse d'Expression Différentielle"),
                 DT::dataTableOutput("de_results_table"),
                 hr(),
                 conditionalPanel(
                   condition = "input.plot_top_markers > 0",
                   h4("Heatmap des Top Marqueurs"),
                   plotOutput("de_top_markers_heatmap", height = "600px"), # Utilise plotOutput pour DoHeatmap
                   helpText("Les marqueurs sont colorés par expression relative, les clusters par leur identité.")
                 )
               ) # Fin mainPanel DE
             ) # Fin sidebarLayout DE
    ), # Fin tabPanel "Analyse d'Expression Différentielle"
    
    # ---------------------------------------------------------------------
    # Onglet "Modification" (Simplifié)
    # ---------------------------------------------------------------------
    tabPanel("Modification",
             icon = icon("scissors"), # Icône de ciseaux ou de modification
             sidebarLayout(
               sidebarPanel(
                 h4("1. Créer un Subset de Clusters"),
                 helpText("Sélectionnez un ou plusieurs clusters à **exclure** de l'analyse."),
                 selectInput("subset_clusters_to_exclude", "Clusters à exclure:", # Mis à jour ici
                             choices = NULL, # Sera rempli dynamiquement
                             multiple = TRUE),
                 textInput("subset_name", "Nom du nouvel objet subset:", value = "subset_data"),
                 actionButton("create_subset_obj", "Créer l'Objet Subset"),
                 hr(),
                 
                 h4("2. Exporter le Subset"),
                 helpText("Exporter l'objet subset actuellement créé sous forme de fichier .rds."),
                 downloadButton("download_subset_rds", "Télécharger le Subset en .rds"),
                 hr()
               ), # Fin sidebarPanel Modification
               
               mainPanel(
                 h3("Statut et Résumé du Subset"),
                 verbatimTextOutput("subset_summary")
               ) # Fin mainPanel Modification
             ) # Fin sidebarLayout Modification
    ) # Fin tabPanel "Modification"
  ) # Fin navbarPage
) # Fin fluidPage

# -------------------------------------------------------------------------
# Définition de la logique du serveur (Server)
# -------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # Variables réactives pour stocker les objets Seurat à différentes étapes
  initial_seurat_obj <- reactiveVal(NULL)      # L'objet Seurat chargé AVANT tout filtrage
  filtered_seurat_obj <- reactiveVal(NULL)     # L'objet Seurat APRES application des filtres QC
  classified_seurat_obj <- reactiveVal(NULL)   # Objet Seurat avec les classifications de doublets, AVANT suppression
  final_seurat_obj <- reactiveVal(NULL)        # L'objet Seurat APRES application (ou non) du filtrage des doublets
  
  # Reactive pour l'objet Seurat utilisé dans la visualisation (avec clustering potentiellement mis à jour)
  # Ce sera l'objet `final_seurat_obj()` par défaut, mais son clustering peut être modifié.
  visual_seurat_obj <- reactiveVal(NULL)
  
  # Reactive pour stocker les résultats de l'analyse DE
  de_results <- reactiveVal(NULL)
  
  # Reactive pour stocker l'objet subset
  subset_seurat_obj <- reactiveVal(NULL)
  
  
  # --- Statut de l'application ---
  output$seurat_status <- renderText({
    if (is.null(initial_seurat_obj())) {
      "Veuillez charger les données (Cell Ranger ou .rds)."
    } else if (is.null(filtered_seurat_obj())) {
      "Données chargées, en attente de filtrage QC."
    } else if (is.null(classified_seurat_obj())) {
      paste("Données QC chargées. Contient", ncol(filtered_seurat_obj()), "cellules et", nrow(filtered_seurat_obj()), "gènes. En attente de détection des doublets.")
    } else if (is.null(final_seurat_obj())) {
      paste("Doublets détectés. Contient", ncol(classified_seurat_obj()), "cellules. Appliquez le filtrage des doublets pour obtenir l'objet final.")
    } else {
      paste("Analyse terminée. Contient", ncol(final_seurat_obj()), "cellules et", nrow(final_seurat_obj()), "gènes après QC et filtrage des doublets.")
    }
  })
  
  # --- 1. Chargement des données Cell Ranger ---
  observeEvent(input$load_cellranger_data, {
    req(input$mtx_file, input$features_file, input$barcodes_file)
    
    withProgress(message = "Chargement des données Seurat (Cell Ranger)...", value = 0, {
      incProgress(0.1, detail = "Lecture des fichiers...")
      
      tryCatch({
        counts <- ReadMtx(
          mtx = input$mtx_file$datapath,
          features = input$features_file$datapath,
          cells = input$barcodes_file$datapath
        )
        incProgress(0.5, detail = "Création de l'objet Seurat et calcul des métriques...")
        
        seurat_obj <- CreateSeuratObject(counts = counts, project = "scRNAseq_QC")
        
        # Calcul du pourcentage de gènes mitochondriaux
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
        # Calcul du pourcentage de gènes ribosomaux
        seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
        
        # Ajouter une colonne 'orig.ident' si elle n'existe pas, nécessaire pour baranal
        if (!"orig.ident" %in% colnames(seurat_obj@meta.data)) {
          seurat_obj$orig.ident <- "sample1"
          showNotification("Colonne 'orig.ident' non trouvée. Création d'une colonne 'sample1' par défaut pour la démo de baranal.", type = "warning")
        }
        
        
        initial_seurat_obj(seurat_obj)
        # Initialiser l'objet filtré et final avec l'objet initial après chargement pour que le workflow puisse commencer
        filtered_seurat_obj(seurat_obj)
        classified_seurat_obj(NULL)
        final_seurat_obj(NULL)
        visual_seurat_obj(NULL) # Réinitialiser l'objet de visualisation
        subset_seurat_obj(NULL) # Réinitialiser l'objet subset
        
        incProgress(1, detail = "Chargement Cell Ranger terminé.")
        showNotification(paste0("Données Cell Ranger chargées. Initialement: ", ncol(seurat_obj), " cellules."), type = "message")
        
      }, error = function(e) {
        showNotification(paste("Erreur lors du chargement des données Cell Ranger:", e$message), type = "error")
        initial_seurat_obj(NULL)
        filtered_seurat_obj(NULL)
        classified_seurat_obj(NULL)
        final_seurat_obj(NULL)
        visual_seurat_obj(NULL)
        subset_seurat_obj(NULL)
      })
    })
  })
  
  # --- 1. Chargement d'un fichier Seurat (.rds) ---
  observeEvent(input$load_rds_data, {
    req(input$rds_file)
    
    withProgress(message = "Chargement de l'objet Seurat (.rds)...", value = 0, {
      incProgress(0.1, detail = "Lecture du fichier RDS...")
      
      tryCatch({
        seurat_obj <- readRDS(input$rds_file$datapath)
        
        if (!inherits(seurat_obj, "Seurat")) {
          stop("Le fichier chargé n'est pas un objet Seurat valide.")
        }
        
        # Vérification si les métriques QC sont présentes, sinon les calculer
        if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
          showNotification("Calcul de percent.mt (non trouvé dans le RDS).", type = "message")
          seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
        }
        if (!"percent.ribo" %in% colnames(seurat_obj@meta.data)) {
          showNotification("Calcul de percent.ribo (non trouvé dans le RDS).", type = "message")
          seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
        }
        
        # Vérifier si 'orig.ident' existe, sinon le créer (nécessaire pour baranal)
        if (!"orig.ident" %in% colnames(seurat_obj@meta.data)) {
          showNotification("Colonne 'orig.ident' non trouvée. Création d'une colonne 'sample1' par défaut pour la démo de baranal.", type = "warning")
          seurat_obj$orig.ident <- "sample1" # Par défaut à un seul échantillon
        }
        
        initial_seurat_obj(seurat_obj)
        # L'objet RDS est considéré comme "initial" et peut être filtré ou non selon son état
        filtered_seurat_obj(seurat_obj)
        classified_seurat_obj(seurat_obj) # Si le RDS est déjà traité, on le considère classifié
        final_seurat_obj(seurat_obj)       # Et final pour la visualisation
        visual_seurat_obj(seurat_obj)      # L'objet de visualisation est maintenant le RDS chargé
        subset_seurat_obj(NULL)            # Réinitialiser l'objet subset
        
        incProgress(1, detail = "Chargement RDS terminé.")
        showNotification(paste0("Objet Seurat (.rds) chargé. Contient: ", ncol(seurat_obj), " cellules."), type = "message") # CORRECTION: type = "success" -> "message"
        
      }, error = function(e) {
        showNotification(paste("Erreur lors du chargement du fichier RDS:", e$message), type = "error")
        initial_seurat_obj(NULL)
        filtered_seurat_obj(NULL)
        classified_seurat_obj(NULL)
        final_seurat_obj(NULL)
        visual_seurat_obj(NULL)
        subset_seurat_obj(NULL)
      })
    })
  })
  
  
  # --- Résumé des données initiales ---
  output$seurat_initial_summary <- renderPrint({
    obj <- initial_seurat_obj()
    req(obj)
    print(obj)
    cat("\n")
    cat("Résumé des métriques initiales:\n")
    cat("nFeature_RNA (Gènes par Cellule):\n")
    print(summary(obj@meta.data$nFeature_RNA))
    cat("nCount_RNA (UMI par Cellule):\n")
    print(summary(obj@meta.data$nCount_RNA))
    cat("percent.mt (% Mitochondries):\n")
    print(summary(obj@meta.data$percent.mt))
    cat("percent.ribo (% Ribosomales):\n")
    print(summary(obj@meta.data$percent.ribo))
  })
  
  # --- Plot des métriques QC initiales ---
  output$qc_metrics_plot_initial <- renderPlotly({
    obj <- initial_seurat_obj()
    req(obj)
    
    plot_nFeature <- ggplot(obj@meta.data, aes(x = nFeature_RNA)) +
      geom_histogram(binwidth = 20, fill = "lightblue", color = "black") +
      ggtitle("Distribution de nFeature_RNA (Gènes par Cellule)") +
      theme_minimal()
    
    plot_nCount <- ggplot(obj@meta.data, aes(x = nCount_RNA)) +
      geom_histogram(binwidth = 500, fill = "lightgreen", color = "black") +
      ggtitle("Distribution de nCount_RNA (UMI par Cellule)") +
      theme_minimal()
    
    plot_percent_mt <- ggplot(obj@meta.data, aes(x = percent.mt)) +
      geom_histogram(binwidth = 0.5, fill = "salmon", color = "black") +
      ggtitle("Distribution de percent.mt (% Mitochondries)") +
      theme_minimal()
    
    subplot(ggplotly(plot_nFeature), ggplotly(plot_nCount), ggplotly(plot_percent_mt),
            nrows = 1, shareY = FALSE, titleX = TRUE) %>%
      layout(title = "Métriques de Contrôle Qualité (Avant Filtrage)")
  })
  
  # --- 2. Logique de filtrage manuel ---
  observeEvent(input$apply_manual_qc, {
    obj <- initial_seurat_obj()
    req(obj)
    showNotification("Application du filtrage manuel...", type = "message")
    
    withProgress(message = "Application du filtrage manuel...", value = 0.5, {
      filtered_obj <- subset(obj,
                             subset = nFeature_RNA > input$min_features &
                               nFeature_RNA < input$max_features &
                               nCount_RNA < input$max_counts &
                               percent.mt < input$max_mito)
      filtered_seurat_obj(filtered_obj)
      classified_seurat_obj(NULL)
      final_seurat_obj(NULL)
      visual_seurat_obj(NULL) # Réinitialiser l'objet de visualisation
      subset_seurat_obj(NULL) # Réinitialiser l'objet subset
    })
    showNotification(paste0("Filtrage manuel appliqué. Cellules restantes: ", ncol(filtered_seurat_obj()), "."), type = "message")
  })
  
  # --- 2. Logique de filtrage ddqcR ---
  observeEvent(input$run_ddqcr, {
    obj <- initial_seurat_obj()
    req(obj)
    showNotification("Lancement du filtrage ddqcR...", type = "message")
    
    withProgress(message = "Exécution de ddqcR...", value = 0, {
      # SIMULATION DE ddqcR pour l'exemple. Décommenter et implémenter la vraie logique ddqcR si nécessaire.
      initial_metadata <- obj@meta.data
      low_feature_threshold <- quantile(initial_metadata$nFeature_RNA, 0.05, na.rm = TRUE)
      high_mito_threshold <- quantile(initial_metadata$percent.mt, 0.95, na.rm = TRUE)
      
      initial_metadata$ddqc_call <- ifelse(
        initial_metadata$nFeature_RNA < low_feature_threshold |
          initial_metadata$percent.mt > high_mito_threshold,
        "Bad_Quality_Simulated",
        "Good_Quality_Simulated"
      )
      obj@meta.data <- initial_metadata
      
      incProgress(0.7, detail = "Filtrage basé sur ddqcR simulé...")
      filtered_obj <- subset(obj, subset = ddqc_call == "Good_Quality_Simulated")
      filtered_seurat_obj(filtered_obj)
      classified_seurat_obj(NULL)
      final_seurat_obj(NULL)
      visual_seurat_obj(NULL) # Réinitialiser l'objet de visualisation
      subset_seurat_obj(NULL) # Réinitialiser l'objet subset
      
      showNotification("Filtrage ddqcR (simulé) appliqué avec succès.", type = "message")
      
      output$ddqcr_results_summary <- renderPrint({
        cat("Répartition des cellules par appel ddqcR (simulé):\n")
        print(table(obj@meta.data$ddqc_call))
      })
    })
  })
  
  # --- Résumé des données après filtrage QC ---
  output$seurat_filtered_summary <- renderPrint({
    obj <- filtered_seurat_obj()
    req(obj)
    print(obj)
    cat("\n")
    cat("Résumé des métriques après filtrage QC:\n")
    cat("Cellules restantes: ", ncol(obj), "\n")
    cat("Gènes restants: ", nrow(obj), "\n")
    cat("nFeature_RNA (Gènes par Cellule):\n")
    print(summary(obj@meta.data$nFeature_RNA))
    cat("nCount_RNA (UMI par Cellule):\n")
    print(summary(obj@meta.data$nCount_RNA))
    cat("percent.mt (% Mitochondries):\n")
    print(summary(obj@meta.data$percent.mt))
    cat("percent.ribo (% Ribosomales):\n")
    print(summary(obj@meta.data$percent.ribo))
  })
  
  # --- Plot des métriques QC après filtrage ---
  output$qc_metrics_plot_filtered <- renderPlotly({
    obj <- filtered_seurat_obj()
    req(obj)
    
    plot_nFeature <- ggplot(obj@meta.data, aes(x = nFeature_RNA)) +
      geom_histogram(binwidth = 20, fill = "lightblue", color = "black") +
      ggtitle("Distribution de nFeature_RNA (Après Filtrage)") +
      theme_minimal()
    
    plot_nCount <- ggplot(obj@meta.data, aes(x = nCount_RNA)) +
      geom_histogram(binwidth = 500, fill = "lightgreen", color = "black") +
      ggtitle("Distribution de nCount_RNA (UMI par Cellule)") +
      theme_minimal()
    
    plot_percent_mt <- ggplot(obj@meta.data, aes(x = percent.mt)) +
      geom_histogram(binwidth = 0.5, fill = "salmon", color = "black") +
      ggtitle("Distribution de percent.mt (% Mitochondries)") +
      theme_minimal()
    
    subplot(ggplotly(plot_nFeature), ggplotly(plot_nCount), ggplotly(plot_percent_mt),
            nrows = 1, shareY = FALSE, titleX = TRUE) %>%
      layout(title = "Métriques de Contrôle Qualité (Après Filtrage)")
  })
  
  # --- 3. Préparation de l'objet Seurat pour la détection des doublets ---
  processed_seurat_for_doublets <- eventReactive(input$run_doublet_detection, {
    obj <- filtered_seurat_obj() # Utilise l'objet après QC
    req(obj)
    
    if (ncol(obj) < 50) {
      showNotification(paste0("Moins de 50 cellules restantes (", ncol(obj), ") après QC. La détection des doublets peut échouer ou être non fiable."), type = "warning")
      validate("Trop peu de cellules pour exécuter la détection des doublets.")
    }
    
    # S'assurer que l'assay par défaut est "RNA"
    DefaultAssay(obj) <- "RNA"
    
    withProgress(message = "Préparation de l'objet Seurat pour détection doublets...", value = 0, {
      incProgress(0.1, detail = "Normalisation des données...")
      obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
      
      incProgress(0.3, detail = "Identification des gènes variables...")
      obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
      if (length(VariableFeatures(obj)) == 0) {
        showNotification("Aucun gène variable trouvé. La PCA/UMAP pourrait échouer.", type = "warning")
        validate("Aucun gène variable détecté. Impossible de procéder.")
      }
      
      incProgress(0.5, detail = "Mise à l'échelle des données...")
      obj <- ScaleData(obj)
      
      incProgress(0.7, detail = "Exécution de la PCA...")
      obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 30)
      
      # Vérification post-PCA
      if (is.null(obj@reductions$pca) || nrow(obj@reductions$pca@cell.embeddings) == 0) {
        showNotification("Erreur: PCA n'a pas pu être exécutée ou n'a pas produit de résultats valides. Vérifiez les données.", type = "error")
        validate("Problème avec la PCA. Impossible de procéder à la détection des doublets.")
      }
      
      num_pcs_for_embedding <- min(30, ncol(obj@reductions$pca@cell.embeddings))
      if (num_pcs_for_embedding == 0) {
        showNotification("Aucune PC valide disponible après PCA. Impossible de procéder au clustering/UMAP pour la détection des doublets.", type = "error")
        validate("Pas de PCs pour la détection des doublets.")
      }
      if (num_pcs_for_embedding < 5 && ncol(obj) > 100) { # Avertir si le nombre de PCs est très faible
        showNotification(paste0("Seulement ", num_pcs_for_embedding, " PCs calculées. Les résultats de détection des doublets pourraient être affectés."), type = "warning")
      }
      
      
      incProgress(0.8, detail = "Clustering des cellules...")
      obj <- FindNeighbors(obj, dims = 1:num_pcs_for_embedding)
      obj <- FindClusters(obj, resolution = 0.5) # Résolution par défaut pour la détection initiale
      
      incProgress(0.9, detail = "Exécution de l'UMAP...")
      obj <- RunUMAP(obj, dims = 1:num_pcs_for_embedding)
      
      incProgress(1, detail = "Préparation terminée.")
    })
    return(obj)
  })
  
  
  # --- 3. Logique de détection des doublets avec scDblFinder (avec conversion SCE) ---
  observeEvent(input$run_doublet_detection, {
    obj <- processed_seurat_for_doublets() # Déclenche la préparation et récupère l'objet
    req(obj)
    
    tryCatch({
      withProgress(message = "Détection des doublets avec scDblFinder...", value = 0, {
        
        incProgress(0.1, detail = "Conversion en SingleCellExperiment...")
        sce <- as.SingleCellExperiment(obj)
        
        # S'assurer que les réductions sont présentes dans l'objet SCE si scDblFinder en a besoin
        # scDblFinder peut utiliser les embeddings existants
        if ("umap" %in% names(obj@reductions) && !("UMAP" %in% reducedDimNames(sce))) {
          reducedDim(sce, "UMAP") <- obj@reductions$umap@cell.embeddings
        }
        if ("pca" %in% names(obj@reductions) && !("PCA" %in% reducedDimNames(sce))) {
          reducedDim(sce, "PCA") <- obj@reductions$pca@cell.embeddings
        }
        
        incProgress(0.5, detail = "Exécution de scDblFinder...")
        sce <- scDblFinder(sce, dbr = input$scdbl_dbr)
        
        # Récupérer les résultats de scDblFinder et les re-stocker dans l'objet Seurat original
        obj@meta.data$scDblFinder.score <- sce$scDblFinder.score
        obj@meta.data$scDblFinder.class <- as.character(sce$scDblFinder.class)
        obj@meta.data$is_doublet <- obj@meta.data$scDblFinder.class == "doublet"
        
        classified_seurat_obj(obj) # Met à jour l'objet avec la classification des doublets
        final_seurat_obj(NULL) # Réinitialise l'objet final, il doit être appliqué manuellement
        visual_seurat_obj(NULL) # Réinitialiser l'objet de visualisation
        subset_seurat_obj(NULL) # Réinitialiser l'objet subset
        
        incProgress(1, detail = "Détection des doublets terminée.")
        showNotification("Détection des doublets avec scDblFinder terminée. Veuillez appliquer le filtrage.", type = "message")
        
        output$doublet_summary <- renderPrint({
          cat("Résumé de la détection des doublets (scDblFinder):\n")
          if (!is.null(obj@meta.data$scDblFinder.class)) {
            print(table(obj@meta.data$scDblFinder.class))
            cat(paste0("\nNombre de cellules total (avant filtrage): ", ncol(obj), "\n"))
            cat(paste0("Nombre de singlets identifiés: ", sum(obj@meta.data$is_doublet == FALSE), "\n"))
            cat(paste0("Nombre de doublets identifiés: ", sum(obj@meta.data$is_doublet == TRUE), "\n"))
          } else {
            cat("Aucun résultat de doublet à afficher.\n")
          }
          cat("\nDistribution des scores scDblFinder:\n")
          print(summary(obj@meta.data$scDblFinder.score))
        })
      })
    }, error = function(e) {
      showNotification(paste("Erreur critique lors de la détection des doublets avec scDblFinder (via SCE):", e$message), type = "error")
      classified_seurat_obj(NULL)
      final_seurat_obj(NULL)
      visual_seurat_obj(NULL)
      subset_seurat_obj(NULL)
    })
  })
  
  # --- Plot des doublets sur UMAP ---
  output$doublet_plot <- renderPlotly({
    req(input$run_doublet_detection > 0)
    obj_to_plot <- classified_seurat_obj()
    req(obj_to_plot)
    
    if (is.null(obj_to_plot@reductions$umap)) {
      showNotification("UMAP non disponible pour la visualisation des doublets. Vérifiez les étapes de préparation.", type = "warning")
      return(NULL)
    }
    if (!"is_doublet" %in% colnames(obj_to_plot@meta.data)) {
      showNotification("La classification des doublets n'est pas disponible pour le plot. Exécutez la détection.", type = "warning")
      return(NULL)
    }
    obj_to_plot@meta.data$is_doublet <- factor(obj_to_plot@meta.data$is_doublet, levels = c("FALSE", "TRUE"))
    
    plot <- DimPlot(obj_to_plot, reduction = "umap", group.by = "is_doublet", cols = c("FALSE" = "blue", "TRUE" = "red")) +
      ggtitle("UMAP: Cellules classées Doublets (Rouge) / Singlets (Bleu)")
    
    ggplotly(plot)
  })
  
  # --- 4. Logique d'application du filtrage des doublets ---
  observeEvent(input$apply_doublet_removal, {
    obj <- classified_seurat_obj()
    req(obj)
    
    withProgress(message = "Application du filtrage des doublets...", value = 0.5, {
      if (input$remove_doublets_checkbox) {
        final_obj <- subset(obj, subset = is_doublet == FALSE)
        showNotification(paste0("Doublets retirés. Cellules restantes: ", ncol(final_obj), "."), type = "message")
      } else {
        final_obj <- obj
        showNotification("Les doublets n'ont pas été retirés (option décochée).", type = "message")
      }
      final_seurat_obj(final_obj)
      visual_seurat_obj(final_obj) # L'objet de visualisation commence avec l'objet final
      subset_seurat_obj(NULL) # Réinitialiser l'objet subset
    })
  })
  
  # --- Résumé de l'objet Seurat final (après filtrage des doublets) ---
  output$seurat_final_summary <- renderPrint({
    obj <- final_seurat_obj()
    req(obj)
    print(obj)
    cat("\n")
    cat("Résumé des métriques de l'objet final (après QC et filtrage des doublets):\n")
    cat("Cellules restantes: ", ncol(obj), "\n")
    cat("Gènes restants: ", nrow(obj), "\n")
    cat("nFeature_RNA (Gènes par Cellule):\n")
    print(summary(obj@meta.data$nFeature_RNA))
    cat("nCount_RNA (UMI par Cellule):\n")
    print(summary(obj@meta.data$nCount_RNA))
    cat("percent.mt (% Mitochondries):\n")
    print(summary(obj@meta.data$percent.mt))
    cat("percent.ribo (% Ribosomales):\n")
    print(summary(obj@meta.data$percent.ribo))
  })
  
  # ---------------------------------------------------------------------
  # Logique de l'onglet "Visualisation"
  # ---------------------------------------------------------------------
  
  # Observer l'objet visual_seurat_obj et mettre à jour les choix de gènes/module scores et de variables de groupement
  observeEvent(visual_seurat_obj(), {
    obj <- visual_seurat_obj()
    req(obj)
    
    # Mettre à jour les gènes disponibles dans le sélecteur
    # On inclut tous les gènes et les colonnes numériques qui ne sont PAS les métriques QC standard
    all_features_and_scores <- rownames(obj)
    
    # Ajouter toutes les colonnes numériques de metadata (qui pourraient être des module scores)
    numeric_meta_cols <- names(obj@meta.data)[sapply(obj@meta.data, is.numeric)]
    # Exclure les métriques QC standard (nFeature_RNA, nCount_RNA, percent.mt, percent.ribo)
    qc_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo", "scDblFinder.score")
    module_score_cols <- setdiff(numeric_meta_cols, qc_metrics)
    
    all_features_and_scores <- c(all_features_and_scores, module_score_cols)
    
    updateSelectizeInput(session, "gene_select", choices = all_features_and_scores, server = TRUE)
    
    # Mettre à jour les variables de groupement disponibles pour VlnPlot/DotPlot
    grouping_candidates <- names(obj@meta.data)[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
    
    # Ajouter 'seurat_clusters' s'il existe et n'est pas déjà inclus, et le mettre en première position
    if ("seurat_clusters" %in% colnames(obj@meta.data)) {
      grouping_candidates <- c("seurat_clusters", grouping_candidates[!grouping_candidates %in% "seurat_clusters"])
    }
    
    if(is.null(obj@meta.data$seurat_clusters) && is.null(grouping_candidates)){
      updateSelectInput(session, "group_by_var", choices = c("Veuillez clusteriser d'abord" = ""), selected = "")
    } else {
      updateSelectInput(session, "group_by_var", choices = unique(grouping_candidates), selected = "seurat_clusters")
    }
  })
  
  # --- Génération de la palette de couleurs ---
  # Cette réactive est déclenchée chaque fois que l'objet visuel (et donc les clusters) change
  cluster_colors <- reactive({
    obj <- visual_seurat_obj()
    req(obj)
    if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
      return(NULL) # Pas de clusters pour générer une palette
    }
    n_clusters <- length(levels(obj@meta.data$seurat_clusters))
    # Utiliser la palette par défaut de ggplot2/Seurat
    set.seed(42) # Pour la reproductibilité des couleurs
    colors <- scales::hue_pal()(n_clusters)
    
    # Stocker les couleurs dans l'objet Seurat pour la fonction baranal si elle en a besoin
    # Note: On doit mettre à jour l'objet réactif lui-même pour propager la modification de @misc
    # C'est un pattern un peu avancé en Shiny, mais utile ici.
    obj_copy <- obj # Créer une copie pour éviter de modifier l'objet réactif directement dans le bloc `reactive`
    obj_copy@misc$cluster_colors <- colors
    visual_seurat_obj(obj_copy) # Met à jour l'objet avec les couleurs stockées
    return(colors)
  })
  
  
  # --- Mise à jour du clustering UMAP ---
  observeEvent(input$update_umap_clustering, {
    obj <- visual_seurat_obj() # Utilise l'objet de visualisation actuel
    req(obj)
    
    if (is.null(obj@reductions$pca)) {
      showNotification("Veuillez d'abord exécuter les étapes de pré-traitement (Normalisation, PCA, etc.) dans l'onglet 'Preprocessing'.", type = "warning")
      return(NULL)
    }
    
    # Vérification robuste du nombre de PCs
    num_pcs_for_embedding <- min(30, ncol(obj@reductions$pca@cell.embeddings))
    if (num_pcs_for_embedding == 0) {
      showNotification("Aucune PC valide disponible pour le clustering/UMAP. Vérifiez les données ou les étapes précédentes.", type = "error")
      return(NULL)
    }
    if (num_pcs_for_embedding < 5 && ncol(obj) > 100) { # Avertir si le nombre de PCs est très faible
      showNotification(paste0("Seulement ", num_pcs_for_embedding, " PCs calculées. Les résultats de clustering/UMAP pourraient être affectés."), type = "warning")
    }
    
    # S'assurer que l'assay par défaut est "RNA"
    DefaultAssay(obj) <- "RNA"
    
    withProgress(message = "Mise à jour du clustering UMAP...", value = 0, {
      incProgress(0.5, detail = "Re-clustering...")
      # Re-cluster avec la nouvelle résolution
      obj <- FindNeighbors(obj, dims = 1:num_pcs_for_embedding)
      obj <- FindClusters(obj, resolution = input$umap_resolution)
      
      incProgress(1, detail = "Mise à jour terminée.")
      visual_seurat_obj(obj) # Met à jour l'objet réactif de visualisation
      showNotification(paste0("Clustering mis à jour avec résolution: ", input$umap_resolution), type = "message")
    })
  })
  
  # --- Calcul des Module Scores ---
  observeEvent(input$calculate_module_scores, {
    obj <- visual_seurat_obj() # On travaille sur l'objet de visualisation actuel
    req(obj)
    
    # S'assurer que l'assay par défaut est "RNA"
    DefaultAssay(obj) <- "RNA"
    
    gene_list_input <- input$module_genes_input
    if (is.null(gene_list_input) || gene_list_input == "") {
      showNotification("Veuillez entrer au moins une liste de gènes pour le calcul des module scores.", type = "warning")
      return(NULL)
    }
    
    lines <- strsplit(gene_list_input, "\n")[[1]]
    
    # Prepare list for AddModuleScore
    modules_to_score <- list()
    module_names <- c()
    
    for (line in lines) {
      parts <- strsplit(line, ",")[[1]] %>% trimws()
      if (length(parts) < 2) {
        showNotification(paste("Ligne ignorée (format incorrect):", line, ". Nécessite au moins un nom de module et un gène."), type = "warning")
        next
      }
      
      module_name <- parts[1]
      genes <- parts[-1]
      
      # Vérifier que les gènes existent dans l'objet Seurat
      genes_present <- intersect(genes, rownames(obj))
      if (length(genes_present) == 0) {
        showNotification(paste0("Aucun gène du module '", module_name, "' n'a été trouvé dans les données. Ce module sera ignoré."), type = "warning")
        next
      } else if (length(genes_present) < length(genes)) {
        showNotification(paste0("Attention: Certains gènes du module '", module_name, "' n'ont pas été trouvés et seront ignorés. (Gènes absents: ", paste(setdiff(genes, genes_present), collapse = ", "), ")"), type = "warning")
      }
      
      modules_to_score[[module_name]] <- genes_present
      module_names <- c(module_names, module_name)
    }
    
    if (length(modules_to_score) == 0) {
      showNotification("Aucun module valide n'a pu être extrait de l'entrée. Veuillez vérifier le format.", type = "error")
      return(NULL)
    }
    
    withProgress(message = "Calcul des Module Scores...", value = 0, {
      # AddModuleScore ajoute les scores dans la métadata sous le nom {name}1, {name}2, etc.
      # Donc, si le nom est "Module1", la colonne sera "Module11".
      # Je vais s'assurer que le nom de la colonne soit simplement "Module1" + le suffixe "_score"
      
      tryCatch({
        # AddModuleScore modifies the object directly
        obj <- AddModuleScore(
          object = obj,
          features = modules_to_score,
          ctrl = input$ctrl_module_score,
          name = module_names
        )
        
        # Renommer les colonnes si nécessaire pour qu'elles correspondent exactement aux noms des modules
        # AddModuleScore nomme les colonnes `name` suivi de 1 (ex: Module11, Module21)
        # On veut qu'elles soient `name_score` pour être plus claires.
        for (i in seq_along(module_names)) {
          old_col_name <- paste0(module_names[i], "1") # Nom par défaut généré par Seurat
          new_col_name <- paste0(module_names[i], "_score") # Nom désiré
          if (old_col_name %in% colnames(obj@meta.data)) {
            obj@meta.data[[new_col_name]] <- obj@meta.data[[old_col_name]]
            obj@meta.data[[old_col_name]] <- NULL # Supprimer l'ancienne colonne si renommée
          }
        }
        
        visual_seurat_obj(obj) # Mettre à jour l'objet réactif de visualisation
        showNotification("Calcul des Module Scores terminé et ajoutés à l'objet Seurat.", type = "message")
        
        # Mettre à jour immédiatement la liste des gènes/features sélectionnables
        # Appeler l'observeEvent de final_seurat_obj() pour re-remplir les selectizeInput
        # On ne peut pas déclencher un observeEvent manuellement, on doit faire la mise à jour ici
        current_choices <- rownames(obj)
        module_score_cols_added <- paste0(module_names, "_score")
        all_features_and_scores <- c(current_choices, module_score_cols_added) # S'assurer que les anciens choix sont aussi là
        
        updateSelectizeInput(session, "gene_select", choices = all_features_and_scores, server = TRUE)
        
      }, error = function(e) {
        showNotification(paste("Erreur lors du calcul des Module Scores:", e$message), type = "error")
      })
    })
  })
  
  
  # --- Plot UMAP avec clustering ---
  output$clustered_umap_plot <- renderPlotly({
    obj <- visual_seurat_obj() # Utilise l'objet avec le clustering potentiellement mis à jour
    req(obj)
    
    if (is.null(obj@reductions$umap) || !("seurat_clusters" %in% colnames(obj@meta.data)) || !is.factor(obj@meta.data$seurat_clusters)) {
      showNotification("UMAP ou Clusters non disponibles ou mal formatés. Assurez-vous d'avoir exécuté la préparation et le clustering.", type = "warning")
      return(NULL)
    }
    
    # Récupérer les couleurs générées
    colors_for_clusters <- cluster_colors() # Appelle la réactive pour s'assurer que les couleurs sont à jour
    
    plot <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
      ggtitle(paste0("UMAP des Cellules par Cluster (Résolution: ", input$umap_resolution, ")")) +
      NoLegend() + # Souvent préférable d'enlever la légende si on a beaucoup de clusters
      scale_color_manual(values = colors_for_clusters) # Appliquer la palette
    
    ggplotly(plot) %>% layout(hovermode = "closest", legend = list(x = 100, y = 0)) # Améliore la légende pour Plotly
  })
  
  # --- Plot de l'expression génique (FeaturePlot, VlnPlot, DotPlot, Pathway, Baranal) ---
  # Utilise uiOutput pour basculer entre plotlyOutput et plotOutput
  output$gene_expression_or_pathway_output <- renderUI({
    # Ces types de plots renvoient des objets ggplot2 standards, donc plotOutput est approprié.
    # FeaturePlot, VlnPlot, DotPlot de scCustomize retournent des objets ggplot qui peuvent être convertis en plotly.
    if (input$plot_type == "baranal_plot" || input$plot_type == "pathway_plot") {
      plotOutput("plot_output_static") # Renvoie un plot statique pour baranal et pathway (qui est ggplot2)
    } else {
      plotlyOutput("plot_output_interactive") # Renvoie un plot interactif pour feature/violin/dot
    }
  })
  
  output$plot_output_interactive <- renderPlotly({
    obj <- visual_seurat_obj()
    req(obj)
    req(input$plot_type)
    
    if (input$plot_type == "feature") {
      req(input$gene_select)
      # Filtrer les gènes/modules réellement présents dans l'objet
      features_to_plot <- intersect(input$gene_select, c(rownames(obj), colnames(obj@meta.data)))
      
      if (length(features_to_plot) == 0) {
        showNotification("Aucun des gènes ou module scores sélectionnés n'est trouvé dans les données.", type = "warning")
        return(NULL)
      }
      plots_list <- list()
      for (feature in features_to_plot) {
        current_p <- scCustomize::FeaturePlot_scCustom(seurat_object = obj, features = feature, reduction = "umap") +
          ggtitle(paste("Expression/Score de", feature, "sur UMAP")) +
          theme_minimal()
        plots_list[[feature]] <- ggplotly(current_p)
      }
      if (length(plots_list) == 1) {
        return(plots_list[[1]])
      } else {
        p <- subplot(plots_list, nrows = ceiling(sqrt(length(plots_list))), shareX = FALSE, shareY = FALSE) %>%
          layout(title = "Visualisation des Features sur UMAP (FeaturePlot_scCustom)")
        return(p)
      }
    } else if (input$plot_type == "violin") {
      req(input$gene_select, input$group_by_var)
      features_to_plot <- intersect(input$gene_select, c(rownames(obj), colnames(obj@meta.data)))
      if (length(features_to_plot) == 0) {
        showNotification("Aucun des gènes ou module scores sélectionnés n'est trouvé dans les données.", type = "warning")
        return(NULL)
      }
      if (!input$group_by_var %in% colnames(obj@meta.data)) {
        showNotification(paste0("Variable de groupement '", input$group_by_var, "' non trouvée."), type = "error")
        return(NULL)
      }
      
      p <- scCustomize::VlnPlot_scCustom(seurat_object = obj, features = features_to_plot, group.by = input$group_by_var) +
        ggtitle(paste("Distribution des Features par", input$group_by_var)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      if (input$group_by_var == "seurat_clusters" && !is.null(cluster_colors())) {
        p <- p + scale_fill_manual(values = cluster_colors())
      }
      return(ggplotly(p))
    } else if (input$plot_type == "dot") {
      req(input$gene_select, input$group_by_var)
      features_to_plot <- intersect(input$gene_select, c(rownames(obj), colnames(obj@meta.data)))
      if (length(features_to_plot) == 0) {
        showNotification("Aucun des gènes ou module scores sélectionnés n'est trouvé dans les données.", type = "warning")
        return(NULL)
      }
      if (!input$group_by_var %in% colnames(obj@meta.data)) {
        showNotification(paste0("Variable de groupement '", input$group_by_var, "' non trouvée."), type = "error")
        return(NULL)
      }
      
      p <- scCustomize::DotPlot_scCustom(seurat_object = obj, features = features_to_plot, group.by = input$group_by_var) +
        ggtitle(paste("Dot Plot des Features par", input$group_by_var)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      if (input$group_by_var == "seurat_clusters" && !is.null(cluster_colors())) {
        p <- p + scale_color_manual(values = cluster_colors())
      }
      return(ggplotly(p))
    }
    return(NULL) # Retourne NULL si le type de plot n'est pas interactif
  })
  
  output$plot_output_static <- renderPlot({
    obj <- visual_seurat_obj()
    req(obj)
    req(input$plot_type)
    
    if (input$plot_type == "pathway_plot") {
      req(input$pathway_select)
      tryCatch({
        p <- get_pathway(seurat_object = obj, pathway_name = input$pathway_select)
        print(p) # Doit utiliser print() pour renderPlot
      }, error = function(e) {
        showNotification(paste("Erreur lors de la génération du plot de pathway:", e$message), type = "error")
        return(NULL)
      })
    } else if (input$plot_type == "baranal_plot") {
      # Assurez-vous que la colonne 'orig.ident' existe pour baranal
      if (!"orig.ident" %in% colnames(obj@meta.data)) {
        showNotification("La colonne 'orig.ident' n'est pas trouvée dans les métadonnées. Impossible de générer le plot de proportion de clusters. Veuillez la créer ou charger un objet Seurat avec cette colonne.", type = "error")
        return(NULL)
      }
      tryCatch({
        # `cluster_colors()` met à jour l'objet `visual_seurat_obj` avec les couleurs dans `@misc`.
        # On doit récupérer la version la plus récente de l'objet pour s'assurer que `baranal` y a accès.
        current_obj_with_colors <- visual_seurat_obj()
        p <- baranal(seurat_object = current_obj_with_colors)
        print(p) # Doit utiliser print() pour renderPlot
      }, error = function(e) {
        showNotification(paste("Erreur lors de la génération du plot de proportion de clusters:", e$message), type = "error")
        return(NULL)
      })
    }
    return(NULL) # Retourne NULL si le type de plot n'est pas statique
  })
  
  # Observe le bouton pour afficher l'expression des gènes
  observeEvent(input$plot_gene_expression, {
    # Ce bouton déclenche simplement le renderPlotly/renderPlot pour les gènes/modules/pathways/baranal
    # Pas besoin de logique supplémentaire ici car les renderFunctions sont réactives aux inputs
  })
  
  # ---------------------------------------------------------------------
  # Logique de l'onglet "Analyse d'Expression Différentielle"
  # ---------------------------------------------------------------------
  
  # Mettre à jour les choix de clusters pour la comparaison 1 vs 1
  observeEvent(visual_seurat_obj(), {
    obj <- visual_seurat_obj()
    req(obj)
    if (!"seurat_clusters" %in% colnames(obj@meta.data) || !is.factor(obj@meta.data$seurat_clusters)) {
      updateSelectInput(session, "de_cluster1", choices = c("Clustering non disponible" = ""))
      updateSelectInput(session, "de_cluster2", choices = c("Clustering non disponible" = ""))
      return(NULL)
    }
    
    cluster_levels <- levels(obj$seurat_clusters)
    updateSelectInput(session, "de_cluster1", choices = cluster_levels, selected = head(cluster_levels, 1))
    # Cluster 2 doit être différent du Cluster 1 par défaut
    if (length(cluster_levels) > 1) {
      updateSelectInput(session, "de_cluster2", choices = cluster_levels, selected = tail(cluster_levels, 1))
    } else {
      updateSelectInput(session, "de_cluster2", choices = cluster_levels, selected = cluster_levels[1]) # Si un seul cluster, sélectionne le même
    }
  })
  
  # --- Lancer l'analyse DE (Tous les clusters) ---
  observeEvent(input$run_all_markers, {
    obj <- visual_seurat_obj()
    req(obj)
    
    if (!"seurat_clusters" %in% colnames(obj@meta.data) || !is.factor(obj@meta.data$seurat_clusters)) {
      showNotification("Le clustering (colonne 'seurat_clusters') n'est pas disponible ou est mal formaté. Veuillez clusteriser d'abord.", type = "error")
      return(NULL)
    }
    
    withProgress(message = "Calcul des marqueurs pour tous les clusters...", value = 0, {
      
      # Set default identity to seurat_clusters if not already set
      DefaultAssay(obj) <- "RNA"
      obj <- SetIdent(obj, value = "seurat_clusters")
      
      incProgress(0.3, detail = "Exécution de FindAllMarkers...")
      tryCatch({
        markers <- FindAllMarkers(
          object = obj,
          assay = "RNA",
          only.pos = input$de_only_pos,
          min.pct = input$de_min_pct,
          logfc.threshold = input$de_logfc_threshold,
          test.use = input$de_method
        )
        
        # Calculate Pct_diff
        markers <- markers %>%
          dplyr::mutate(pct_diff = abs(pct.1 - pct.2))
        
        de_results(markers)
        showNotification("Analyse FindAllMarkers terminée.", type = "message")
      }, error = function(e) {
        showNotification(paste("Erreur lors de FindAllMarkers:", e$message), type = "error")
        de_results(NULL)
      })
    })
  })
  
  # --- Lancer l'analyse DE (Cluster vs Cluster) ---
  observeEvent(input$run_cluster_vs_cluster_markers, {
    obj <- visual_seurat_obj()
    req(obj)
    
    if (!"seurat_clusters" %in% colnames(obj@meta.data) || !is.factor(obj@meta.data$seurat_clusters)) {
      showNotification("Le clustering (colonne 'seurat_clusters') n'est pas disponible ou est mal formaté. Veuillez clusteriser d'abord.", type = "error")
      return(NULL)
    }
    req(input$de_cluster1, input$de_cluster2)
    
    if (input$de_cluster1 == input$de_cluster2) {
      showNotification("Veuillez sélectionner deux clusters différents pour la comparaison.", type = "warning")
      return(NULL)
    }
    
    withProgress(message = paste0("Calcul des marqueurs: Cluster ", input$de_cluster1, " vs Cluster ", input$de_cluster2, "..."), value = 0, {
      
      # Set default identity to seurat_clusters
      DefaultAssay(obj) <- "RNA"
      obj <- SetIdent(obj, value = "seurat_clusters")
      
      incProgress(0.3, detail = "Exécution de FindMarkers...")
      tryCatch({
        markers <- FindMarkers(
          object = obj,
          assay = "RNA",
          ident.1 = input$de_cluster1,
          ident.2 = input$de_cluster2,
          min.pct = input$de_min_pct,
          logfc.threshold = input$de_logfc_threshold,
          test.use = input$de_method
        )
        
        # Add gene names as a column and calculate Pct_diff
        markers <- markers %>%
          tibble::rownames_to_column(var = "gene") %>%
          dplyr::mutate(pct_diff = abs(pct.1 - pct.2))
        
        # Add a dummy 'cluster' column for compatibility with FindAllMarkers output structure
        markers$cluster <- paste0(input$de_cluster1, "_vs_", input$de_cluster2)
        
        de_results(markers)
        showNotification(paste0("Analyse FindMarkers terminée: Cluster ", input$de_cluster1, " vs Cluster ", input$de_cluster2, "."), type = "message")
      }, error = function(e) {
        showNotification(paste("Erreur lors de FindMarkers:", e$message), type = "error")
        de_results(NULL)
      })
    })
  })
  
  # --- Traitement et filtrage des résultats DE pour l'affichage ---
  filtered_de_results <- reactive({
    df <- de_results()
    req(df)
    
    # Filtrer par p-value ajustée
    df <- df %>%
      dplyr::filter(p_val_adj < input$de_pvalue_cutoff)
    
    # Filtrer par log2FC (en prenant la valeur absolue pour les comparaisons bilatérales)
    df <- df %>%
      dplyr::filter(abs(avg_log2FC) >= input$de_display_logfc_cutoff)
    
    # Trier les résultats
    sort_column <- input$de_sort_by
    if (!sort_column %in% colnames(df)) {
      showNotification(paste0("Colonne de tri '", sort_column, "' non trouvée. Tri par défaut (p_val_adj)."), type = "warning")
      sort_column <- "p_val_adj"
    }
    
    if (input$de_sort_desc) {
      df <- df %>% arrange(desc(!!sym(sort_column)))
    } else {
      df <- df %>% arrange(!!sym(sort_column))
    }
    
    return(df)
  })
  
  # --- Affichage du tableau de résultats DE ---
  output$de_results_table <- DT::renderDataTable({
    df <- filtered_de_results()
    req(df)
    DT::datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # --- Téléchargement des résultats DE ---
  output$download_de_results <- downloadHandler(
    filename = function() {
      paste("DE_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(filtered_de_results(), file, row.names = FALSE)
    }
  )
  
  # --- Plot Heatmap des Top Marqueurs ---
  output$de_top_markers_heatmap <- renderPlot({
    obj <- visual_seurat_obj()
    req(obj)
    df <- filtered_de_results()
    req(df)
    
    if (!"seurat_clusters" %in% colnames(obj@meta.data) || !is.factor(obj@meta.data$seurat_clusters)) {
      showNotification("La colonne 'seurat_clusters' n'est pas disponible ou est mal formatée dans l'objet Seurat. Le heatmap ne peut pas être généré.", type = "error")
      return(NULL)
    }
    
    top_n_markers <- input$de_num_top_markers
    
    if (nrow(df) == 0) {
      showNotification("Aucun marqueur significatif trouvé avec les filtres actuels pour la heatmap.", type = "warning")
      return(NULL)
    }
    
    # Sélectionner les N top marqueurs par cluster si FindAllMarkers,
    # ou les N top marqueurs globaux si Cluster vs Cluster
    if ("cluster" %in% colnames(df)) {
      # This means FindAllMarkers or Cluster vs Cluster with 'cluster' col added
      # For FindAllMarkers, take top N per cluster
      if (length(unique(df$cluster)) > 1) { # Implies FindAllMarkers
        top_genes <- df %>%
          group_by(cluster) %>%
          top_n(top_n_markers, wt = avg_log2FC) %>% # Use avg_log2FC as weight for top genes
          pull(gene) %>%
          unique()
      } else { # Implies Cluster vs Cluster
        top_genes <- df %>%
          top_n(top_n_markers, wt = avg_log2FC) %>% # Use avg_log2FC as weight for top genes
          pull(gene) %>%
          unique()
      }
    } else {
      # Should not happen with current logic, but fallback for safety
      top_genes <- df %>%
        top_n(top_n_markers, wt = avg_log2FC) %>%
        pull(gene) %>%
        unique()
    }
    
    # Ensure selected genes are actually present in the Seurat object
    genes_in_object <- intersect(top_genes, rownames(obj))
    if (length(genes_in_object) == 0) {
      showNotification("Aucun des top marqueurs sélectionnés n'est présent dans l'objet Seurat pour la heatmap.", type = "warning")
      return(NULL)
    }
    
    # Check if there are enough features for DoHeatmap
    if (length(genes_in_object) < 2) { # DoHeatmap needs at least 2 features to draw a heatmap
      showNotification("Pas assez de gènes (moins de 2) pour générer une heatmap. Veuillez ajuster les filtres DE ou le nombre de top marqueurs.", type = "warning")
      return(NULL)
    }
    
    # Re-order cells by cluster for better visualization in heatmap
    cell_order <- order(obj$seurat_clusters)
    
    plot <- DoHeatmap(obj, features = genes_in_object, group.by = "seurat_clusters", cells = colnames(obj)[cell_order]) +
      ggtitle(paste0("Top ", length(genes_in_object), " Marqueurs d'Expression Différentielle")) +
      theme(axis.text.y = element_text(size = 8))
    
    print(plot) # Must use print() for plotOutput
  })
  
  # Observe the button to trigger heatmap plot
  observeEvent(input$plot_top_markers, {
    # This button triggers the renderPlot directly
  })
  
  # ---------------------------------------------------------------------
  # Logique de l'onglet "Modification" (Simplifié)
  # ---------------------------------------------------------------------
  
  # Mettre à jour les choix de clusters disponibles pour l'exclusion (UI)
  observeEvent(visual_seurat_obj(), {
    obj <- visual_seurat_obj()
    req(obj)
    if (!"seurat_clusters" %in% colnames(obj@meta.data) || !is.factor(obj@meta.data$seurat_clusters)) {
      updateSelectInput(session, "subset_clusters_to_exclude", choices = c("Clustering non disponible ou mal formaté" = ""))
      return(NULL)
    }
    
    cluster_levels <- levels(obj$seurat_clusters)
    # Au lieu de tout sélectionner par défaut, on ne sélectionne rien (aucun cluster à exclure initialement)
    updateSelectInput(session, "subset_clusters_to_exclude", choices = cluster_levels, selected = NULL)
  })
  
  # --- Création de l'objet subset ---
  observeEvent(input$create_subset_obj, {
    obj <- visual_seurat_obj()
    req(obj)
    req(input$subset_name) # On requiert le nom du subset
    
    # input$subset_clusters_to_exclude peut être vide si l'utilisateur ne veut rien exclure
    # C'est géré par la logique suivante.
    
    if (trimws(input$subset_name) == "") {
      showNotification("Veuillez donner un nom à l'objet subset.", type = "warning")
      return(NULL)
    }
    
    if (!"seurat_clusters" %in% colnames(obj@meta.data) || !is.factor(obj@meta.data$seurat_clusters)) {
      showNotification("La colonne 'seurat_clusters' n'est pas disponible ou est mal formatée dans l'objet principal. Impossible de subsetter par cluster.", type = "error")
      subset_seurat_obj(NULL)
      return(NULL)
    }
    
    # Intersecter les clusters sélectionnés pour exclusion avec les niveaux réels de l'objet
    clusters_to_exclude <- intersect(input$subset_clusters_to_exclude, levels(obj$seurat_clusters))
    
    if (length(clusters_to_exclude) == length(levels(obj$seurat_clusters))) {
      showNotification("Attention : Tous les clusters sont sélectionnés pour exclusion. Le subset résultant sera vide.", type = "warning")
      subset_seurat_obj(NULL)
      output$subset_summary <- renderPrint({
        cat("L'objet subset créé est vide car tous les clusters ont été exclus.\n")
      })
      return(NULL)
    }
    
    withProgress(message = "Création de l'objet subset...", value = 0.5, {
      tryCatch({
        DefaultAssay(obj) <- "RNA"
        obj$seurat_clusters <- as.factor(obj$seurat_clusters)
        obj <- SetIdent(obj, value = "seurat_clusters")
        
        # --- Début du débogage avancé pour le problème "'arg' should be one of..." ---
        print("--- Début du Debugging Subset ---")
        print(paste("Nom de l'objet Seurat (dans create_subset_obj):", deparse(substitute(obj))))
        print(paste("Classe de l'objet Seurat:", class(obj)))
        print(paste("Nombre de cellules (avant subset):", ncol(obj)))
        print(paste("Nombre de gènes (avant subset):", nrow(obj)))
        
        current_ident <- Idents(obj)
        print(paste("Identité active avant subset:", if (is.null(current_ident)) "NULL" else paste(unique(current_ident)[1:min(5, length(unique(current_ident)))], collapse = ", "), "..."))
        
        if ("seurat_clusters" %in% colnames(obj@meta.data)) {
          print(paste("Classe de obj$seurat_clusters dans metadata:", class(obj@meta.data$seurat_clusters)))
          print(paste("obj$seurat_clusters est un facteur? ", is.factor(obj@meta.data$seurat_clusters)))
          print(paste("Nombre de niveaux de seurat_clusters:", length(levels(obj@meta.data$seurat_clusters))))
          print("Niveaux de seurat_clusters:")
          print(levels(obj@meta.data$seurat_clusters))
          print("Table des seurat_clusters (premières 10 entrées):")
          print(head(table(obj@meta.data$seurat_clusters), 10))
        } else {
          print("La colonne 'seurat_clusters' n'existe pas dans obj@meta.data (debugging print).")
        }
        
        print(paste("Assay par défaut avant subset:", DefaultAssay(obj)))
        print(paste("Clusters à exclure (clusters_to_exclude):", paste(clusters_to_exclude, collapse = ", ")))
        print(paste("Longueur de clusters_to_exclude:", length(clusters_to_exclude)))
        print(paste("Type de clusters_to_exclude:", typeof(clusters_to_exclude)))
        
        if (is.factor(obj@meta.data$seurat_clusters)) {
          invalid_to_exclude <- setdiff(clusters_to_exclude, levels(obj@meta.data$seurat_clusters))
          if (length(invalid_to_exclude) > 0) {
            print(paste("ATTENTION (debugging): Clusters sélectionnés pour exclusion non trouvés dans les niveaux existants:", paste(invalid_to_exclude, collapse = ", ")))
          }
        }
        print("Tentative de subset...")
        # --- Fin du débogage avancé ---
        
        
        if (length(clusters_to_exclude) == 0) {
          # Si aucun cluster n'est sélectionné pour exclusion, l'objet subset est l'objet original
          subset_obj <- obj
          showNotification("Aucun cluster n'a été sélectionné pour exclusion. L'objet subset est une copie de l'objet principal.", type = "message") # CORRECTION: type = "info" -> "message"
        } else {
          # Sinon, on exclut les clusters sélectionnés
          subset_obj <- subset(obj, idents = clusters_to_exclude, invert = TRUE)
          showNotification(paste0("Clusters ", paste(clusters_to_exclude, collapse = ", "), " exclus du subset."), type = "message")
        }
        
        if (ncol(subset_obj) < 1) {
          showNotification("Le subset créé est vide. Aucun cellule n'a été conservée avec les clusters restants.", type = "error")
          subset_seurat_obj(NULL)
          output$subset_summary <- renderPrint({
            cat("Aucun subset créé (vide).\n")
          })
          return(NULL)
        }
        
        subset_seurat_obj(subset_obj)
        showNotification(paste0("Objet subset '", input$subset_name, "' créé avec succès. Contient ", ncol(subset_obj), " cellules."), type = "message") # CORRECTION: type = "success" -> "message"
        
        output$subset_summary <- renderPrint({
          cat(paste0("Résumé de l'objet subset: '", input$subset_name, "'\n"))
          print(subset_obj)
          if (length(clusters_to_exclude) > 0) {
            cat("\nClusters exclus:\n")
            print(paste(clusters_to_exclude, collapse = ", "))
          }
          cat("\nClusters conservés:\n")
          print(table(subset_obj$seurat_clusters))
        })
      }, error = function(e) {
        showNotification(paste("Erreur interne lors de la création du subset. Message: ", e$message), type = "error")
        subset_seurat_obj(NULL)
      })
    })
  })
  
  # --- Téléchargement du subset en .rds ---
  output$download_subset_rds <- downloadHandler(
    filename = function() {
      # Assurez-vous que le nom de fichier est valide et se termine par .rds
      req(subset_seurat_obj()) # CRITICAL: Assure que l'objet subset existe
      req(input$subset_name)   # Assure qu'un nom a été donné
      
      clean_name <- gsub("[^A-Za-z0-9_.-]", "", input$subset_name) # Supprime les caractères spéciaux
      paste0(clean_name, "_", Sys.Date(), ".rds")
    },
    content = function(file) {
      obj <- subset_seurat_obj()
      req(obj) # CRITICAL: Répète la vérification dans `content` pour la robustesse du téléchargement
      
      withProgress(message = "Exportation du subset en .rds...", value = 0.5, {
        tryCatch({
          saveRDS(obj, file = file)
          showNotification("Subset exporté avec succès en .rds.", type = "message") # CORRECTION: type = "success" -> "message" (bien que "success" n'existe pas, celle-ci était une ajout pour la robustesse, la correction standard s'applique)
        }, error = function(e) {
          showNotification(paste("Erreur lors de l'exportation du subset:", e$message), type = "error")
          # En cas d'erreur, ne pas écrire le fichier ou écrire un fichier vide pour éviter le HTML.
          # Le tryCatch avec showNotification suffit généralement.
        })
      })
    }
  )
  
} # Fin du serveur

shinyApp(ui = ui, server = server)