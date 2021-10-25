
##PLEASE CHANGE THE WORKING DIRECTORY###

setwd("/data/2623287c/Project1/upload_app")

source("tcell_libs_raw_dash.R", local = TRUE)

##Creating the ui object for the user interface

ui <- tagList(
  dashboardPage(title = "Upload Data Single Cell Analysis",
                dashboardHeader(title = "Upload Data Single Cell Analysis",
                                tags$li(class = "dropdown",
                                        tags$style(".main-header {max-height: 90px}"),
                                        tags$style(".main-header .logo {height: 90px}"),
                                        tags$a(href='https://www.gla.ac.uk/', target="_blank",
                                               tags$img(src='uog_logo.png', height = 60, width = 180))
                                ),
                                titleWidth = 380),
                dashboardSidebar(
                  tags$style(".left-side, .main-sidebar {padding-top: 90px}"),
                  width = 250,
                  sidebarMenu(
                    id = 'sidebar',
                    menuItem("About", tabName = "About", icon = icon("door-open")),
                    menuItem("Upload Data", tabName = "add_data", icon = icon("folder-plus")),
                    uiOutput("men_cluster"),
                    uiOutput("men_de"),
                    uiOutput("men_pseudotime")
                  )
                ),
                dashboardBody(
                  tags$head(
                    tags$link(rel = "stylesheet", type = "text/css", href = "tcell_amp_custom_dash.css"),
                    includeHTML("tcell_js.htm")
                  ),
                  # Boxes need to be put in a row (or column)
                  tabItems(
                    # First tab content
                    # source(file.path("tcell_introduction_page_dash.R"), local = TRUE)$value,
                    tabItem(tabName = "About",
                            source(file.path("tcell_introduction_page_dash.R"), local = TRUE)$value),
                    tabItem(tabName = "add_data", useShinyjs(),
                            tabsetPanel(
                              id = "hidden_tabs",
                              type = "hidden",
                              tabPanelBody("PAN_UPLOAD",
                                           h1("Upload Data"),
                                           tabsetPanel(
                                             id ="hidden_tabs_upload",
                                             type = "hidden",
                                             tabPanelBody("panel1",
                                                          wellPanel(
                                                            tags$b("Please select a uniprot file"), 
                                                            fluidRow(
                                                              column(
                                                                2,
                                                                actionButton("upload_uniprot_file",label = "Browse") #from: https://stackoverflow.com/questions/51191701/r-shiny-fileinput-large-files
                                                                
                                                              ),
                                                              column(
                                                                8,
                                                                textOutput("uniprot_filechosen")
                                                              )
                                                            ),
                                                            hr(),
                                                            shinyjs::hidden(checkboxGroupInput("radio_uniprot_columns", "Choose gene name and uniprot link columns (assumes gene names is first - you may need to edit your choices", choices = NULL, inline = TRUE)),
                                                            shinyjs::hidden(actionButton("choose_uniprot_columns", label = "Click to choose the uniprot columns (Gene Name and Link"))
                                                          ),
                                                          
                                                          wellPanel(
                                                            h4("Seurat files uploaded:"),
                                                            textOutput("files_so_far_text")
                                                          ),
                                                          
                                                          wellPanel(
                                                            tags$b("Please select a Seurat file"), 
                                                            fluidRow(
                                                              column(
                                                                2,
                                                                actionButton("filechoose",label = "Browse") #from: https://stackoverflow.com/questions/51191701/r-shiny-fileinput-large-files
                                                                
                                                              ),
                                                              column(
                                                                8,
                                                                textOutput("filechosen")
                                                              )
                                                            ),
                                                            hr(),
                                                            
                                                            
                                                            textInput("sample_text", "Sample name:", value = "", placeholder = "i.e. WT1 or KO1"),
                                                            textInput("group_text", "Group:", value = "", placeholder = "i.e. WT or KO"),
                                                            checkboxInput("check_addmtfeat", "Add percent.mt and n_feature RNA now:", value = F),
                                                            checkboxInput("check_addnumberofcells", "Add if you want to sample the cells:", value = F)
                                                            
                                                          )
                                                          
                                             ),
                                             tabPanelBody("panel2",
                                                          wellPanel(h4("The plots will appear here:"),
                                                                    withSpinner(plotOutput("nfeat_plot")),
                                                                    withSpinner(plotOutput("percent_plot"))
                                                          )),
                                             conditionalPanel(
                                               condition = "input.check_addmtfeat == 1",
                                               wellPanel(
                                                 numericInput("for_mt", "percent.mt:", 1, min = 0, max = 100),
                                                 numericInput("for_featmin", "n_feature RNA minimum", 1, min = 0, max = 1000),
                                                 numericInput("for_featmax", "n_feature RNA maximum", 1, min = 0, max = 5000),
                                                 shinyjs::hidden(actionButton("add_mtandnfeat", "Click to subset"))
                                               )
                                             ),
                                             conditionalPanel(
                                               condition = "input.check_addnumberofcells == 1",
                                               wellPanel(
                                                 numericInput("for_samp", "Sample number:", 100, min = 0, max = 6000)
                                               )
                                             ),
                                             
                                             wellPanel(
                                               fluidRow(
                                                 column(
                                                   4,
                                                   actionButton("upload_button", "Upload")),
                                                 column(
                                                   4,
                                                   shinyjs::disabled(actionButton("add_ano_button", "Add another"))),
                                                 column(
                                                   4,
                                                   shinyjs::disabled(actionButton("fin_button", "Finished")))
                                               )
                                             )
                                           )
                                           
                              ),
                              tabPanelBody("PAN_INTEGRATE", 
                                           wellPanel(
                                             h1("Integrate")
                                           ),
                                           
                                           # wellPanel(
                                           #   h4("Seurat files to integrate:"),
                                           #   textOutput("files_so_far_text")
                                           # ),
                                           wellPanel(
                                             actionButton("but_int", "Integrate")
                                           )
                              ),
                              
                              
                              tabPanelBody("PAN_DIM", 
                                           wellPanel(
                                             h1("Choose a Dimension")
                                           ),
                                           tabsetPanel(
                                             id = "hidden_tabs_dim",
                                             type = "hidden",
                                             tabPanelBody("pan_dim_1",
                                                          wellPanel(
                                                            checkboxInput("check_adddim", "Add dim now", value = F)
                                                          )
                                             ),
                                             tabPanelBody("pan_dim_2",
                                                          wellPanel(
                                                            withSpinner(plotOutput("elbow_plot"))
                                                          )
                                             ),
                                             
                                             conditionalPanel(
                                               condition = "input.check_adddim == 1",
                                               wellPanel(
                                                 numericInput("for_dim", "Number of dimensions:", 1, min = 0, max = 100),
                                                 shinyjs::hidden(actionButton("but_dclust", "Click to add dimensions")),
                                                 shinyjs::hidden(actionButton("but_drep", "Click to add dimensions"))
                                                 
                                               )
                                             ),
                                             wellPanel(
                                               actionButton("but_dimclust", "Choose Dimensions"),
                                               shinyjs::hidden(actionButton("but_dimrep", "Choose Dimensions"))
                                             )
                                           )
                              ),
                              
                              tabPanelBody("PAN_SUBCLUST", 
                                           wellPanel(
                                             h1("Subset the clusters")
                                           ),
                                           tabsetPanel(
                                             id = "hidden_tabs_subclust",
                                             type = "hidden",
                                             tabPanelBody("pan_clus_1",
                                                          wellPanel(
                                                            checkboxInput("check_clus", "Add clusters to keep now", value = F)
                                                          )
                                             ),
                                             tabPanelBody("pan_clus_2",
                                                          withSpinner(plotOutput("umi_plot")),
                                                          withSpinner(plotOutput("gene_plot")),
                                                          withSpinner(plotOutput("dim_plot"))
                                             ),
                                             
                                             
                                             conditionalPanel(
                                               condition = "input.check_clus == 1",
                                               wellPanel(
                                                 checkboxGroupInput("for_clust", "Choose clusters to keep:", choices = c(0,1,2,3,4))),
                                               shinyjs::hidden(actionButton("but_clus", "Click to add clusters"))
                                             ),
                                             wellPanel(
                                               actionButton("but_sub", "Subset clusters"),
                                               shinyjs::hidden(actionButton("but_cont", "Continue")),
                                               shinyjs::hidden(actionButton("but_rep", "Repeat dimension and subset"))
                                             )
                                           )
                              ),
                              tabPanelBody("PAN_LOOPS", 
                                           wellPanel(
                                             h1("Run Analysis Code"),
                                             h4("The code is going to run to produce different objects for the different resolutions.
                                        This will take a while to run."),
                                             h4("Please choose the parameters you wish to use. 
                                        The default values will be used if nothing selected. This will impact the results."),
                                             br(),
                                             h4("Click to change the parameters for each section:"),
                                             
                                             checkboxInput("check_res", "Click to add resolutions (default: 0.15, 0.55):", value = T),
                                             conditionalPanel(condition = "input.check_res == 1",
                                                              h5("Please select start and end resolutions for the loops. 
                                                              The resolution loops will be run in increments of 0.1"),
                                                              sliderInput(inputId = "get_res", label = strong("Louvain algorithm resolution"), value = c(0.1, 0.5), min = 0.1, max = 0.5, step = 0.1, round = F, post = "5"),
                                                              # numericInput("get_minres", "Minimum resolution:", 0.1, min = 0.05, max = 0.85, value = 0.15),
                                                              # numericInput("get_maxres", "Maximum resolution:", 0.1, min = 0.05, max = 0.85, value = 0.55)
                                             ),
                                             checkboxInput("check_slingshot", strong("Slingshot"), value = F),
                                             conditionalPanel(condition = "input.check_slingshot == 1",
                                                              wellPanel(
                                                                h5("Click to choose the parameters you wish to add now and change the values"),
                                                                checkboxInput("check_slingstart", "Start Cluster(s):", value = F),
                                                                conditionalPanel(condition = "input.check_slingstart == 1",
                                                                                 h5("Please select start cluster(s) you want to use for slingshot"),
                                                                                 selectizeInput(inputId = "get_slingstart", label = strong("Choose start cluster(s):"),choices = NULL, multiple = T)
                                                                                 # checkboxGroupInput("get_slingstart", "Choose start cluster:", choices = levels(combined$seurat_clusters), inline = TRUE)  
                                                                ),
                                                                checkboxInput("check_slingend", "End Cluster(s)", value = F),
                                                                conditionalPanel(condition = "input.check_slingend == 1",
                                                                                 h5("Please select end cluster(s) you want to use for slingshot"),
                                                                                 selectizeInput(inputId = "get_slingend", label = strong("Choose end cluster(s):"),choices = NULL, multiple = T)
                                                                                 # checkboxGroupInput("get_slingend", "Choose end cluster:", choices = levels(combined$seurat_clusters), inline = TRUE)  
                                                                ))),
                                             
                                             checkboxInput("check_tradeseq", strong("tradeSeq"), value = F),
                                             conditionalPanel(condition = "input.check_tradeseq == 1",
                                                              wellPanel(
                                                                h5("Click to choose the parameters you wish to add now and change the values"),
                                                                checkboxInput("check_tradeknots", "Number of Knots", value = F),
                                                                conditionalPanel(condition = "input.check_tradeknots == 1",
                                                                                 h5("Please select the number of knots you want to use for tradeSeq"),
                                                                                 numericInput("get_tradenknots", "Number of knots:", 1, min = 2, max = 8, value = 3)
                                                                )
                                                              )),
                                             
                                             checkboxInput("check_mon2", strong("Monocle 2"), value = F),
                                             conditionalPanel(condition = "input.check_mon2 == 1",
                                                              wellPanel(
                                                                h5("Click to choose the parameters you wish to add now and change the values"),
                                                                checkboxInput("check_mon2start", "Start cluster", value = F),
                                                                conditionalPanel(condition = "input.check_mon2start == 1",
                                                                                 h5("Please select the starting cluster to use for monocle 2"),
                                                                                 selectizeInput(inputId = "get_mon2start", label = strong("Choose start cluster:"),choices = NULL, multiple = F)
                                                                                 # checkboxGroupInput("get_mon2start", "Choose start cluster:", choices = 1:levels(combined$seurat_clusters), inline = TRUE)
                                                                )
                                                              )),
                                             
                                             checkboxInput("check_mon3", strong("Monocle 3"), value = F),
                                             conditionalPanel(condition = "input.check_mon3 == 1",
                                                              wellPanel(
                                                                checkboxInput("check_mon3start", "Click to choose to select your own root cluster", value = F),
                                                                conditionalPanel(condition = "input.check_mon3start == 1",
                                                                                 h5("Please select start cluster(s) you want to use for monocle 3"),
                                                                                 selectizeInput(inputId = "get_mon3start", label = strong("Choose start cluster(s):"),choices = NULL, multiple = T)
                                                                                 # checkboxGroupInput("get_mon3start", "Choose start cluster:", choices = 0:15, inline = TRUE)  
                                                                )
                                                              )
                                                              
                                                              
                                             )
                                           ),
                                           actionButton("but_loop", "Run Code")
                              ),
                              tabPanelBody("PAN_WHENRUNNING",
                                           source(file.path("upload_loops_page.R"), local = TRUE)$value),
                              tabPanelBody("PAN_AFTERRUN",
                                           h1("Completed!"),
                                           h3("Explore the tabs to see the single cell and pseudotime analysis"),
                                           actionButton("but_analysis", "View Analysis Tabs")
                              )
                            )
                    ),
                    tabItem(tabName = "all_cluster_res",
                            
                            wellPanel(fluidRow(
                              
                              column(
                                12,
                                tabBox(
                                  title = "UMAP",
                                  # The id lets us use input$tabset1 on the server to find the current tab
                                  #id = "vis_de", 
                                  #height = "250px",
                                  width = 12,
                                  tabPanel("All cells", 
                                           wellPanel(
                                             # sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res1, min = res1, max = res2, step = diff_res, round = F),
                                             sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = 0.1, min = 0.1, max = 0.5, step = 0.1, round = F, post = "5"),
                                             withSpinner(plotOutput("labelled_umap", width = "600px", height = "400px")),
                                             wellPanel(
                                               h4("Download specifications"),
                                               flowLayout(numericInput("lumap_height", "Plot height (cm):", value = 14),
                                                          numericInput("lumap_width", "Plot width (cm):", value = 14),
                                                          radioButtons("lumap_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                               downloadButton('dwnl_lumap','Download Plot')
                                             ),
                                             # )
                                             
                                             tags$hr(),
                                             uiOutput("cluster_annot"),
                                             wellPanel(style = "background:#385A4F",
                                                       tags$hr(),
                                                       tags$p(style = "font-family:Arial;color:white",
                                                              paste("Adjusting resolution (0.15, 0.25, 0.35, 0.45 or 0.55) to set the clustering granularity. This option allows the subdivision of clusters further into sub-populations and the subsequent labelling and interrogation of differential expression.")
                                                              
                                                              
                                                       )
                                             )
                                             
                                           )
                                  ),
                                  tabPanel("Sample comparison", withSpinner(plotOutput("all_groups", height = "400px")),
                                           wellPanel(
                                             h4("Download specifications"),
                                             flowLayout(numericInput("grps_height", "Plot height (cm):", value = 7),
                                                        numericInput("grps_width", "Plot width (cm):", value = 30),
                                                        radioButtons("grps_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                             downloadButton('dwnl_grps','Download Plot')
                                           )
                                  ), 
                                  tabPanel("Group comparison", withSpinner(plotOutput("splitby_group", height = "400px")),
                                           wellPanel(
                                             h4("Download specifications"),
                                             flowLayout(numericInput("grps_height", "Plot height (cm):", value = 7),
                                                        numericInput("grps_width", "Plot width (cm):", value = 30),
                                                        radioButtons("grps_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                             downloadButton('dwnl_splitby_groups', 'Download Plot')
                                           )
                                  )
                                  
                                )
                                
                              )
                              
                              
                            )
                            )
                    ),
                    tabItem(tabName = "grps_cluster_res",
                            
                            wellPanel(fluidRow(
                              
                              column(
                                12, 
                                tabBox(
                                  title = "Cluster markers",
                                  # The id lets us use input$tabset1 on the server to find the current tab
                                  id = "vis_de", 
                                  #height = "250px",
                                  width = 12,
                                  tabPanel("Marker table", 
                                           wellPanel(
                                             # box(
                                             # width = NULL,
                                             # solidHeader = TRUE,
                                             uiOutput("dyn_clusters"),
                                             # uiOutput("topclgenes"),
                                             withSpinner(DTOutput("top_conserved_genes"))
                                           )
                                           
                                           # )
                                  ),
                                  tabPanel("Marker feature plots", 
                                           uiOutput("top_markers_umap"),
                                           withSpinner(plotOutput("conserved_markers_umap")),
                                           wellPanel(
                                             h4("Download specifications"),
                                             flowLayout(numericInput("markers_height", "Plot height (cm):", value = 20),
                                                        numericInput("markers_width", "Plot width (cm):", value = 20),
                                                        radioButtons("markers_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                             downloadButton('dwnl_markers','Download Plot')
                                           )
                                  )
                                )
                                
                                
                              )
                            ),
                            uiOutput("box_2_2")
                            )
                            
                            
                    ),
                    tabItem(tabName = "ge_vis_gene",
                            #fluidRow(
                            wellPanel(
                              fluidRow(
                                column(
                                  12,
                                  h4("Single gene DE visualization"),
                                  selectizeInput(inputId = "de_genes", label = strong("Choose gene:"),choices = NULL, multiple = F)
                                )
                              ),
                              fluidRow(
                                tabsetPanel(
                                  tabPanel("Violin plot", withSpinner(plotOutput("de_stim_vs_ctrl_vp")),
                                           wellPanel(
                                             h4("Download specifications"),
                                             flowLayout(numericInput("violp_height", "Plot height (cm):", value = 12),
                                                        numericInput("violp_width", "Plot width (cm):", value = 30),
                                                        radioButtons("violp_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                             downloadButton('dwnl_violp','Download Plot')
                                           )
                                  ),
                                  tabPanel("UMAP feature plot", withSpinner(plotOutput("de_stim_vs_ctrl_um")),
                                           wellPanel(
                                             h4("Download specifications"),
                                             flowLayout(numericInput("featurep_height", "Plot height (cm):", value = 10),
                                                        numericInput("featurep_width", "Plot width (cm):", value = 30),
                                                        radioButtons("featurep_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                             downloadButton('dwnl_featurep','Download Plot')
                                           )
                                  )
                                ),
                                
                                
                              ),
                              tags$hr(),
                              uiOutput("box_1_1")
                              
                            )
                    ),
                    tabItem(tabName = "m_ge_vis_gene",
                            fluidRow(
                              wellPanel(
                                h4("Dotplot for multiple gene DE visualization"),
                                selectizeInput(inputId = "select_markers_dotplot", label = strong("Choose gene:"), choices = NULL, multiple = T),
                                withSpinner(plotOutput("marker_dotplot")),
                                wellPanel(
                                  h4("Download specifications"),
                                  flowLayout(numericInput("dotp_height", "Plot height (cm):", value = 10),
                                             numericInput("dotp_width", "Plot width (cm):", value = 30),
                                             radioButtons("dotp_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                  downloadButton('dwnl_dotp','Download Plot')
                                ),
                                tags$hr(),
                                uiOutput("box_1_2")
                                
                                # )
                              )
                            )
                    ),
                    tabItem(tabName = "ge_vis_cell",
                            
                            wellPanel(
                              
                              fluidRow(
                                column(
                                  4,
                                  uiOutput("cluster_ids"),
                                  ##Comparing conditions new addition
                                  selectInput(inputId = "ra_conds", label = strong("Choose conditions to compare:"),choices = conds)
                                )
                              ),
                              fluidRow(
                                
                                column(6,
                                       wellPanel(
                                         h4("DE scatterplot"),
                                         fluidRow(
                                           withSpinner(plotOutput("cell_type_plot", click = clickOpts(id ="plot_click"))),
                                           dataTableOutput("click_info"),
                                           wellPanel(
                                             h4("Download specifications"),
                                             flowLayout(numericInput("scatter_height", "Plot height (cm):", value = 14),
                                                        numericInput("scatter_width", "Plot width (cm):", value = 14),
                                                        radioButtons("scatter_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                             downloadButton('dwnl_scatter','Download Plot')
                                           ),
                                           
                                           
                                           
                                         ),
                                         tags$hr(),
                                         uiOutput("box_1_3a")
                                         
                                       )
                                ),
                                
                                column(
                                  6,
                                  wellPanel(
                                    h4("DE table"),
                                    
                                    # sliderInput(inputId = "top_genes", label = strong("Number of top DE genes:"), value = 100, min = 1, max = dim(cluster)[1], step = 1),
                                    # uiOutput("topdegenes"),
                                    withSpinner(DTOutput("top_de_genes")),
                                    tags$hr(),
                                    uiOutput("box_1_3b")
                                    
                                    # )
                                  )
                                )
                              )
                              
                            )
                            
                    ),
                    tabItem(tabName = "slingshot_tab", 
                            h1("Slingshot"),
                            wellPanel(
                              h3("Click the tabs to explore the different pseudotime plots"),
                              tabsetPanel(
                                tabPanel("UMAP",
                                         withSpinner(plotOutput("sling_UMAP_plot")),
                                         withSpinner(plotOutput("sling_UMAP_leg")),
                                         wellPanel(
                                           h4("Download specifications"),
                                           flowLayout(numericInput("slingumap_height", "Plot height (cm):", value = 14),
                                                      numericInput("slingumap_width", "Plot width (cm):", value = 14),
                                                      radioButtons("slingumap_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                           downloadButton('dwnl_slingumap','Download Plot')
                                         )),                                      
                                tabPanel("Phate",
                                         withSpinner(plotOutput("sling_PHATE_plot")),
                                         withSpinner(plotOutput("sling_PHATE_leg")),
                                         wellPanel(
                                           h4("Download specifications"),
                                           flowLayout(numericInput("slingphate_height", "Plot height (cm):", value = 20),
                                                      numericInput("slingphate_width", "Plot width (cm):", value = 20),
                                                      radioButtons("slingphate_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                           downloadButton('dwnl_slingphate','Download Plot')
                                         )),
                                tabPanel("Heatmap",
                                         wellPanel(
                                           radioButtons("select_sling_reduc", "Choose reduction:", c("UMAP" = "UMAP", "PHATE" = "PHATE"),
                                                        selected = "UMAP", inline = TRUE),
                                           withSpinner(plotOutput("sling_HEAT_PLOT")),
                                           withSpinner(DTOutput("sling_heat_info"))
                                         ),
                                         wellPanel(
                                           h4("Download specifications"),
                                           flowLayout(numericInput("slingheat_height", "Plot height (cm):", value = 20),
                                                      numericInput("slingheat_width", "Plot width (cm):", value = 20),
                                                      radioButtons("slingheat_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                           downloadButton('dwnl_slingheat','Download Plot')
                                         )) 
                              )),
                            uiOutput("box_3_1")
                    ),
                    tabItem(tabName = "trade_tab",                            
                            h1("tradeSeq Downstream of Slingshot"),
                            wellPanel(
                              h3("Click the tabs to explore the different pseudotime plots"),
                              tabsetPanel(
                                tabPanel("UMAP",
                                         withSpinner(plotOutput("trade_UMAP_plot")),
                                         wellPanel(
                                           h4("Download specifications"),
                                           flowLayout(numericInput("tradeumap_height", "Plot height (cm):", value = 20),
                                                      numericInput("tradeumap_width", "Plot width (cm):", value = 20),
                                                      radioButtons("tradeumap_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                           downloadButton('dwnl_tradeumap','Download Plot')
                                         )),                                      
                                tabPanel("Phate",
                                         withSpinner(plotOutput("trade_PHATE_plot")),
                                         wellPanel(
                                           h4("Download specifications"),
                                           flowLayout(numericInput("tradephate_height", "Plot height (cm):", value = 20),
                                                      numericInput("tradephate_width", "Plot width (cm):", value = 20),
                                                      radioButtons("tradephate_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                           downloadButton('dwnl_tradephate','Download Plot')
                                         )),
                                tabPanel("Gene Temporal Expression",
                                         radioButtons("select_gene_reduc", "Choose reduction:", c("UMAP" = "UMAP", "PHATE" = "PHATE"),
                                                      selected = "UMAP", inline = TRUE),
                                         selectInput(inputId = "gene_trade", label = strong("Choose gene:"), choices = NULL, multiple = F),
                                         withSpinner(plotOutput("trade_GENE_plot")),
                                         withSpinner(plotOutput("trade_SMOOTH_plot"))),
                                tabPanel("Heatmap",
                                         h5("The heatmaps may take around 10 minutes to load"),
                                         wellPanel(
                                           radioButtons("select_trade_reduc", "Choose reduction:", c("UMAP" = "UMAP", "PHATE" = "PHATE"),
                                                        selected = "UMAP", inline = TRUE),
                                           withSpinner(plotOutput("trade_HEAT_PLOT")),
                                           withSpinner(DTOutput("trade_heat_info"))
                                         ),
                                         wellPanel(
                                           h4("Download specifications"),
                                           flowLayout(numericInput("tradeheat_height", "Plot height (cm):", value = 20),
                                                      numericInput("tradeheat_width", "Plot width (cm):", value = 20),
                                                      radioButtons("tradeheat_format", "File format:", choices = list("PNG" = ".png", "PDF" = ".pdf", "JPEG" = ".jpeg"), selected = ".png", inline = T)),
                                           downloadButton('dwnl_tradeheat','Download Plot')
                                         )) 
                              )),
                            uiOutput("box_3_2")
                    ),
                    tabItem(tabName = "mon2_tab", 
                            wellPanel(
                              h1("Click the tabs to explore the different pseudotime plots"),
                              tabsetPanel(
                                tabPanel("Trajectory",
                                         radioButtons("colour_by_mon2", "Colour plot by:", c("Clusters" = "seurat_clusters", "Pseudotime" = "Pseudotime"),
                                                      selected = "seurat_clusters", inline = TRUE),
                                         withSpinner(plotOutput("mon2_TRA_plot"))),                                      
                                tabPanel("Heatmap",
                                         withSpinner(plotOutput("mon2_HEAT_plot")),
                                         withSpinner(DTOutput("mon2_heat_info"))
                                )
                              ) 
                            ),
                            uiOutput("box_3_3")
                    ),
                    tabItem(tabName = "mon3_tab", 
                            wellPanel(
                              h1("Click the tabs to explore the different pseudotime plots"),
                              tabsetPanel(
                                tabPanel("UMAP",
                                         radioButtons("colour_by_mon", "Colour plot by:", c("Clusters" = "seurat_clusters", "Pseudotime" = "pseudotime"),
                                                      selected = "seurat_clusters", inline = TRUE),
                                         withSpinner(plotOutput("mon_UMAP_plot"))),                                      
                                tabPanel("Heatmap",
                                         withSpinner(plotOutput("mon3_HEAT_plot")),
                                         withSpinner(DTOutput("mon3_heat_info"))
                                )
                              ) 
                            ),
                            uiOutput("box_3_4")
                    )
                  )
                )
  ),
  includeHTML("tcell_footer_dash.htm")
)

# alveri = readRDS("alveri.rds")

##server function to compute the outputs
server = function(input, output, session) {
  #remove
  # updateTabsetPanel(session, "hidden_tabs", selected = "PAN_LOOPS")
  # updateSelectizeInput(session = session, inputId = 'get_slingstart', choices = levels(Combined$seurat_clusters), server = TRUE)
  # updateSelectizeInput(session = session, inputId = 'get_slingend', choices = levels(Combined$seurat_clusters), server = TRUE)
  # updateSelectizeInput(session = session, inputId = 'get_mon2start', choices = levels(Combined$seurat_clusters), server = TRUE)
  # updateSelectizeInput(session = session, inputId = 'get_mon3start', choices = levels(Combined$seurat_clusters), server = TRUE)
  
  
  #Selecting the Seurat object
  umap_clusters = reactive({
    if(input$but_analysis == 1){
      tcells_combined_umap_list_res_skinny = readRDS("tcells_combined_umap_list_res_skinny.rds")
      tcells_combined_umap_list_res_skinny[[(input$clusters_res * 10)]]
    }else{
      NULL
    }
  })
  
  
  ##Upload uniprot####
  uniprot_path <- reactiveVal(value = NULL)
  
  observeEvent(input$upload_uniprot_file, {
    tryPath <- tryCatch(
      file.choose(), error = function(e){e}
    )
    
    if(inherits(tryPath, "error")){
      uniprot_path(NULL)
    } else {
      uniprot_path(tryPath)
    }
    shinyjs::show("radio_uniprot_columns")
    shinyjs::show("choose_uniprot_columns")
    updateCheckboxGroupInput(session, "radio_uniprot_columns", choices = colnames(fread(uniprot_path(), stringsAsFactors = F)))
    
  })
  
  observeEvent(input$tcells_uniprot_file, {
    tryPath <- tryCatch(
      file.choose(), error = function(e){e}
    )
    
    if(inherits(tryPath, "error")){
      uniprot_path(NULL)
    } else {
      uniprot_path(tryPath)
    }
    shinyjs::show("radio_uniprot_columns")
    shinyjs::show("choose_uniprot_columns")
    updateCheckboxGroupInput(session, "radio_uniprot_columns", choices = colnames(fread(uniprot_path(), stringsAsFactors = F)))
    
  })
  
  output$uniprot_filechosen <- renderText({
    if(is.null(uniprot_path())){
      "No file selected"
    } else {
      uniprot_path()     
    }    
  })
  
  
  observeEvent(input$choose_uniprot_columns, {
    uniprot_info <<- fread(uniprot_path(), stringsAsFactors = F, select = input$radio_uniprot_columns) #reads only the necessary columns
    shinyjs::hide("radio_uniprot_columns")
    shinyjs::hide("choose_uniprot_columns")
    shinyjs::disable("upload_uniprot_file")
    output$uniprot_filechosen <- renderText({
      "File successfully uploaded"
    })
  })
  
  
  
  
  ##Upload data####
  i = 1 #for file number 
  object <<- vector(mode = "list")
  Combined <<- NULL
  
  path <- reactiveVal(value = NULL)
  
  observeEvent(input$filechoose, {
    tryPath <- tryCatch(
      file.choose(), error = function(e){e}
    )
    
    if(inherits(tryPath, "error")){
      path(NULL)
    } else {
      path(tryPath)
    }
    
  })
  
  output$filechosen <- renderText({
    if(is.null(path())){
      "No file selected"
    } else {
      path() 
    }
  })
  
  
  observeEvent(input$upload_button,{
    shinyjs::disable("sample_text")
    shinyjs::disable("group_text")
    shinyjs::disable("upload_button")
    shinyjs::disable("fin_button")
    shinyjs::disable("for_samp")
    shinyjs::disable("filechoose")

    if(input$check_addmtfeat == 1){ #if checkbox clicked
      shinyjs::disable("for_featmin")
      shinyjs::disable("for_featmax")
      shinyjs::disable("for_mt")
      
      object[[i]] = readRDS(path())
      object[[i]]$sample = input$sample_text 
      object[[i]]$group = input$group_text 
      object[[i]][["percent.mt"]] <- PercentageFeatureSet(object = object[[i]], pattern = "^mt-")
      
      l1 <<-input$for_featmin
      l2 <<-input$for_featmax
      l3 <<-input$for_mt
      t1 <<- input$for_samp
      
      if(input$check_addnumberofcells == 1){ #if checkbox clicked to sample genes 
        object[[i]] <- subset(object[[i]], cells = sample(Cells(object[[i]]), t1), subset = nFeature_RNA > l1 & nFeature_RNA < l2 & percent.mt < l3)
      }else{ #if not to sample genes 
        object[[i]] <- subset(object[[i]], subset = nFeature_RNA > l1 & nFeature_RNA < l2 & percent.mt < l3)
      }
      
      object[[i]][["percent.mt"]] <- PercentageFeatureSet(object = object[[i]], pattern = "^mt-")
      
      object[[i]] = SCTransform(object[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
      assign(paste("all_genes_",i), rownames(object[[i]])) #makes a new variable with the gene names
      print(object[[i]])

      object[[i]] <<- object[[i]]
      
      print(object)
      
      updateCheckboxInput(session, "check_addmtfeat", value =0)
      shinyjs::enable("add_ano_button")
      shinyjs::enable("fin_button")
      
      output$files_so_far_text = renderText({
        ptext = c()
        for(t in 1:length(object)){
          ptext[t] = paste(t, ": ", unique(object[[t]]$sample, " |"))
        }
        ptext #to print the seurat files loaded so far
      })
      
      i <<- i+1
      return()
    } 
    updateTabsetPanel(session, "hidden_tabs_upload", selected = "panel2")
    object[[i]] = readRDS(path())
    object[[i]]$sample = input$sample_text 
    object[[i]]$group = input$group_text 
    object[[i]][["percent.mt"]] = PercentageFeatureSet(object = object[[i]], pattern = "^mt-")
    
    output$nfeat_plot = renderPlot(
      FeatureScatter(object = object[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    )
    output$percent_plot = renderPlot( 
      FeatureScatter(object = object[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
    )
    updateCheckboxInput(session, "check_addmtfeat", value =1)
    shinyjs::show("add_mtandnfeat")
    object[[i]] <<- object[[i]]
  })
  
  
  observeEvent(input$add_mtandnfeat,{
    shinyjs::disable("for_featmin")
    shinyjs::disable("for_featmax")
    shinyjs::disable("for_mt")
    shinyjs::hide("add_mtandnfeat")
    l1 = input$for_featmin
    l2 = input$for_featmax
    l3 = input$for_mt
    
    if(input$check_addnumberofcells == 1){ #if checkbox clicked to sample genes 
      object[[i]] <- subset(object[[i]], cells = sample(Cells(object[[i]]), input$for_samp), subset = nFeature_RNA > l1 & nFeature_RNA < l2 & percent.mt < l3)
    }else{ #if not to sample genes 
      object[[i]] <- subset(object[[i]], subset = nFeature_RNA > l1 & nFeature_RNA < l2 & percent.mt < l3)
    }        
    
    object[[i]][["percent.mt"]] <- PercentageFeatureSet(object = object[[i]], pattern = "^mt-")
    object[[i]] <<- SCTransform(object[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
    
    assign(paste("all_genes_",i), rownames(object[[i]])) #makes a new variable with the gene names
    updateCheckboxInput(session, "check_addmtfeat", value =0)
    updateCheckboxInput(session, "check_addnumberofcells", value =0)
    shinyjs::enable("add_ano_button")
    shinyjs::enable("fin_button")
    updateTabsetPanel(session, "hidden_tabs_upload", selected = "panel1")
    
    
    output$files_so_far_text = renderText({
      ptext = c()
      for(t in 1:length(object)){
        ptext[t] = paste(t, ": ", unique(object[[t]]$sample, " |"))
      }
      ptext #to print the seurat files loaded so far
    })
    i <<- i+1
  })
  
  
  
  observeEvent(input$add_ano_button,{
    #reset all the values to add a new seurat object
    path(NULL)
    updateTextInput(session, "filechosen", value = "")
    updateTextInput(session, "sample_text", value = "")
    updateTextInput(session, "group_text", value = "")
    shinyjs::enable("sample_text")
    shinyjs::enable("group_text")
    shinyjs::enable("for_featmin")
    shinyjs::enable("for_featmax")
    shinyjs::enable("for_mt")
    shinyjs::enable("for_samp")
    shinyjs::enable("filechoose")
    

    updateNumericInput(session, "for_featmin", value =1)
    updateNumericInput(session, "for_featmax", value =1)
    updateNumericInput(session, "for_mt", value =1)
    updateCheckboxInput(session, "check_addmtfeat", value =0)
    updateCheckboxInput(session, "check_addnumberofcells", value =0)
    shinyjs::enable("upload_button")
    shinyjs::disable("add_ano_button")
    shinyjs::enable("fin_button")
    updateTabsetPanel(session, "hidden_tabs_upload", selected = "panel1")
    
    
  })
  
  observeEvent(input$fin_button,{
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_INTEGRATE")
  })
  
  observeEvent(input$but_int,{
    shinyjs::disable("but_int")
    TC.anchors = FindIntegrationAnchors(object.list = object, dims = 1:23)
    # list(for(t in 1:length(object)){object[[t]]})
    Combined <- IntegrateData(anchorset = TC.anchors, dims = 1:23)
    remove(TC.anchors)
    gc()
    Combined <<- Combined
    
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_DIM")
    updateTabsetPanel(session, "hidden_tabs_dim", selected = "pan_dim_1")
    
  })
  
  
  
  observeEvent(input$but_dimclust,{ #this tab is run first 
    #button to direct to cluster tab after dimension code
    updateTabsetPanel(session, "hidden_tabs_dim", selected = "pan_dim_2")
    shinyjs::disable("but_dimclust")
    
    if(input$check_adddim == 1){
      shinyjs::disable("for_dim")
      # remove(object)
      gc()
      DefaultAssay(Combined) <- "integrated"
      Combined <- ScaleData(Combined, verbose = FALSE)
      Combined <- RunPCA(Combined, npcs = 50, verbose = FALSE)
      
      dim=input$for_dim
      Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
      Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
      Combined <- FindClusters(Combined, resolution = 1)
      Combined <- FindClusters(Combined, resolution = 0.1)
      Combined[["UMI"]] <-  Combined$nCount_RNA  # Why divided by 100
      Combined[["genes"]] <-  Combined$nFeature_RNA
      
      Combined <<- Combined
      
      shinyjs::hide("but_dimclust")
      shinyjs::show("but_dimrep")
      
      updateTabsetPanel(session, "hidden_tabs", selected = "PAN_SUBCLUST")
      updateTabsetPanel(session, "hidden_tabs_subclust", selected = "pan_clus_1")
      updateCheckboxGroupInput(session, "for_clust", choices = levels(Combined$seurat_clusters))
      shinyjs::enable("for_dim")
      
      
      return()
    }
    
    DefaultAssay(Combined) <- "integrated"
    Combined <- ScaleData(Combined, verbose = FALSE)
    Combined <- RunPCA(Combined, npcs = 50, verbose = FALSE)
    
    output$elbow_plot = renderPlot(
      ElbowPlot(object = Combined, ndims = 50)
    )
    updateCheckboxInput(session, "check_adddim", value = 1)
    shinyjs::show("but_dclust")
    Combined <<- Combined
    
  }) 
  
  observeEvent(input$but_dclust,{
    shinyjs::disable("but_dclust")
    shinyjs::disable("for_dim")
    DefaultAssay(Combined) <- "integrated"
    Combined <- ScaleData(Combined, verbose = FALSE)
    Combined <- RunPCA(Combined, npcs = 50, verbose = FALSE)
    
    dim=input$for_dim
    Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
    Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
    Combined <- FindClusters(Combined, resolution = 0.1)
    Combined[["UMI"]] <-  Combined$nCount_RNA  # Why divided by 100
    Combined[["genes"]] <-  Combined$nFeature_RNA
    Combined <<- Combined
    
    shinyjs::hide("but_dimclust")
    shinyjs::show("but_dimrep")
    shinyjs::hide("but_dclust")
    
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_SUBCLUST")
    updateTabsetPanel(session, "hidden_tabs_subclust", selected = "pan_clus_1")
    updateCheckboxGroupInput(session, "for_clust", choices = levels(Combined$seurat_clusters))
    shinyjs::enable("for_dim")
    
    
  })
  
  
  
  

  observeEvent(input$but_dimrep,{ 
    #button to direct to repeat page after dimension code
    updateTabsetPanel(session, "hidden_tabs_dim", selected = "pan_dim_2")
    shinyjs::disable("but_dimrep")
    
    if(input$check_adddim == 1){
      
      # remove(object)
      gc()
      DefaultAssay(Combined) <- "integrated"
      Combined <- ScaleData(Combined, verbose = FALSE)
      Combined <- RunPCA(Combined, npcs = 50, verbose = FALSE)
      
      dim=input$for_dim
      Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
      Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
      Combined <- FindClusters(Combined, resolution = 0.1)
      Combined[["UMI"]] <-  Combined$nCount_RNA  # Why divided by 100
      Combined[["genes"]] <-  Combined$nFeature_RNA
      
      Combined <<- Combined
      
      #to show repeat option
      updateTabsetPanel(session, "hidden_tabs", selected = "PAN_SUBCLUST")
      updateTabsetPanel(session, "hidden_tabs_subclust", selected = "pan_clus_2")
      shinyjs::show("but_cont")
      shinyjs::show("but_rep")
      shinyjs::hide("but_sub")
      updateCheckboxInput(session, "check_clus", value =0)
      output$umi_plot = renderPlot({
        FeaturePlot(Combined, features = "UMI")    
      })
      output$gene_plot = renderPlot({
        FeaturePlot(Combined, features = "genes")    
      })
      output$dim_plot = renderPlot({
        DimPlot(Combined, reduction = "umap", label = TRUE)
      })
      return()
    }
    
    DefaultAssay(Combined) <- "integrated"
    Combined <- ScaleData(Combined, verbose = FALSE)
    Combined <- RunPCA(Combined, npcs = 50, verbose = FALSE)
    
    output$elbow_plot = renderPlot(
      ElbowPlot(object = Combined, ndims = 50)
    )
    updateCheckboxInput(session, "check_adddim", value = 1)
    shinyjs::show("but_drep")
    Combined <<- Combined
    
  }) 
  
  observeEvent(input$but_drep,{
    shinyjs::disable("but_drep")
    
    dim=input$for_dim
    Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
    Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
    Combined <- FindClusters(Combined, resolution = 0.1)
    Combined[["UMI"]] <-  Combined$nCount_RNA  # Why divided by 100
    Combined[["genes"]] <-  Combined$nFeature_RNA
    Combined <<- Combined
    
    
    #to show repeat option
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_SUBCLUST")
    updateTabsetPanel(session, "hidden_tabs_subclust", selected = "pan_clus_2")
    shinyjs::show("but_cont")
    shinyjs::show("but_rep")
    shinyjs::hide("but_sub")
    updateCheckboxInput(session, "check_clus", value =0)
    output$umi_plot = renderPlot({
      FeaturePlot(Combined, features = "UMI")    
    })
    output$gene_plot = renderPlot({
      FeaturePlot(Combined, features = "genes")    
    })
    output$dim_plot = renderPlot({
      DimPlot(Combined, reduction = "umap", label = TRUE)
    })
  })
  
  
  
  observeEvent(input$but_sub,{
    shinyjs::disable("but_sub")
    if(input$check_clus == 1){
      iden = input$for_clust
      
      Combined <- subset(Combined, idents = dput(as.character(iden)), invert = FALSE)
      DefaultAssay(object = Combined) <- "integrated"
      Combined <<- Combined
      
      
      updateTabsetPanel(session, "hidden_tabs", selected = "PAN_DIM")
      updateTabsetPanel(session, "hidden_tabs_dim", selected = "pan_dim_1")
      shinyjs::show("but_drep")
      shinyjs::enable("but_drep")
      
      shinyjs::show("but_dimrep")
      shinyjs::enable("but_dimrep")
      
      
      shinyjs::hide("but_dclust")
      shinyjs::hide("but_dimclust")
      
      updateNumericInput(session, "for_dim", value = 1)
      updateCheckboxInput(session, "check_adddim", value = 0)
      
      
      return()
    }
    
    updateTabsetPanel(session, "hidden_tabs_subclust", selected = "pan_clus_2")
    
    output$umi_plot = renderPlot({
      FeaturePlot(Combined, features = "UMI")    
    })
    output$gene_plot = renderPlot({
      FeaturePlot(Combined, features = "genes")    
    })
    output$dim_plot = renderPlot({
      DimPlot(Combined, reduction = "umap", label = TRUE)
    })
    updateCheckboxInput(session, "check_clus", value =1)
    shinyjs::show("but_clus")
  })
  
  
  observeEvent(input$but_clus,{
    shinyjs::disable("but_clus")
    iden = input$for_clust
    
    Combined <- subset(Combined, idents = dput(as.character(iden)), invert = FALSE)
    DefaultAssay(object = Combined) <- "integrated"
    Combined <<- Combined
    
    
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_DIM")
    updateTabsetPanel(session, "hidden_tabs_dim", selected = "pan_dim_1")
    shinyjs::show("but_drep")
    shinyjs::enable("but_drep")
    
    shinyjs::show("but_dimrep")
    shinyjs::enable("but_dimrep")
    
    
    shinyjs::hide("but_dclust")
    shinyjs::hide("but_dimclust")
    
    updateNumericInput(session, "for_dim", value = 1)
    updateCheckboxInput(session, "check_adddim", value = 0)
  })
  
  
  
  observeEvent(input$but_rep,{ 
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_SUBCLUST")
    updateTabsetPanel(session, "hidden_tabs_subclust", selected = "pan_clus_1")
    
    shinyjs::hide("but_dclust")
    shinyjs::enable("but_drep")
    shinyjs::hide("but_drep")
    shinyjs::enable("but_dimrep")
    shinyjs::hide("but_dimclust")
    
    updateNumericInput(session, "for_dim", value = 1)
    updateCheckboxInput(session, "check_adddim", value = 0)
    
    shinyjs::hide("but_clus")
    shinyjs::enable("but_clus")
    shinyjs::enable("but_sub")
    shinyjs::show("but_sub")
    
    shinyjs::hide("but_cont")
    shinyjs::hide("but_rep")
    updateCheckboxInput(session, "check_clus", value =0)
    
    updateCheckboxGroupInput(session, "for_clust", choices = levels(Combined$seurat_clusters))
  })
  
  observeEvent(input$but_cont,{
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_LOOPS")
    updateSelectizeInput(session = session, inputId = 'get_slingstart', choices = levels(Combined$seurat_clusters), server = TRUE)
    updateSelectizeInput(session = session, inputId = 'get_slingend', choices = levels(Combined$seurat_clusters), server = TRUE)
    updateSelectizeInput(session = session, inputId = 'get_mon2start', choices = levels(Combined$seurat_clusters), server = TRUE)
    updateSelectizeInput(session = session, inputId = 'get_mon3start', choices = levels(Combined$seurat_clusters), server = TRUE)
  })
  
  
  
  observeEvent(input$but_loop,{
    #calls the functions to run the raw code:
    updateTabsetPanel(session, "hidden_tabs", selected = "PAN_WHENRUNNING")
    
  })   
  
  observe({ 
    if(input$hidden_tabs == "PAN_WHENRUNNING"){ 
      #remove
      saveRDS(Combined, "Combined_newSeuratObjects.rds") #allows the user to save the Combined Seurat object 
      # Combined <- readRDS("/data/2623287c/Project1/tcelltest/Combined.filt.rds")

      #upload
      if(input$check_res == 1){
        print("single cell analysis")
        print(input$get_res[1]+0.05)
        tryCatch(
          upload_loops(Combined, as.numeric(input$get_res[1]+0.05), as.numeric(input$get_res[2]+0.05)), error = function(e){e} #add 0.5 because get_res 5 is only added to the label and not to the value in get_res
        )
      }else{
        tryCatch(
          upload_loops(Combined, 0.15, 0.55), error = function(e){e}
        )
      }

      #slingshot
      if(input$check_slingstart == 1){ #if user has selected start cluster(s)
        if(input$check_slingend == 1){ #if user has selected end cluster(s)
          print("slingshot")
          tryCatch(
            slingshot_loops(input$get_slingstart, input$get_slingend), error = function(e){e}
          )
        }else{
          tryCatch(
            slingshot_loops(input$get_slingstart, NULL), error = function(e){e}
          )
        }
      }else{
        if(input$check_slingend == 1){
          tryCatch(
            slingshot_loops(NULL, input$get_slingend), error = function(e){e}
          )
        }else{
          tryCatch(
            slingshot_loops(NULL, NULL), error = function(e){e}
          )
        }
      }

      #tradeSeq
      if(input$check_tradeknots == 1){
        print("tradeSeq")
        print(input$get_tradenknots)
        tryCatch(
          tradeseq_loops(input$get_tradenknots), error = function(e){e}
        )
      }else{
        tryCatch(
          tradeseq_loops(3), error = function(e){e}
        )
      }

      #monocle 2
      if(input$check_mon2start==1){
        print("monocle 2")
        tryCatch(
          monocle2_loops(input$get_mon2start), error = function(e){e}
        )
      }else{
        tryCatch(
          monocle2_loops(NULL), error = function(e){e}
        )
      }

      #monocle 3
      if(input$check_mon3start == 1){ #if user has selected start cluster(s)
        print("monocle 3")
        tryCatch(
          monocle3_loops(input$get_mon3start), error = function(e){e}
        )
      }else{
        tryCatch(
          monocle3_loops(NULL), error = function(e){e}
        )
      }

      Sys.time()
      
      updateTabsetPanel(session, "hidden_tabs", selected = "PAN_AFTERRUN")
      #reads in the files created as global variables available to the rest of the code 
      if(file.exists("tcells_combined_clusters_tables_res.rds")){tcells_combined_clusters_tables_res<<-readRDS("tcells_combined_clusters_tables_res.rds")}
      if(file.exists("tcells_combined_de_tables.rds")){tcells_combined_de_tables<<-readRDS("tcells_combined_de_tables.rds")}
      if(file.exists("tcells_combined_de_ggplots_table.rds")){tcells_combined_de_ggplots_table<<-readRDS("tcells_combined_de_ggplots_table.rds")}
      if(file.exists("all_genes_common_in_all_groups.rds")){all_genes_common_in_all_groups<<-readRDS("all_genes_common_in_all_groups.rds")}
      
      
      ##Reading in Pseudotime
      #UMAP
      if(file.exists("sds.rds")){sds <<- readRDS("sds.rds")}
      if(file.exists("sce.rds")){sce <<- readRDS("sce.rds")}
      # if(file.exists("counts.rds")){counts <<- readRDS("counts.rds")}
      # if(file.exists("clusters.rds")){clusters <<- readRDS("clusters.rds")}
      
      # PHATE
      if(file.exists("sdsPhate.rds")){sdsPhate <<- readRDS("sdsPhate.rds")}
      if(file.exists("scePhate.rds")){scePhate <<- readRDS("scePhate.rds")}
      if(file.exists("clustersPhate.rds")){clustersPhate <<- readRDS("clustersPhate.rds")}
      
      if(file.exists("cell_colors_clust.rds")){cell_colors_clust <<- readRDS("cell_colors_clust.rds")}
      
      if(file.exists("slingUMAPHeat_ALL.rds")){slingUMAPHeat_ALL <<- readRDS("slingUMAPHeat_ALL.rds")}
      if(file.exists("slingPHATEHeat_ALL.rds")){slingPHATEHeat_ALL <<- readRDS("slingPHATEHeat_ALL.rds")}
      if(file.exists("tradeUMAPHeat_ALL.rds")){tradeUMAPHeat_ALL <<- readRDS("tradeUMAPHeat_ALL.rds")}
      if(file.exists("tradePHATEHeat_ALL.rds")){tradePHATEHeat_ALL <<- readRDS("tradePHATEHeat_ALL.rds")}
      
      
      # MONOCLE 2
      
      if(file.exists("cds2.rds")){cds2 <<- readRDS("cds2.rds")}
      if(file.exists("new_monocle2_heatmap.rds")){new_monocle2_heatmap <<- readRDS("new_monocle2_heatmap.rds")}
      
      # MONOCLE3
      
      if(file.exists("cds3.rds")){cds3 <<- readRDS("cds3.rds")}
      if(file.exists("new_monocle3_heatmap.rds")){new_monocle3_heatmap <<- readRDS("new_monocle3_heatmap.rds")}
      
    }
  })
  
  counts = reactive({
    if(input$but_analysis == 1){
      if(file.exists("counts.rds")){readRDS("counts.rds")}
    }
  })
  
  clusters = reactive({
    if(input$but_analysis == 1){
      if(file.exists("clusters.rds")){readRDS("clusters.rds")}
    }
  })
  
  
  
  
  observeEvent(input$but_analysis,{
    
    output$men_cluster<- renderUI({
      
      modify_stop_propagation(menuItem("Cluster exploration", tabName = "cluster_res", icon = icon("puzzle-piece"),
                                       menuSubItem("UMAP", tabName = "all_cluster_res"),
                                       menuSubItem("Cluster markers", tabName = "grps_cluster_res"),
                                       startExpanded = T
      ))
      
    })
    
    output$men_de<- renderUI({
      modify_stop_propagation(menuItem("Differential expression (DE)", tabName = "ra_de", icon = icon("balance-scale"),
                                       modify_stop_propagation(menuItem("Gene view", tabName = "vis", icon = icon("braille"),
                                                                        menuSubItem("Single gene view", tabName = "ge_vis_gene"),
                                                                        menuSubItem("Multiple gene view", tabName = "m_ge_vis_gene"), 
                                                                        startExpanded = T
                                       )),
                                       menuItem("Cell population view", tabName = "ge_vis_cell", icon = icon("braille")), 
                                       startExpanded = T))
    })
    
    output$men_pseudotime<- renderUI({
      modify_stop_propagation(menuItem("Pseudotime", tabName = "Pseudotime", icon = icon("hourglass-start"),
                                       menuSubItem("Slingshot", tabName = "slingshot_tab"),
                                       menuSubItem("tradeSeq", tabName = "trade_tab"), 
                                       menuSubItem("Monocle 2", tabName = "mon2_tab"),
                                       menuSubItem("Monocle 3", tabName = "mon3_tab"),
                                       startExpanded = T
      ))
    })
    
    updateSelectizeInput(session = session, inputId = 'de_genes', choices = all_genes_common_in_all_groups, selected = all_genes_common_in_all_groups[1], server = TRUE)
    updateSelectizeInput(session = session, inputId = 'select_markers_dotplot', choices = all_genes_common_in_all_groups, selected = fav_genes, server = TRUE)
    updateSelectizeInput(session = session, inputId = 'gene_trade', choices = all_genes_common_in_all_groups, selected = all_genes_common_in_all_groups[1], server = TRUE)
    updateSliderInput(session = session, inputId = 'clusters_res', value = input$get_res[1], min = input$get_res[1], max = input$get_res[2])
    
    
    
    
  })
  
  
  
  
  
  ######################
  ### Labelling clusters under subtitle 2.3
  ######################
  ##Generating dynamic fields for labelling UMAP clusters and initializing the fields with placeholders
  
  ##Preparing and plotting UMAP cluster markers for annotating cell types
  
  umap_cluster_modified_rna = reactive({
    if(input$but_analysis == 1){
      umap_cluster_modified_ul = umap_clusters()
      DefaultAssay(umap_cluster_modified_ul) = "RNA"
      umap_cluster_modified_ul
    }else{
      NULL
    }
  })
  
  
  output$cluster_annot <- renderUI({
    
    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster_names)){
      do.call(flowLayout, 
              lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
                textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = cluster_names[x+1])
              })
      )
      
    } else {
      do.call(flowLayout,
              lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
                textInput(inputId = paste("labeller", x), label = strong(paste("Input cluster", x, "label")), value = x)
              })
      )
    }
  })
  
  
  ##Storing names decided on by the researchers based on optimal clustering to start off the differential expression visualization
  annotation = reactiveValues(annot = cluster_names)
  
  ##Observer to allow updating of cluster names dynamically as typed in
  observe({
    
    req(unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    })))
    annotation$annot = unlist(lapply(0:(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]))-1), function(x) {
      new_cluster_name = input[[paste("labeller",x)]]
    }))
  })
  
  ##Dynamic input for selecting celltypes (clusters) for diffential expression visualization
  output$cluster_ids <- renderUI({
    umap_names = annotation$annot
    #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)
    
    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      umap_names = annotation$annot
      selectInput(inputId = "select_cell_type", label = strong(paste("Select cell population to compare gene expression between",conditions[1],"and", conditions[2],":")),choices = umap_names, multiple = F)
      
    } else {
      selectInput(inputId = "select_cell_type", label = strong("Select cell population to compare gene expression across conditions:"),choices = unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
    }
    
  })
  
  ##Renaming clusters
  umap_cluster_modified_ren_reo = reactive({
    umap_names = annotation$annot
    umap_cluster_modified1 = umap_cluster_modified_rna()
    if(length(unique(umap_cluster_modified1$seurat_clusters)) == length(umap_names)){
      names(umap_names) <- levels(umap_cluster_modified1)
      umap_cluster_modified1 <- RenameIdents(umap_cluster_modified1, umap_names)
      
    } else {
      umap_cluster_modified1
    }
    
  })
  
  ##Plotting labelled umap
  labelled_umap_r = reactive({
    
    if(length(unique(umap_cluster_modified_ren_reo()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster.colours)){
      DimPlot(umap_cluster_modified_ren_reo(), order = T,  pt.size = 1, label = TRUE, label.size = 6)#, cols = cluster.colours)
      
    } else {
      DimPlot(umap_cluster_modified_ren_reo(), order = T, pt.size = 1, label = TRUE, label.size = 6)
    }
    
  })
  output$labelled_umap = renderPlot({
    labelled_umap_r()
  })
  
  
  output$dwnl_lumap <- downloadHandler(
    filename = function(){paste("labelled_umap",input$lumap_format,sep="")},
    content = function(file){
      ggsave(file,plot=labelled_umap_r(), width = input$lumap_width, height = input$lumap_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  ######################
  ###END Labelling clusters under subtitle 2.3
  ######################
  
  
  ######################
  ### Selection of seurat object and plotting UMAP plots under subtitle 2.1
  ######################
  
  ##Plotting UMAP plots for clustering
  #plot split by sample
  umap_p_split = reactive({
    withProgress(message = 'Plotting',
                 detail = 'Please wait...',
                 value = 0.8,
                 {
                   if (length(unique(umap_cluster_modified_ren_reo()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster.colours)) {
                     DimPlot(
                       umap_cluster_modified_ren_reo(),
                       reduction = "umap",
                       split.by = "sample",
                       pt.size = 1,
                       label = TRUE,
                       label.size = 6
                     )#, cols = cluster.colours)
                     
                   } else {
                     DimPlot(
                       umap_cluster_modified_ren_reo(),
                       reduction = "umap",
                       split.by = "sample",
                       pt.size = 1,
                       label = TRUE,
                       label.size = 6
                     )
                   }
                   
                 })
  })
  
  output$all_groups = renderPlot({
    umap_p_split()
  })
  output$dwnl_grps <- downloadHandler(
    filename = function() {
      paste("umap_split_by_samples", input$grps_format, sep = "")
    },
    content = function(file) {
      ggsave(
        file,
        plot = umap_p_split(),
        width = input$grps_width,
        height = input$grps_height,
        units = "cm",
        dpi = 300
      )
    },
    contentType = "image"
  )
  
  #plot split by group
  umap_p_splitgroup = reactive({
    withProgress(message = 'Plotting',
                 detail = 'Please wait...',
                 value = 0.8,
                 {
                   if (length(unique(umap_cluster_modified_ren_reo()[["seurat_clusters"]][["seurat_clusters"]])) == length(cluster.colours)) {
                     DimPlot(
                       umap_cluster_modified_ren_reo(),
                       reduction = "umap",
                       split.by = "group",
                       pt.size = 1,
                       label = TRUE,
                       label.size = 6
                     )#, cols = cluster.colours)
                     
                   } else {
                     DimPlot(
                       umap_cluster_modified_ren_reo(),
                       reduction = "umap",
                       split.by = "group",
                       pt.size = 1,
                       label = TRUE,
                       label.size = 6
                     )
                   }
                   
                 })
  })
  
  output$splitby_group = renderPlot({
    umap_p_splitgroup()
  })
  
  output$dwnl_splitby_groups <- downloadHandler(
    filename = function() {
      paste("umap_split_by_groups", input$grps_format, sep = "")
    },
    content = function(file) {
      ggsave(
        file,
        plot = umap_p_splitgroup(),
        width = input$grps_width,
        height = input$grps_height,
        units = "cm",
        dpi = 300
      )
    },
    contentType = "image"
  )
  ######################
  ### END Selection of seurat object and plotting UMAP plots under subtitle 2.1
  ######################
  
  ######################
  ### Cluster markers plots and tables under subtitle 2.2
  ######################
  ##Dynamic input field for selecting cluster to plot table of markers
  output$dyn_clusters <- renderUI({
    selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_clusters()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
  })
  
  ##Dynamic input for selecting celltypes (clusters) for diffential expression visualization
  output$dyn_clusters <- renderUI({
    umap_names = annotation$annot
    #FindMarkers(stim_markers(), ident.1 = paste(input$select_cell_type, "KO", sep = "_"), ident.2 = paste(input$select_cell_type, "WT", sep = "_"), verbose = FALSE)
    
    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      umap_names = annotation$annot
      selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = umap_names, multiple = F)
      
    } else {
      selectInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), choices = unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]]), multiple = F)
    }
    
  })
  
  ##Displaying table of cluster markers for annotating cell types
  cluster_markers = reactive({
    req(input$marker_genes_cluster)
    umap_names = annotation$annot
    
    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      withProgress(message = 'Tabulating', {
        setProgress(detail = 'Please wait...', value = 0.4)
        marker_tb = tcells_combined_clusters_tables_res[[(input$clusters_res * 10)]][[as.numeric(match(input$marker_genes_cluster,umap_names))]] %>%
          rownames_to_column(var = 'gene') %>% dplyr::left_join(x = ., y = uniprot_info, by = c("gene" = colnames(uniprot_info)[1])) %>%
          dplyr::distinct(., .keep_all = T) %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>%
          dplyr::select(contains(c("avg", "adj")), colnames(uniprot_info)[2]) %>%
          mutate_if(is.numeric, ~sprintf("%.3f", .))        
        
        setProgress(detail = 'Please wait...', value = 0.8)
        return(marker_tb)
      })
    } else {
      # ra_macrophage_combined_de_tables_full[[1]][[(as.numeric(input$select_cell_type) + 1)]]
      withProgress(message = 'Tabulating', {
        setProgress(detail = 'Please wait...', value = 0.4)
        marker_tb = tcells_combined_clusters_tables_res[[(input$clusters_res * 10)]][[(as.numeric(input$marker_genes_cluster) + 1)]]%>%
          rownames_to_column(var = 'gene') %>% dplyr::left_join(x = ., y = uniprot_info, by = c("gene" = colnames(uniprot_info)[1])) %>%
          dplyr::distinct(., .keep_all = T) %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>%
          dplyr::select(contains(c("avg", "adj")), colnames(uniprot_info)[2]) %>%
          mutate_if(is.numeric, ~sprintf("%.3f", .))
        setProgress(detail = 'Please wait...', value = 0.8)
        return(marker_tb)
      })
    }
    
  })
  
  
  output$top_conserved_genes = DT::renderDataTable({
    # req(input$topclgenes_i)
    #numeric_cols =  colnames(data.frame(cluster_markers()))[which_numeric_cols(data.frame(cluster_markers()))]
    
    #   # Javascript-enabled table.
    #   datatable(
    DT::datatable(
      cluster_markers(),
      escape = F,
      selection = "single",
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller"),
      options = list(
        deferRender = TRUE,
        scrollY = 350,
        #scroller = TRUE,
        #lengthMenu = FALSE,
        autoWidth = FALSE,
        dom = "Blfrtip",
        buttons =
          list(list(
            extend = "collection",
            buttons = c("csv", "pdf"),
            text = "Download"
          )),  # end of buttons customization
        
        # customize the length menu
        lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
        pageLength = 10
      ), fillContainer = TRUE
    )
  }, server = TRUE)
  
  
  output$top_markers_umap <- renderUI({
    req(cluster_markers())
    selectInput(inputId = "select_markers_umap", label = strong("Select marker to visualize in clusters:"), choices = cluster_markers()[,1], multiple = T, selected = head(cluster_markers()[,1], n=4))
    
  })
  conserved_markers_umap_r = reactive({
    req(input$select_markers_umap)
    req(umap_cluster_modified_ren_reo())
    withProgress(message = 'Plotting',{
      shiny::setProgress( detail = 'Please wait...', value = 0.4)
      fp_conserved = FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
      shiny::setProgress( detail = 'Please wait...', value = 0.8)
      return(fp_conserved)
    })
  })
  
  output$conserved_markers_umap = renderPlot({
    req(input$select_markers_umap)
    req(umap_cluster_modified_ren_reo())
    conserved_markers_umap_r()
  })
  
  
  output$dwnl_markers <- downloadHandler(
    filename = function(){paste(input$select_markers_umap, "_feature_plot",input$markers_format,sep="")},
    content = function(file){
      ggsave(file,plot=conserved_markers_umap_r(), width = input$markers_width, height = input$markers_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  ##information box
  output$box_2_2 =renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     paste("Listing top cluster marker which can subsequently be used in labelling the cluster.Here, markers are genes highly expressed in a cluster as compared to all other clusters in both", conditions[1], "and", conditions[2],".")
              ),
              tags$hr()
    )
  })
  
  ######################
  ### END Cluster markers plots and tables under subtitle 2.2
  ######################
  
  
  ######################
  ### Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
  ######################
  stim_markers = reactive({
    
    umap_cluster_modified = umap_cluster_modified_ren_reo()
    umap_cluster_modified$celltype.group <- interaction(Idents(umap_cluster_modified), umap_cluster_modified$group, sep = "_")
    umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
    Idents(umap_cluster_modified) <- "celltype.group"
    umap_cluster_modified
  })
  
  
  #Functions to update differentially expressed genes
  de_stim_vs_ctrl_um_r = eventReactive(input$de_genes,{
    withProgress(message = 'Plotting',{
      shiny::setProgress(detail = 'Please wait...', value = 0.5)
      fp_umap = FeaturePlot(stim_markers(), features = input$de_genes, split.by = "group", max.cutoff = 3,cols = c("grey", "red"))
      shiny::setProgress(detail = 'Please wait...', value = 0.8)
      return(fp_umap)
    })
  })
  
  output$de_stim_vs_ctrl_um = renderPlot({
    req(input$de_genes)
    req(de_stim_vs_ctrl_um_r())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      
      de_stim_vs_ctrl_um_r()
    })
  })
  
  output$dwnl_featurep <- downloadHandler(
    filename = function(){paste(input$de_genes, "_feature_plot",input$featurep_format,sep="")},
    content = function(file){
      ggsave(file,plot=de_stim_vs_ctrl_um_r(), width = input$featurep_width, height = input$featurep_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  
  de_stim_vs_ctrl_vp_r = reactive({
    # plots <- VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE, multi.group = T, cols = group.cols, assay = "RNA")
    # for(i in 1:length(plots)) {
    #     plots[[i]] <- plots[[i]] + stat_summary(fun.y= median, geom='point', size = 2, colour = "black", position = position_dodge(0.9)) + scale_fill_manual(values=group.cols)
    # }
    # CombinePlots(plots)
    VlnPlot(stim_markers(), features = input$de_genes, split.by = "group", group.by = "celltype", pt.size = 0, combine = FALSE, cols = group.cols, assay = "RNA")
    
  })
  
  output$de_stim_vs_ctrl_vp = renderPlot({
    req(input$de_genes)
    req(stim_markers())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.7, {
      
      de_stim_vs_ctrl_vp_r()
    })
  })
  
  output$dwnl_violp <- downloadHandler(
    filename = function(){paste(input$de_genes,"_violin_plot",input$violp_format,sep="")},
    content = function(file){
      ggsave(file,plot=de_stim_vs_ctrl_vp_r(), width = input$violp_width, height = input$violp_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  ##Information box
  output$box_1_1 <- renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     
                     paste("Comparison of", input$de_genes, "expression between", cond," across all clusters using violin plots and umap feature plots.")
                     
              ),
              tags$hr()
    )  })
  
  ######################
  ### END Differential expression dynamic using UMAP and Violin plots under subtitle 1.1
  ######################
  
  
  ######################
  ### Differential expression using dotplot under subtitle 1.2
  ######################
  ##Dotplot for DE comparison between KO and WT across cell types
  marker_dotplot_r = eventReactive(input$select_markers_dotplot,{
    req(input$select_markers_dotplot)
    req(umap_cluster_modified_ren_reo())
    withProgress(message = 'Plotting',{
      shiny::setProgress(detail = 'Please wait...', value = 0.3)
      
      umap_cluster_modified_ren_reo = umap_cluster_modified_ren_reo()
      umap_cluster_modified_ren_reo@meta.data$grp_od <- umap_cluster_modified_ren_reo@meta.data$group
      umap_cluster_modified_ren_reo@meta.data <- umap_cluster_modified_ren_reo@meta.data %>% mutate(grp_od = case_when(grp_od == "Healthy" ~ 4,grp_od == "UPA" ~ 3,grp_od == "Naive RA" ~ 2,grp_od == "Resistant RA" ~ 1,grp_od == "Remission RA" ~ 0))
      shiny::setProgress(detail = 'Please wait...', value = 0.5)
      
      umap_cluster_modified_ren_reo@meta.data <- dplyr::arrange(umap_cluster_modified_ren_reo@meta.data, umap_cluster_modified_ren_reo@meta.data$grp_od)
      
      ##or **(not the '-' sign)**
      umap_cluster_modified_ren_reo@meta.data <- dplyr::arrange(umap_cluster_modified_ren_reo@meta.data, -umap_cluster_modified_ren_reo@meta.data$grp_od)
      dp = DotPlot(umap_cluster_modified_ren_reo, features = input$select_markers_dotplot, cols = group.cols, dot.scale = 6, split.by = "group") + RotatedAxis()
      shiny::setProgress(detail = 'Please wait...', value = 0.7)
      return(dp)
      
    })
  })
  
  output$marker_dotplot = renderPlot({
    marker_dotplot_r()
  })
  
  
  output$dwnl_dotp <- downloadHandler(
    filename = function(){paste(input$select_markers_dotplot,"dotplot",input$dotp_format,sep="_")},
    content = function(file){
      ggsave(file,plot=marker_dotplot_r(), width = input$dotp_width, height = input$dotp_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  ##Information box
  output$box_1_2 <- renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     
                     paste("Comparison of gene expression between", conditions[1], "and", conditions[2], "cells across clusters using a dotplot. The genes are on the y-axis and the clusters on the x-axis. Red is for", conditions[1], "and Blue for", conditions[2], "cells with the increase in intensity of the respective colour (from grey to blue/red) correlating with the increase in the average level of gene expression across all cells in the cluster. The size of the dot corresponds to the percentage of cells in the cluster expressing the gene.")
                     
              ),
              tags$hr()
    )})
  ######################
  ###END Differential expression using dotplot under subtitle 1.2
  ######################
  
  
  ######################
  ### Differential expression using ggplot and tables under subtitle 1.3
  ######################
  
  ##Retrieving table for DE expression from precomputed list
  genes_in_de_order = reactive({
    umap_names = annotation$annot
    
    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      tcells_combined_de_tables[[(input$clusters_res * 10)]][[as.numeric(match(input$select_cell_type,umap_names))]][[as.numeric(match(input$ra_conds,conds))]]
      
    } else {
      # ra_macrophage_combined_de_tables_full[[1]][[(as.numeric(input$select_cell_type) + 1)]]
      tcells_combined_de_tables[[(input$clusters_res * 10)]][[(as.numeric(input$select_cell_type) + 1)]][[as.numeric(match(input$ra_conds,conds))]]
      
    }
    
  })
  
  ##Retrieving table for DE expression from precomputed list
  top_de_g = reactive({
    req(input$select_cell_type, input$ra_conds)
    withProgress(message = 'Tabulating',{
      setProgress(detail = 'Please wait...', value = 0.4)
      t_d_g = genes_in_de_order() %>% rownames_to_column(var = 'gene') %>% inner_join(x=., y = uniprot_info, by = c("gene" = colnames(uniprot_info)[1])) %>%
        dplyr::distinct(., gene, .keep_all = T) %>% filter(p_val_adj <= 0.05) %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)) %>% mutate_if(is.numeric, ~sprintf("%.3f", .)) %>%
        select(gene, p_val, avg_logFC, pct.1, pct.2, p_val_adj, colnames(uniprot_info)[2])
      setProgress(detail = 'Please wait...', value = 0.8)
      return(t_d_g)
    })
  })
  
  output$top_de_genes = DT::renderDataTable({
    DT::datatable(
      top_de_g(),
      selection = "single",
      escape = F,
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller"),
      options = list(
        deferRender = TRUE,
        scrollY = 350,
        #scroller = TRUE,
        #lengthMenu = FALSE,
        autoWidth = FALSE,
        dom = "Blfrtip",
        buttons =
          list(list(
            extend = "collection",
            buttons = c("csv", "pdf"),
            text = "Download"
          )),  # end of buttons customization
        
        # customize the length menu
        lengthMenu = list( c(10, 20, -1), c(10, 20, "All")), # end of lengthMenu customization
        pageLength = 10
      )
    )
    # DT::formatSignif(columns = numeric_cols, digits = 3)
  }, server = TRUE)
  
  
  ##Allowing for download of DE table
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("differentially_expressed_genes_in",input$select_cell_type,input$ra_conds, ".csv", sep = "_")
      
    },
    content = function(file) {
      write.csv(genes_in_de_order() %>% rownames_to_column(var = 'genes') %>% filter(p_val_adj <= 0.05)  %>% mutate_at(vars(matches("p_val|pval") ), ~formatC(., format = "e", digits = 2)), file) %>% mutate_if(is.numeric, ~sprintf("%.3f", .))
    }
  )
  
  ##Retrieving table for DE scatterplotfrom precomputed list
  cell_type_de = reactive({
    req(input$select_cell_type, input$ra_conds)
    umap_names = annotation$annot
    
    if(length(unique(umap_cluster_modified_rna()[["seurat_clusters"]][["seurat_clusters"]])) == length(umap_names)){
      tcells_combined_de_ggplots_table[[(input$clusters_res * 10)]][[as.numeric(match(input$select_cell_type,umap_names))]]
      
    } else {
      tcells_combined_de_ggplots_table[[(input$clusters_res * 10)]][[(as.numeric(input$select_cell_type) + 1)]]
    }
    
  })
  
  
  ##preparing ggplot for average DE expression for genes above
  cell_type_de_plot_no_grb = reactive({
    req(input$select_cell_type, input$ra_conds)
    theme_set(theme_cowplot())
    ggplot(data=cell_type_de(), aes_string(paste("`",unlist(str_split(input$ra_conds[1], " VS "))[1],"`", sep=""), paste("`",unlist(str_split(input$ra_conds[1], " VS "))[2],"`", sep=""))) + geom_point() + ggtitle(input$select_cell_type) + theme_bw()
    
    
  })
  
  cell_type_de_plot = reactive({
    req(input$select_cell_type, input$ra_conds)
    #theme_set(theme_cowplot())
    grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
                              gp=gpar(col="red", fontsize=13, fontface="italic")))
    cell_type_de_plot_no_grb() + annotation_custom(grob)
    
    
  })
  
  ##plotting ggplot for average DE expression for genes above
  output$cell_type_plot = renderPlot({
    #theme_set(theme_cowplot())
    cell_type_de_plot()
  })
  
  output$dwnl_scatter <- downloadHandler(
    filename = function(){paste("scatter_plot_of_average_expression_in_", input$ra_conds,"_among_", input$select_cell_type, input$scatter_format,sep="")},
    content = function(file){
      ggsave(file,plot=cell_type_de_plot_no_grb(), width = input$scatter_width, height = input$scatter_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  ##Displaying further details upon clicking points
  displayed_text <- reactive({
    req(input$plot_click)
    nearPoints(cell_type_de(), input$plot_click)
    
  })
  
  ##Displaying table with gene details upon click of point in DE scatterplot
  output$click_info <- renderDataTable({
    req(displayed_text())
    displayed_text()
    # DT::datatable(displayed_text(),
    #               extensions=c('Scroller'),
    #               options = list(dom = 'Bfrtip',
    #                              scroller = TRUE,
    #                              scrollX=TRUE))
  }, escape = F)
  
  ##Information box
  output$box_1_3a <- renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     paste("Comparison of average gene expression between", input$ra_conds,"in", input$select_cell_type, "cells using a scatter plot.")
              ),
              tags$hr()
    )
  })
  
  output$box_1_3b <- renderUI({
    
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     paste("Listing differentially expressed genes (adjusted P value <0.05) between" , input$ra_conds, "in", input$select_cell_type,"cells.")
              ),
              tags$hr()
    )
  })
  ######################
  ### END Differential expression using ggplot and tables under subtitle 1.3
  ######################
  
  
  ##Plot Pseudotime####
  
  ##Plot Pseudotime####
  
  ##slingshot####
  
  #UMAP
  sling_make_legend_plot = reactive({
    x = plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    x = legend("topleft",title = "Clusters", legend = levels(unique(umap_cluster_modified_rna()$seurat_clusters)), 
               col = hue_pal()(length(unique(umap_cluster_modified_rna()$seurat_clusters))), ncol = 6, pch = 16)
  })
  
  sling_make_UMAP = reactive({
    # req(input$UMAP_phate)
    sling_toplot = plot(reducedDim(sds[[input$clusters_res * 10]]), col = cell_colors_clust[[input$clusters_res * 10]], lwd = 1) #UMAP
    sling_toplot = lines(sds[[input$clusters_res * 10]], lwd = 2, type = 'lineages', col = 'black')
    return(sling_toplot)
  })
  
  
  output$sling_UMAP_plot = renderPlot({
    # req(sling_make_UMAP())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      sling_make_UMAP()
    })
  })
  
  output$sling_UMAP_leg = renderPlot({
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      sling_make_legend_plot()
    })
  })
  
  
  output$dwnl_slingumap <- downloadHandler(
    filename = function(){paste("sling_umap",input$slingumap_format,sep="")},
    content = function(file){
      ggsave(file,plot=cell_type_de_plot_no_grb(), width = input$slingumap_width, height = input$slingumap_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  #PHATE
  phate_make_legend_plot = reactive({
    x = plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    x = legend("topleft",title = "Clusters", legend = levels(unique(umap_cluster_modified_rna()$seurat_clusters)), 
               col = hue_pal()(length(unique(umap_cluster_modified_rna()$seurat_clusters))), ncol = 6, pch = 16)  
  })
  
  sling_make_PHATE = reactive({
    # req(input$UMAP_phate)
    sling_toplot = plot(reducedDim(sdsPhate[[input$clusters_res * 10]]), col = cell_colors_clust[[input$clusters_res * 10]], lwd = 1) #UMAP
    sling_toplot = lines(sdsPhate[[input$clusters_res * 10]], lwd = 2, col = 'black')
    return(sling_toplot)
  })
  
  
  output$sling_PHATE_plot = renderPlot({
    # req(sling_make_UMAP())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      sling_make_PHATE()
    })
  })
  
  output$sling_PHATE_leg = renderPlot({
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      phate_make_legend_plot()
    })
  })
  
  
  output$dwnl_slingphate<- downloadHandler(
    filename = function(){paste("slingphate_plot",input$slingphate_format,sep="")},
    content = function(file){
      ggsave(file,plot=sling_make_PHATE(), width = input$slingphate_width, height = input$slingphate_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  
  #HEATMAP 
  sling_make_HEAT = reactive({
    # req(input$UMAP_phate)
    if(input$select_sling_reduc == "UMAP"){
      return(slingUMAPHeat_ALL[[input$clusters_res * 10]])    
    }
    return(slingPHATEHeat_ALL[[input$clusters_res * 10]])
    
  })
  
  output$sling_HEAT_PLOT = renderPlot({
    # req(sling_make_UMAP())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      sling_make_HEAT()
    })
  })
  
  
  sling_for_heat_table = reactive(
    if(input$select_sling_reduc == "UMAP"){
      data.frame(slingUMAPHeat_ALL[[input$clusters_res * 10]]@row_names_param[["labels"]][row_order(slingUMAPHeat_ALL[[input$clusters_res * 10]])]) %>% rename("Genes" = colnames(.))  %>% 
        dplyr::left_join(x = ., y = uniprot_info, by = c("Genes" = "Gene names  (primary )")) %>% dplyr::distinct(., Genes, .keep_all = T)  %>% dplyr::select(Genes, uniprot)
    }else{
      data.frame(slingPHATEHeat_ALL[[input$clusters_res * 10]]@row_names_param[["labels"]][row_order(slingPHATEHeat_ALL[[input$clusters_res * 10]])]) %>% rename("Genes" = colnames(.))  %>% 
        dplyr::left_join(x = ., y = uniprot_info, by = c("Genes" = "Gene names  (primary )")) %>% dplyr::distinct(., Genes, .keep_all = T)  %>% dplyr::select(Genes, uniprot)
      
    }
  )
  
  
  output$sling_heat_info <- renderDataTable({
    datatable(
      sling_for_heat_table(),
      selection = "single",
      escape = F,
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller")
      
    )
  }, server = TRUE) 
  
  output$dwnl_slingheat<- downloadHandler(
    filename = function(){paste("slingheat_plot",input$slingheat_format,sep="")},
    content = function(file){
      ggsave(file,plot=sling_make_HEAT(), width = input$slingheat_width, height = input$slingheat_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  
  
  #Slingshot information box
  output$box_3_1 =renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     paste("UMAP and PHATE plots produced by slingshot pseudotime analysis. Heatmap plot produced using variable features and calculation of q value. Heatmap ordered by q value")
              ),
              tags$hr()
    )
  })
  
  
  
  ##tradeseq####
  #UMAP
  trade_make_UMAP = reactive({
    # req(input$UMAP_phate)
    trade_toplot = plotGeneCount(curve = sds[[input$clusters_res * 10]],counts = counts()[[input$clusters_res * 10]],clusters = clusters()[[input$clusters_res * 10]],models = sce[[input$clusters_res * 10]])
    return(trade_toplot)
  })
  
  output$trade_UMAP_plot = renderPlot({
    # req(sling_make_plot())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      trade_make_UMAP()
    })
  })
  
  output$dwnl_tradeumap<- downloadHandler(
    filename = function(){paste("tradeumap_plot",input$tradeumap_format,sep="")},
    content = function(file){
      ggsave(file,plot=trade_make_UMAP(), width = input$tradeumap_width, height = input$tradeumap_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  
  #PHATE
  trade_make_PHATE = reactive({
    # req(input$UMAP_phate)
    trade_toplot = plotGeneCount(curve = sdsPhate[[input$clusters_res * 10]],counts = counts()[[input$clusters_res * 10]],clusters = clustersPhate[[input$clusters_res * 10]],models = scePhate[[input$clusters_res * 10]])
    return(trade_toplot)
  })
  
  output$trade_PHATE_plot = renderPlot({
    # req(sling_make_plot())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      trade_make_PHATE()
    })
  })
  
  output$dwnl_tradephate<- downloadHandler(
    filename = function(){paste("tradephate_plot",input$tradephate_format,sep="")},
    content = function(file){
      ggsave(file,plot=trade_make_PHATE(), width = input$tradephate_width, height = input$tradephate_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  
  #GENE PLOT
  observeEvent(input$select_gene_reduc,{
    if(input$but_analysis == 1){
      if(input$select_gene_reduc == "UMAP"){
        updateSelectizeInput(session = session, inputId = 'gene_trade', 
                             choices =tradeUMAPHeat_ALL[[input$clusters_res * 10]]@row_names_param[["labels"]], 
                             selected = choice_gene, server = TRUE)
        
      }else{
        updateSelectizeInput(session = session, inputId = 'gene_trade', 
                             choices = tradePHATEHeat_ALL[[input$clusters_res * 10]]@row_names_param[["labels"]], 
                             selected = choice_gene, server = TRUE)
        
      }
    }
  })
  
  
  
  trade_make_GENE = reactive({
    if(input$select_gene_reduc == "UMAP"){
      trade_toplot = plotGeneCount(sds[[input$clusters_res * 10]], counts()[[input$clusters_res * 10]], gene = input$gene_trade)
      return(trade_toplot)
    }
    trade_toplot = plotGeneCount(sdsPhate[[input$clusters_res * 10]], counts()[[input$clusters_res * 10]], gene = input$gene_trade)
    return(trade_toplot)
  })
  
  
  
  output$trade_GENE_plot = renderPlot({
    # req(sling_make_plot())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      trade_make_GENE()
    })
  })
  
  
  #SMOOTHERS PLOT
  trade_make_SMOOTH = reactive({
    if(input$select_gene_reduc == "UMAP"){
      return(plotSmoothers(sce[[input$clusters_res * 10]], counts()[[input$clusters_res * 10]], gene = input$gene_trade, lwd = 2))
    }
    return(plotSmoothers(scePhate[[input$clusters_res * 10]], counts()[[input$clusters_res * 10]], gene = input$gene_trade, lwd = 2))
  })
  
  
  output$trade_SMOOTH_plot = renderPlot({
    # req(sling_make_plot())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      trade_make_SMOOTH()
    })
  })
  
  
  #HEATMAP 
  trade_make_HEAT = reactive({
    # req(input$UMAP_phate)
    if(input$select_trade_reduc == "UMAP"){
      return(tradeUMAPHeat_ALL[[input$clusters_res * 10]]) 
    }
    return(tradePHATEHeat_ALL[[input$clusters_res * 10]]) 
  })
  
  output$trade_HEAT_PLOT = renderPlot({
    # req(trade_make_UMAP())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      trade_make_HEAT()
    })
  })
  
  
  trade_for_heat_table = reactive(
    if(input$select_trade_reduc == "UMAP"){
      data.frame(tradeUMAPHeat_ALL[[input$clusters_res * 10]]@row_names_param[["labels"]][row_order(tradeUMAPHeat_ALL[[input$clusters_res * 10]])]) %>% rename("Genes" = colnames(.))  %>%
        dplyr::left_join(x = ., y = uniprot_info, by = c("Genes" = "Gene names  (primary )")) %>% dplyr::distinct(., Genes, .keep_all = T)  %>% dplyr::select(Genes, uniprot) 
    }else{
      data.frame(tradePHATEHeat_ALL[[input$clusters_res * 10]]@row_names_param[["labels"]][row_order(tradePHATEHeat_ALL[[input$clusters_res * 10]])]) %>% rename("Genes" = colnames(.))  %>%
        dplyr::left_join(x = ., y = uniprot_info, by = c("Genes" = "Gene names  (primary )")) %>% dplyr::distinct(., Genes, .keep_all = T)  %>% dplyr::select(Genes, uniprot) 
    }
  )
  
  output$trade_heat_info <- renderDataTable({
    datatable(
      trade_for_heat_table(),
      selection = "single",
      escape = F,
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller")
      
    )
  }, server = TRUE) 
  
  #tradeSeq information box
  output$box_3_2 =renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     paste("UMAP and PHATE plots produced by tradeSeq pseudotime analysis. Heatmap plot produced using tradeSeq association test. tradeSeq is used downstream from slingshot")
              ),
              tags$hr()
    )
  })
  
  
  output$dwnl_tradeheat<- downloadHandler(
    filename = function(){paste("tradeheat_plot",input$tradeheat_format,sep="")},
    content = function(file){
      ggsave(file,plot=trade_make_HEAT(), width = input$tradeheat_width, height = input$tradeheat_height, units = "cm", dpi = 300)
    },
    contentType = "image"
  )
  
  
  
  #monocle3
  #UMAP
  mon_make_UMAP = reactive({
    # req(input$UMAP_phate)
    mon_toplot =  plot_cells(cds3[[input$clusters_res * 10]],
                             color_cells_by = input$colour_by_mon,
                             label_cell_groups=FALSE,
                             label_leaves=TRUE,
                             label_branch_points=TRUE,
                             graph_label_size=1.5)
    
    return(mon_toplot)
  })
  
  output$mon_UMAP_plot = renderPlot({
    # req(trade_make_plot())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      mon_make_UMAP()
    })
  })
  
  
  #HEATMAP
  mon3_make_HEAT = reactive({
    # req(input$UMAP_phate)
    return(new_monocle3_heatmap[[input$clusters_res * 10]])
    
  })
  
  output$mon3_HEAT_plot = renderPlot({
    # req(sling_make_plot())
    mon3_make_HEAT()
  })
  
  mon_3_for_heat_table = reactive(
    data.frame(new_monocle3_heatmap[[input$clusters_res * 10]]@row_names_param$labels) %>% rename("Genes" = colnames(.))  %>% 
      dplyr::left_join(x = ., y = uniprot_info, by = c("Genes" = "Gene names  (primary )")) %>% 
      dplyr::distinct(., Genes, .keep_all = T)  %>% dplyr::select(Genes, uniprot)
    
  )
  
  output$mon3_heat_info <- renderDataTable({
    datatable(
      mon_3_for_heat_table(),
      selection = "single",
      escape = F,
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller")
      
    )
  }, server = TRUE) 
  
  #monocle 3 information box
  output$box_3_4 =renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     paste("Monocle 3 UMAP plot and heatmap")
              ),
              tags$hr()
    )
  })
  
  
  #monocle2
  #TRAJECTORY
  mon_make_TRA = reactive({
    # req(input$UMAP_phate)
    mon_toplot =  plot_cell_trajectory(cds2[[input$clusters_res * 10]], color_by = input$colour_by_mon2)
    return(mon_toplot)
  })
  
  output$mon2_TRA_plot = renderPlot({
    # req(sling_make_plot())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      mon_make_TRA()
    })
  })
  
  #HEATMAP
  mon_make_HEAT = reactive({
    # req(input$UMAP_phate)
    return(new_monocle2_heatmap[[input$clusters_res * 10]])
  })
  
  output$mon2_HEAT_plot = renderPlot({
    # req(sling_make_plot())
    withProgress(message = 'Plotting', detail = 'Please wait...', value = 0.9, {
      mon_make_HEAT()
    })
  })
  
  mon_2_for_heat_table = reactive(
    data.frame(new_monocle2_heatmap[[input$clusters_res * 10]]@row_names_param[["labels"]]) %>% rename("Genes" = colnames(.))  %>%
      dplyr::left_join(x = ., y = uniprot_info, by = c("Genes" = "Gene names  (primary )")) %>% dplyr::distinct(., Genes, .keep_all = T)  %>% dplyr::select(Genes, uniprot)
  )
  
  output$mon2_heat_info <- renderDataTable({
    datatable(
      mon_2_for_heat_table(),
      selection = "single",
      escape = F,
      rownames = FALSE,
      filter = list(position = "top", plain = TRUE),
      style = "default",
      extensions = c("Buttons","Scroller")
      
    )
  }, server = TRUE) 
  
  #monocle 2 information box
  output$box_3_3 =renderUI({
    wellPanel(style = "background:#385A4F",
              tags$hr(),
              tags$p(style = "font-family:Arial;color:white",
                     paste("Monocle 2 trajectory plot and heatmap")
              ),
              tags$hr()
    )
  })
  
}


shinyApp(ui = ui, server = server)

