

library("shiny")
library("shinyjqui")
library("GUniFrac")
library("vegan")
library("ade4")
library("phangorn")
library("cluster")
library("fpc")
library("compare")
library("plotrix")
library("PerformanceAnalytics")
library("reshape")
library("ggplot2")
library("gridExtra")
library("grid")
library("ggrepel")
library("gtable")
library("Matrix")
library("cowplot")
library("Hmisc")
library("corrplot")
library("shinycssloaders")
library("fs")
##### Load and install shiny libraries alternative (not on shinyio) #####
# Check if required packages are already installed, and install if missing
packages <-
  c("shinyjqui",
    "GUniFrac",
    "vegan",
    "ade4",
    "phangorn",
    "cluster",
    "fpc",
   "compare",
"plotrix",
"PerformanceAnalytics",
"reshape",
"ggplot2",
"gridExtra",
"grid","ggrepel","gtable","Matrix","cowplot","Hmisc","corrplot","shinycssloaders")

# # Function to check whether the package is installed
InsPack <- function(pack)
{
  if ((pack %in% installed.packages()) == FALSE) {
    install.packages(pack, repos = "http://cloud.r-project.org/")
  }
}

# # Applying the installation on the list of packages
lapply(packages, InsPack)

# # Make the libraries
lib <- lapply(packages, require, character.only = TRUE)
#
# # Check if it was possible to install all required libraries
# flag <- all(as.logical(lib))
#####
flag <- T
#####
ui <- function(req) {
  #### Page layout
  fluidPage(
    #### Page Head
    tags$head(
      #### Javascript script that disables tabs + custom javascript function to enable them again when getting signal
      tags$script(
        "
        window.onload = function() {
        $('#mynavlist a:contains(\"Alpha-Diversity\")').attr('data-toggle','none').parent().addClass('disabled');
        $('#mynavlist a:contains(\"Beta-Diversity\")').attr('data-toggle','none').parent().addClass('disabled');
        $('#mynavlist a:contains(\"Taxonomic-Binning\")').attr('data-toggle','none').parent().addClass('disabled');
        $('#mynavlist a:contains(\"Create-Input-Tables\")').attr('data-toggle','none').parent().addClass('disabled');
        $('#mynavlist a:contains(\"Serial-Group-Comparisons\")').attr('data-toggle','none').parent().addClass('disabled');
        $('#mynavlist a:contains(\"Correlations\")').attr('data-toggle','none').parent().addClass('disabled');
        };
        
        Shiny.addCustomMessageHandler('activeNavs', function(nav_label) {
        $('#mynavlist a:contains(\"' + nav_label + '\")').attr('data-toggle','tab').parent().removeClass('disabled');
        });

        Shiny.addCustomMessageHandler('deactiveNavs', function(nav_label) {
        $('#mynavlist a:contains(\"' + nav_label + '\")').attr('data-toggle','none').parent().addClass('disabled');
        });
        "
      )
      ),
    ##### Main page
    navbarPage(
      title = "Rhea Pipeline",
      id = "mynavlist",
      # First tab
      tabPanel(title = "Normalization",
               # Define it is sidebar layout
               sidebarLayout(
                 # Sidebar (left side of screen) Used for inputs (Not obliged: help message)
                 sidebarPanel(
                   
                   fileInput(
                     inputId = "OTUtab",
                     label = "OTU table",
                     accept = ".tab"
                   ),
                   
                   radioButtons(
                     inputId = "method",
                     label = "Method",
                     choiceNames = c(
                       "No random subsampling, no rounding",
                       "Random subsampling with rounding"
                     ),
                     choiceValues = c(0, 1)
                   ),
                   
                   numericInput(
                     inputId = "labelCutoff",
                     label = "Label cutoff to present independently",
                     value = 5
                   ),
                   
                   actionButton(
                     inputId = "Submit1",
                     label = "Run normalization",
                     icon = icon("refresh")
                   ),
                   
                   tags$br(),
                   
                   tags$br(),
                   
                   actionButton(
                     inputId = "Help1",
                     label = "Help",
                     icon = icon("question-circle")
                   ),
                   
                   htmlOutput(outputId = "Helptext1")
                   
                 ),
                 # Main panel (right side of screen) used for outputs (Not obliged)
                 mainPanel(
                   conditionalPanel(
                     condition = "input.Submit1",
                     
                     htmlOutput(outputId = "Error1"),
                     
                     plotOutput(outputId = "distPlot") %>% withSpinner(),
                     
                     plotOutput(outputId = "distPlot2"),
                     
                     htmlOutput(outputId = "Files1")
                     
                   )
                 )
               )),
      
      # Second tab
      tabPanel(title = "Alpha-Diversity",
               
               sidebarLayout(
                 sidebarPanel(
                   actionButton(
                     inputId = "Submit2",
                     label = "Run alpha-diversity",
                     icon = icon("refresh")
                   ),
                   
                   tags$br(),
                   
                   tags$br(),
                   
                   actionButton(
                     inputId = "Help2",
                     label = "Help",
                     icon = icon("question-circle")
                   ),
                   
                   htmlOutput(outputId = "Helptext2")
                 ),
                 
                 mainPanel(
                   conditionalPanel(condition = "input.Submit2",
                                    
                                    htmlOutput("Files2") %>% withSpinner())
                   
                 )
                 
               )),
      # Third tab
      tabPanel(title = "Beta-Diversity",
               sidebarLayout(
                 sidebarPanel(
                   fileInput(
                     inputId = "input_meta",
                     label = "Metadata file",
                     accept = ".tab"
                   ),
                   fileInput(
                     inputId = "input_tree",
                     label = "Phylogenetic tree",
                     accept = '.tre'
                   ),
                   selectInput(
                     inputId = "group_name",
                     label = "Group from metadat file to analyze",
                     choices = c("")
                   ),
                   radioButtons(
                     inputId = "label_samples",
                     "Label samples:",
                     choiceNames = c(
                       "Samples are not labeled in the MDS/NMDS plots",
                       "All Samples are labed in the MDS/NMDS plots"
                     ),
                     choiceValues = c(0, 1)
                   ),
                   selectizeInput(
                     inputId = "label_id",
                     "Samples that should be labeled",
                     choices = c(),
                     multiple = T
                   ),
                   sliderInput(
                     inputId = "kmers_limit",
                     label = "Perform De-novo clustering for this many samples:",
                     min = 1,
                     max = 10,
                     value = 10,
                     step = 1
                   ),
                   actionButton(
                     inputId = "Submit3",
                     label = 'Run beta-diversity',
                     icon = icon("refresh")
                   ),
                   tags$br(),
                   tags$br(),
                   actionButton(
                     inputId = "Help3",
                     label = "Help",
                     icon = icon("question-circle")
                   ),
                   htmlOutput(outputId = "Helptext3")
                 ),
                 mainPanel(
                   conditionalPanel(
                     condition = "input.Submit3",
                     htmlOutput(outputId = "Error2"),
                     htmlOutput(outputId = "bdtreetext"),
                     plotOutput(outputId = "bdtree") %>% withSpinner(),
                     htmlOutput(outputId = "bdmdstext"),
                     plotOutput(outputId = "bdmds"),
                     htmlOutput(outputId = "bdnmdstext"),
                     plotOutput(outputId = "bdnmds"),
                     htmlOutput(outputId = "pdfviewer1text"),
                     htmlOutput(outputId = 'pdfviewer1'),
                     htmlOutput(outputId = "pdfviewer2text"),
                     htmlOutput(outputId = 'pdfviewer2'),
                     htmlOutput(outputId = "optclustertext"),
                     plotOutput(outputId = 'optcluster'),
                     htmlOutput(outputId = "Files3")
                   )
                 )
               )),
      # Fourth tab
      tabPanel(
        "Taxonomic-Binning",
        sidebarLayout(
          sidebarPanel(
            actionButton(
              inputId = "Submit4",
              label = 'Run taxonomic-binning',
              icon = icon("refresh")
            ),
            tags$br(),
            tags$br(),
            actionButton(
              inputId = "Help4",
              label = "Help",
              icon = icon("question-circle")
            ),
            htmlOutput(outputId = "Helptext4")
          ),
          mainPanel(
            htmlOutput(outputId = "Error3"),
            tabsetPanel(
              type = "tabs",
              tabPanel(
                title = "Kingdom",
                conditionalPanel(condition = "input.Submit4",plotOutput(outputId = "kingdom") %>% withSpinner())
              ),
              tabPanel(
                title = "Phyla",
                conditionalPanel(condition = "input.Submit4", plotOutput(outputId = "phyla") %>% withSpinner())
              ),
              tabPanel(
                title = "Class",
                conditionalPanel(condition = "input.Submit4", plotOutput(outputId = "class") %>% withSpinner())
              ),
              tabPanel(
                title = "Order",
                conditionalPanel(condition = "input.Submit4", plotOutput(outputId = "order") %>% withSpinner())
              ),
              tabPanel(
                title = "Family",
                conditionalPanel(condition = "input.Submit4", plotOutput(outputId = "family") %>% withSpinner())
              ),
              tabPanel(
                title = "Genus",
                conditionalPanel(condition = "input.Submit4", plotOutput(outputId = "genus") %>% withSpinner())
              )
            ),
            htmlOutput(outputId = "Files4")
          )
        )
      ),
      # Fifth tab
      tabPanel(
        "Create-Input-Tables",
        sidebarLayout(
          sidebarPanel(
            actionButton(
              inputId = "Submit5",
              label = "Create input files",
              icon = icon("refresh")
            ),
            tags$br(),
            tags$br(),
            actionButton(
              inputId = "Help5",
              label = "Help",
              icon = icon("question-circle")
            ),
            htmlOutput(outputId = "Helptext5")
          ),
          mainPanel(
            htmlOutput(outputId = "Error4"),
            conditionalPanel(condition = "input.Submit5", htmlOutput(outputId = "Files5") %>% withSpinner())
          )
        )
      ),
      # Sixth tab
      tabPanel(
        "Serial-Group-Comparisons",
        sidebarLayout(
          sidebarPanel(
            selectInput(
              inputId = "sgcgroup",
              "Name of independent variable that analysis will be performed on:",
              choices = c("")
            ),
            uiOutput(outputId = "grorder"),
            numericInput(
              inputId = "abundance_cutoff",
              label = "Cutoff of relative abundance:",
              value = 0.5,
              min = 0,
              max = 100,
              step = 0.1
            ),
            numericInput(
              inputId = "prevalence_cutoff",
              label = "Prevalence cutoff",
              value = 0.3,
              min = 0,
              max = 100,
              step = 0.1
            ),
            numericInput(
              inputId = "max_median_cutoff",
              label = "Minimum median abundance value that must be observed in at least one group before statistical test is performed",
              value = 1,
              min = 0,
              max = 100,
              step = 0.1
            ),
            checkboxInput(
              inputId = "repzero",
              "Replace 0 Value with NA",
              value = TRUE
            ),
            radioButtons(
              "plotOption",
              "Graphical output parameter",
              choiceNames = c(
                "without individual values as dots",
                "with individual values as dots",
                "with individual values as dots and with sample names"
              ),
              choiceValues = c(1, 2, 3)
            ),
            numericInput(
              inputId = "sig.cutoff",
              label = "Significance cutoff level",
              value = 0.05,
              min = 0,
              max = 1,
              step = 0.01
            ),
            actionButton(
              inputId = "Submit6",
              label = "Run Serial group comparison",
              icon = icon("refresh")
            ),
            tags$br(),
            tags$br(),
            actionButton(
              inputId = "Help6",
              label = "Help",
              icon = icon("question-circle")
            ),
            htmlOutput(outputId = "Helptext6")
          ),
          mainPanel(
            conditionalPanel(
              condition = "input.Submit6",
              htmlOutput(outputId = "Error5"),
              htmlOutput(outputId = "box") %>% withSpinner(),
              htmlOutput(outputId = "boxpoint"),
              htmlOutput(outputId = "violin"),
              htmlOutput(outputId = "Files6")
            )
          )
        )
      ),
      
      
      # Not implemented tab
      
      # tabPanel("Over-Time-Serial-Comparisons",
      #          sidebarLayout(
      #            sidebarPanel(
      #              selectInput(inputId = "sgcgroupot","Name of independent variable that analysis will be performed on:",choices = c("")),
      #              uiOutput(outputId = "grorderot"),
      #              numericInput(inputId = "abundance_cutoffot",label="Cutoff of relative abundance:", value = 0.5, min=0,max=100,step = 0.1),
      #              numericInput(inputId = "prevalence_cutoffot",label="Prevalence cutoff", value = 0.3, min=0,max=100,step = 0.1),
      #              numericInput(inputId = "max_median_cutoffot",label="Minimum median abundance value that must be observed in at least one 
      #              group before statistical test is performed", value = 1, min=0,max=100,step = 0.1),
      #              checkboxInput(inputId = "repzeroot","Replace 0 Value with NA",value = TRUE),
      #              radioButtons("plotOptionot","Graphical output parameter",choiceNames = c("without individual values as dots",
      #             "with individual values as dots","with individual values as dots and with sample names"),choiceValues = c(1,2,3)),
      #              numericInput(inputId = "sig.cutoffot",label = "Significance cutoff level",value = 0.05,min = 0,max = 1,step = 0.01),
      #              actionButton(inputId = "Submit7", label= "Run over-time-serial-comparisons", icon=icon("refresh"))
      #            ),
      #            mainPanel(
      #
      #            )
      #          )
      # ),
      
      
      # Seventh tab
      tabPanel("Correlations",
               sidebarLayout(
                 sidebarPanel(
                   numericInput(
                     inputId = "signf_cutoff",
                     label = "cutoff for significance:",
                     value = 0.05,
                     min = 0,
                     max = 1,
                     step = 0.01
                   ),
                   radioButtons(
                     "includeTax",
                     label = "Calculate correlation among taxonomic variables",
                     choiceNames = c(
                       "Calculate correlations within OTUs or taxa",
                       "No test within taxonomic variables"
                     ),
                     choiceValues = c(1, 0),
                     selected = 0
                   ),
                   radioButtons(
                     "includeMeta",
                     label = "Calculate correlation among meta-variables",
                     choiceNames = c(
                       "Calculate correlations within meta-variables",
                       "No test within meta-variables"
                     ),
                     choiceValues = c(1, 0),
                     selected = 0
                   ),
                   radioButtons(
                     "fill_NA",
                     label = "Handling of missing values for meta-variables",
                     choiceNames = c(
                       "Missing values are filled with the mean for the corresponding variable",
                       "NO imputation (replacing missing data with substituted values)"
                     ),
                     choiceValues = c(1, 0),
                     selected = 0
                   ),
                   radioButtons(
                     "replace_zeros",
                     label = "Treat zeros in taxonomic variables as missing values",
                     choiceNames = c(
                       "Consider taxonomic zeros as missing values",
                       "Keep zeros for the calculation of correlations"
                     ),
                     choiceValues = c(1, 0),
                     selected = 1
                   ),
                   numericInput(
                     inputId = "prevalence_exclusion",
                     label = "Cutoff for the minimum number of values (prevalence) for a given taxonomic variable to be considered for calculation",
                     value = 0.3,
                     min = 0,
                     max = 1,
                     step = 0.01
                   ),
                   sliderInput(
                     inputId = "min_pair_support",
                     label = "Cutoff for the minimal number of pairs observations required for calculation of correlations",
                     value = 4,
                     min = 0,
                     max = 10,
                     step = 1
                   ),
                   numericInput(
                     inputId = "plot_pval_cutoff",
                     label = "Set a significance cutoff for graphical output",
                     value = 0.05,
                     min = 0,
                     max = 1,
                     step = 0.01
                   ),
                   numericInput(
                     inputId = "plot_corr_cutoff",
                     label = "Correlation coefficient cutoff for graphical output",
                     value = 0.5,
                     min = 0,
                     max = 1,
                     step = 0.01
                   ),
                   actionButton(
                     inputId = "Submit8",
                     label = "Run correlations",
                     icon = icon("refresh")
                   ),
                   tags$br(),
                   tags$br(),
                   actionButton(
                     inputId = "Help7",
                     label = "Help",
                     icon = icon("question-circle")
                   ),
                   htmlOutput(outputId = "Helptext7")
                 ),
                 mainPanel(
                   conditionalPanel(
                     condition = "input.Submit8",
                     htmlOutput(outputId = "Error6"),
                     plotOutput(outputId = "corrplot") %>% withSpinner(),
                     htmlOutput(outputId = "Finalpdf")
                   )
                 )
               )),
      # Eight tab
      tabPanel(
        "Download",
        downloadButton(outputId = "Download", label = "Download all files as .zip"),
        tags$br(),
        tags$br(),
        actionButton(
          inputId = "Help8",
          label = "Help",
          icon = icon("question-circle")
        ),
        htmlOutput(outputId = "Helptext8")
      )
    )
  )
}

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  ##### Randomized string function #####
  randomString <- function(n = 5000) {
    a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
    paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
  }
  
  
  ##### Each user gets new temp directory, also in www directory #####
  wd <- randomString(1)[1]
  dir.create(wd)
  dir.create(paste0("www/", wd))
  
  
  ##### remove temp directory on end of script #####
  session$onSessionEnded(function() {
    system(paste("rm -Rf", wd))
    system(paste("rm -Rf", paste0("www/", wd)))
  })
  
  ##### Help Messages #####
  
  observeEvent(input$Help1, {
    output$Helptext1 <- renderUI({
      tags$iframe(src = "Normalization Script ReadMe.pdf",
                  height = 600,
                  width = "100%")
    })
  })
  
  observeEvent(input$Help2, {
    output$Helptext2 <- renderUI({
      tags$iframe(src = "Alpha Diversity Script ReadMe.pdf",
                  height = 600,
                  width = "100%")
    })
  })
  
  observeEvent(input$Help3, {
    output$Helptext3 <- renderUI({
      tags$iframe(src = "Beta Diversity Script ReadMe.pdf",
                  height = 600,
                  width = "100%")
    })
  })
  
  observeEvent(input$Help4, {
    output$Helptext4 <- renderUI({
      tags$iframe(src = "Taxonomic Binning Script ReadMe.pdf",
                  height = 600,
                  width = "100%")
    })
  })
  
  observeEvent(input$Help5, {
    output$Helptext5 <- renderText({
      "Intermediate step that creates some files needed for serial group comparisson"
    })
  })
  
  observeEvent(input$Help6, {
    output$Helptext6 <- renderUI({
      tags$iframe(src = "Serial Group Comparisons Script ReadMe.pdf",
                  height = 600,
                  width = "100%")
    })
  })
  
  observeEvent(input$Help7, {
    output$Helptext7 <- renderUI({
      tags$iframe(src = "Correlations Script ReadMe.pdf",
                  height = 600,
                  width = "100%")
    })
  })
  
  observeEvent(input$Help8, {
    #struct <<- dir_tree(wd,recurse = T)
    output$Helptext8 <- renderText({
      bdgroup <-
        ifelse(input$group_name != "",
               input$group_name,
               "{chosen beta-diveristy group name}")
      sgcgroup <-
        ifelse(
          input$sgcgroup != "",
          input$sgcgroup,
          "{chosen serial group comparisson group name}"
        )
      
      paste0(
        #paste(struct,collapse = "<br>"),
        "<pre>Creates a zip file that has 5 folders with the output files created so far.
        
        Folder structure:
        1.Normalization
          . OTUs_Table-norm.tab
          . OTUs_Table-norm-rel.tab
          . OTUs_Table-norm-rel-tax.tab
          . OTUs_Table-norm-tax.tab
          . RarefactionCurve.pdf
          . RarefactionCurve.tab
        2.Alpha-Diversity
          . alpha-diversity.tab
          . OTUs_Table-norm.tab
        3.Beta-Diversity
          . ",bdgroup,
            "* ",bdgroup,"_beta-diversity.pdf
            * ",bdgroup,"_distance-matrix-gunif.tab
            * phylogram.pdf
          . de-novo-clustering.pdf
          . distance-matrix-gunif.tab
          . OTUs_table-norm.tab
          . samples-Tree.nwk
        4.Taxonomic-Binning
          . 0.Kingdom.all.tab
          . 1.Phyla.all.tab
          . 2.Classes.all.tab
          . 3.Orders.all.tab
          . 4.Families.all.tab
          . 5.Genera.all.tab
          . tax.summary.all.tab
          . taxonomic-overview.pdf
        5.Serial-Group-Comparisons
          . ",sgcgroup,"_OTUsCombined_2019-12-10
            * ",sgcgroup,"plotbox.pdf
            * my_analysis_log.txt
            * OTUsCombined-",sgcgroup,"-modified.txt
            * ",sgcgroup,"plotboxpoint.pdf
            * OTUsCombined-",sgcgroup,"-FisherTestAll.tab
            * OTUsCombined-",sgcgroup,"-pvalues.tab
            * ",sgcgroup,"plotviolin.pdf
            * OTUsCombined-",sgcgroup,"-FisherTestPairWise.tab
            * OTUsCombined-",sgcgroup,"-sign_pairs.tab
          . mapping_file.tab
          . OTUsCombined.tab
          . TaxaCombined.tab
          . alpha-diversity.tab
          . OTUs_Table-norm-rel.tab
          . tax.summary.all.tab
        6.Correlations
          . OTUsCombined_Corr_input_table.tab
          . OTUsCombined_Corr_input_table_2019-12-10
            * corrplot.pdf
            * linear_sign_pairs.pdf
            * OTUsCombined_Corr_input_table_2019-12-10
              - correlation-table.tab
              - cutoff-pairs-corr-sign.tab
              - plotted-pairs-stat.tab
              - pval-table.tab
              - support-table.tab
              - transformed.tab</pre>"
      )
    })
    })
  
  ##### Actual script #####
  
  observeEvent(input$Submit1, {
    setwd(wd)
    # Test if file is actually submitted
    if (!is.null(input$OTUtab)) {
      # Create output folders
      tryCatch(error=function(e){
        output$Error1 <- renderUI(
                                  tags$pre(
                                    paste0("Error received:",
                                           as.character(e),
                                           "Something might be wrong with the input file")
                                    ,style="color:red")
                                  )
        output$distPlot <- renderPlot(NULL)
        
        output$distPlot2 <- renderPlot(NULL)
        
        output$Files1 <- renderUI(NULL)
        
        session$sendCustomMessage('deactiveNavs', 'Alpha-Diversity')
        session$sendCustomMessage('deactiveNavs', 'Beta-Diversity')
        session$sendCustomMessage('deactiveNavs', 'Taxonomic-Binnin')
        session$sendCustomMessage('deactiveNavs', 'Create-Input-Tables')
        session$sendCustomMessage('deactiveNavs', 'Serial-Group-Comparisons')
        session$sendCustomMessage('deactiveNavs', 'Correlations')
        
        },
               warning=function(w) NA,
               expr = {
      dir.create("1.Normalization", showWarnings = F)
      dir.create("2.Alpha-Diversity", showWarnings = F)
      dir.create("3.Beta-Diversity", showWarnings = F)
      dir.create("4.Taxonomic-Binning", showWarnings = F)
      dir.create("5.Serial-Group-Comparisons", showWarnings = F)
      dir.create("6.Correlations", showWarnings = F)
      
      ###################       Read input                         ####################
      output$Error1 <- renderText("")
      file_name<-isolate(input$OTUtab$datapath)
      
      ###################       Read all required input files      ####################
      
      # Load the tab-delimited file containing the values to be be checked (rownames in the first column)
      otu_table <-  read.table (
        file_name,
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
      
      # Clean table from empty lines
      otu_table <-
        otu_table[!apply(is.na(otu_table) |
                           otu_table == "", 1, all), ]
      
      ####################       Normalize OTU Table          ###################
      
      
      # Save taxonomy information in vector
      taxonomy <- as.vector(otu_table$taxonomy)
      
      # Delete column with taxonomy information in dataframe
      otu_table$taxonomy <- NULL
      
      # Determine cutoff (Filter for inpossible values)
      labelCutoff <-
        max(1, min(isolate(input$labelCutoff), ncol(otu_table)))
      updateNumericInput(session, "labelCutoff", value = labelCutoff)
      
      # Calculate the minimum sum of all columns/samples
      min_sum <- min(colSums(otu_table))
      
      if (isolate(input$method) == 0) {
        # Divide each value by the sum of the sample and multiply by the minimal sample sum
        norm_otu_table <-
          t(min_sum * t(otu_table) / colSums(otu_table))
      } else {
        # Rarefy the OTU table to an equal sequencing depth
        norm_otu_table <- Rarefy(t(otu_table), depth = min_sum)
        norm_otu_table <-
          t(as.data.frame(norm_otu_table$otu.tab.rff))
      }
      
      # Calculate relative abundances for all OTUs over all samples
      # Divide each value by the sum of the sample and multiply by 100
      rel_otu_table <- t(100 * t(otu_table) / colSums(otu_table))
      
      # Re-insert the taxonomy information in normalized counts table
      norm_otu_table_tax <- cbind(norm_otu_table, taxonomy)
      
      # Reinsert the taxonomy information in relative abundance table
      rel_otu_table_tax <- cbind(rel_otu_table, taxonomy)
      
      ################################################################################
      # Generate a twosided pdf with a rarefaction curve for all samples and a curve
      pdf(file = "1.Normalization/RarefactionCurve.pdf")
      
      # Plot the rarefaction curve for all samples
      rarefactionCurve <- rarecurve(
        data.frame(t(otu_table)),
        step = 20,
        col = "black",
        lty = "solid",
        label = F,
        xlab = "Number of Reads",
        ylab = "Number of Species",
        main = "Rarefaction Curves of All Samples"
      )
      
      # Generate empy vectors for the analysis of the rarefaction curve
      slope = vector()
      SampleID = vector()
      
      # Iterate through all samples
      for (i in seq_along(rarefactionCurve)) {
        # If the sequencing depth is greater 100 the difference between the last and last-100 richness is calcualted
        richness <-
          ifelse(
            length(rarefactionCurve[[i]]) >= 100,
            rarefactionCurve[[i]][length(rarefactionCurve[[i]])] - rarefactionCurve[[i]][length(rarefactionCurve[[i]]) -
                                                                                           100],
            1000
          )
        slope <- c(slope, richness)
        SampleID <- c(SampleID, as.character(names(otu_table)[i]))
      }
      
      # Generate the output table for rarefaction curve
      curvedf <- cbind(SampleID, slope)
      order <- order(curvedf[, 2], decreasing = TRUE)
      # Order the table
      curvedf <- curvedf[order(curvedf[, 2], decreasing = TRUE), ]
      
      # Generates a graph with all samples
      # Underestimated cases are shown in red
      for (i in 1:labelCutoff) {
        N <- attr(rarefactionCurve[[order[i]]], "Subsample")
        lines(N, rarefactionCurve[[order[i]]], col = "red")
      }
      
      # Determine the plotting width and height
      Nmax <-
        sapply(rarefactionCurve, function(x)
          max(attr(x, "Subsample")))
      Smax <- sapply(rarefactionCurve, max)
      
      # Creates an empty plot for rarefaction curves of underestimated cases
      plot(
        c(1, max(Nmax)),
        c(1, max(Smax)),
        xlab = "Number of Reads",
        ylab = "Number of Species",
        type = "n",
        main = paste(labelCutoff, "- most undersampled cases")
      )
      
      for (i in 1:labelCutoff) {
        N <- attr(rarefactionCurve[[order[i]]], "Subsample")
        lines(N, rarefactionCurve[[order[i]]], col = "red")
        text(max(attr(rarefactionCurve[[order[i]]], "Subsample")), max(rarefactionCurve[[order[i]]]), curvedf[i, 1], cex =
               0.6)
      }
      
      dev.off()
      
      #################################################################################
      ######                        Write Output Files                           ######
      #################################################################################
      
      # Write the normalized table in a file and copy in directories alpha-diversity and beta-diversity if existing
      write.table(
        norm_otu_table,
        "1.Normalization/OTUs_Table-norm.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE
      )
      suppressWarnings (try(write.table(
        norm_otu_table,
        "2.Alpha-Diversity/OTUs_Table-norm.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE
      ),
      silent = TRUE)
      )
      suppressWarnings (try(write.table(
        norm_otu_table,
        "3.Beta-Diversity/OTUs_Table-norm.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE)
        ,
        silent =  TRUE)
      )
      
      # Write the normalized table with taxonomy in a file
      write.table(
        norm_otu_table_tax,
        "1.Normalization/OTUs_Table-norm-tax.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE)
      
      # Write the normalized relative abundance table in a file and copy in directory Serial-Group-Comparisons if existing
      write.table(
        rel_otu_table,
        "1.Normalization/OTUs_Table-norm-rel.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE
      )
      suppressWarnings (try(write.table(
        rel_otu_table,
        "5.Serial-Group-Comparisons/OTUs_Table-norm-rel.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE
      ),
      silent = TRUE)
      )
      
      # Write the normalized relative abundance with taxonomy table in a file and copy in directory Taxonomic-Binning if existing
      write.table(
        rel_otu_table_tax,
        "1.Normalization/OTUs_Table-norm-rel-tax.tab",
        sep =  "\t",
        col.names = NA,
        quote = FALSE)
      suppressWarnings (try(write.table(
        rel_otu_table_tax,
        "4.Taxonomic-Binning/OTUs_Table-norm-rel-tax.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE
      ),
      silent = TRUE)
      )
      
      # Write the rarefaction table
      write.table(
        curvedf,
        "1.Normalization/RarefactionCurve.tab",
        sep =  "\t",
        quote = FALSE,
        row.names = FALSE)
      
      
      
      
      
      output$distPlot <-  renderPlot({
        # Plot the rarefaction curve for all samples
        rarefactionCurve <- rarecurve(
          data.frame(t(otu_table)),
          step = 20,
          col = "black",
          lty = "solid",
          label = F,
          xlab = "Number of Reads",
          ylab = "Number of Species",
          main = "Rarefaction Curves of All Samples"
        )
        
        # Generate empy vectors for the analysis of the rarefaction curve
        slope = vector()
        SampleID = vector()
        
        # Iterate through all samples
        for (i in seq_along(rarefactionCurve)) {
          # If the sequencing depth is greater 100 the difference between the last and last-100 richness is calcualted
          richness <-
            ifelse(
              length(rarefactionCurve[[i]]) >= 100,
              rarefactionCurve[[i]][length(rarefactionCurve[[i]])] - rarefactionCurve[[i]][length(rarefactionCurve[[i]]) -
                                                                                             100],
              1000
            )
          slope <- c(slope, richness)
          SampleID <- c(SampleID, as.character(names(otu_table)[i]))
        }
        
        # Generate the output table for rarefaction curve
        curvedf <- cbind(SampleID, slope)
        order <- order(curvedf[, 2], decreasing = TRUE)
        # Order the table
        curvedf <- curvedf[order(curvedf[, 2], decreasing = TRUE), ]
        for (i in 1:labelCutoff) {
          N <- attr(rarefactionCurve[[order[i]]], "Subsample")
          lines(N, rarefactionCurve[[order[i]]], col = "red")
        }
      })
      
      output$distPlot2 <-  renderPlot({
        plot(
          c(1, max(Nmax)),
          c(1, max(Smax)),
          xlab = "Number of Reads",
          ylab = "Number of Species",
          type = "n",
          main = paste(labelCutoff, "- most undersampled cases")
        )
        
        for (i in 1:labelCutoff) {
          N <- attr(rarefactionCurve[[order[i]]], "Subsample")
          lines(N, rarefactionCurve[[order[i]]], col = "red")
          text(max(attr(
            rarefactionCurve[[order[i]]], "Subsample"
          )), max(rarefactionCurve[[order[i]]]), curvedf[i, 1], cex =
            0.6)
        }
        
        
      })
      
      output$Files1 <-
        renderUI({
          HTML(
            "<pre>The script generates five tab-delimited files and one pdf file<br>
            1. Normalized counts with taxonomy information<br>
            2. Normalized counts without taxonomy information<br>
            3. Normalized relative abundances with taxonomy information<br>
            4. Normalized relative abundances without taxonomy information<br>
            5. Rarefaction curves for all samples and the most undersampled ones (default 5 cases) as PDF<br>
            6. Slope of the Rarefaction curve as species per 100 reads</pre>"
          )
        })
      session$sendCustomMessage('activeNavs', 'Alpha-Diversity')
      })
      
    }
    
    # Error message
    if (!flag) {
      stop(
        "
        It was not possible to install all required R libraries properly.
        Please check the installation of all required libraries manually.\n
        Required libaries: GUniFrac"
      )
    }
    setwd("../")
    
    })
  
  observeEvent(input$Submit2, {
    ##################################################################################
    ######                        Diversity Functions                           ######
    ##################################################################################
    setwd(wd)
    
    
    # Calculate the species richness in a sample
    Species.richness <- function(x)
    {
      # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
      count = sum(x[x > 0.5] ^ 0)
      return(count)
    }
    
    # Calculate the Shannon diversity index
    Shannon.entropy <- function(x)
    {
      total = sum(x)
      se = -sum(x[x > 0] / total * log(x[x > 0] / total))
      return(se)
    }
    
    # Calculate the effective number of species for Shannon
    Shannon.effective <- function(x)
    {
      total = sum(x)
      se = round(exp(-sum(x[x > 0] / total * log(x[x > 0] / total))), digits =
                   2)
      return(se)
    }
    
    # Calculate the Simpson diversity index
    Simpson.concentration <- function(x)
    {
      total = sum(x)
      si = sum((x[x > 0] / total) ^ 2)
      return(si)
    }
    
    # Calculate the effective number of species for Simpson
    Simpson.effective <- function(x)
    {
      total = sum(x)
      si = round(1 / sum((x[x > 0] / total) ^ 2), digits = 2)
      return(si)
    }
    
    ##################################################################################
    ######                             Main Script                              ######
    ##################################################################################
    
    # Read a normalized OTU-table without taxonomy
    otu_table <-
      read.table (
        "2.Alpha-Diversity/OTUs_Table-norm.tab",
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1
      )
    
    # Clean table from empty lines
    otu_table <-
      otu_table[!apply(is.na(otu_table) |
                         otu_table == "", 1, all),]
    
    # Order and transpose OTU-table
    my_otu_table <- otu_table[, order(names(otu_table))]
    my_otu_table <- data.frame(t(my_otu_table))
    
    # Apply diversity functions to table
    otus_div_stats <- data.frame(my_otu_table[, 0])
    otus_div_stats$Richness <-
      apply(my_otu_table, 1, Species.richness)
    otus_div_stats$Shannon <-
      apply(my_otu_table, 1, Shannon.entropy)
    otus_div_stats$Shannon.effective <-
      apply(my_otu_table, 1, Shannon.effective)
    otus_div_stats$Simpson <-
      apply(my_otu_table, 1, Simpson.concentration)
    otus_div_stats$Simpson.effective <-
      apply(my_otu_table, 1, Simpson.effective)
    otus_div_stats$Evenness <-
      otus_div_stats$Shannon / log(otus_div_stats$Richness, 2)
    
    
    # Write the results in a file and copy in directory "Serial-Group-Comparisons" if existing
    write.table(
      otus_div_stats,
      "2.Alpha-Diversity/alpha-diversity.tab",
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    suppressWarnings (try(write.table(
      otus_div_stats[c(1, 3, 5)],
      "5.Serial-Group-Comparisons/alpha-diversity.tab",
      sep = "\t",
      col.names = NA,
      quote = FALSE
    ),
    silent = TRUE))
    
    ##################################################################################
    ######                          End of Script                               ######
    ##################################################################################
    output$Files2 <-
      renderUI({
        HTML("<pre>alpha-diversity.tab has been created</pre>")
      })
    
    updateSliderInput(
      session,
      "kmers_limit",
      max = ncol(otu_table),
      value = ncol(otu_table)
    )
    updateSliderInput(session, "min_pair_support", max = max(4, ncol(otu_table)))
    updateSelectizeInput(session, "label_id", choices = colnames(otu_table))
    setwd("../")
    session$sendCustomMessage('activeNavs', 'Beta-Diversity')
  })
  
  observeEvent(input$input_meta, {
    setwd(wd)
    tryCatch(error=function(e){
      output$Error2 <- renderUI(
        tags$pre(
          paste0("Error received:",
                 as.character(e),
                 "Something is wrong with metadata file")
          ,style="color:red")
      )
      output$bdtreetext <- renderUI(NULL)
      output$bdtree <- renderPlot(NULL) 
      output$bdmdstext <- renderUI(NULL) 
      output$bdmds <- renderPlot(NULL) 
      output$bdnmdstext <- renderUI(NULL) 
      output$bdnmds <- renderPlot(NULL) 
      output$pdfviewer1text <- renderUI(NULL)
      output$pdfviewer1 <- renderUI(NULL)
      output$pdfviewer2text <- renderUI(NULL) 
      output$pdfviewer2 <- renderUI(NULL) 
      output$optclustertext <- renderUI(NULL) 
      output$optcluster <- renderPlot(NULL)
      output$Files3 <- renderUI(NULL)
      
      session$sendCustomMessage('deactiveNavs', 'Taxonomic-Binning')
      session$sendCustomMessage('deactiveNavs', 'Create-Input-Tables')
      session$sendCustomMessage('deactiveNavs', 'Serial-Group-Comparisons')
      session$sendCustomMessage('deactiveNavs', 'Correlations')
      
    },
    warning=function(w) NA,
    expr = {
    # Load the mapping file containing individual sample information (sample names in the first column)
    meta_file <-
      read.table (
        file = input$input_meta$datapath,
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
    
    # Save the column names of the mapping file
    mappingVar <- names(meta_file)
    
    # Clean table from empty lines
    meta_file <-
      data.frame(meta_file[!apply(is.na(meta_file) |
                                    meta_file == "", 1, all), ], row.names = row.names(meta_file))
    
    updateSelectInput(session, "group_name", choices = names(Filter(is.factor, meta_file)))
    updateSelectInput(session, "sgcgroup", choices = names(Filter(is.factor, meta_file)))}
    )
    setwd("../")
  })
  
  observeEvent(input$Submit3, {
    setwd(wd)
    tryCatch(error=function(e){
      if (is.null(input$input_meta) | is.null(input$input_tree)) {
        errm <- "You have to upload a metadata file and a phylogenetic tree!"
      } 
      else {
        errm <- "Something might be wrong with one of the input files"
      }
      output$Error2 <- renderUI(
        tags$pre(
          paste0("Error received:",
                 as.character(e),
                 errm)
          ,style="color:red")
      )
      output$bdtreetext <- renderUI(NULL)
      output$bdtree <- renderPlot(NULL) 
      output$bdmdstext <- renderUI(NULL) 
      output$bdmds <- renderPlot(NULL) 
      output$bdnmdstext <- renderUI(NULL) 
      output$bdnmds <- renderPlot(NULL) 
      output$pdfviewer1text <- renderUI(NULL)
      output$pdfviewer1 <- renderUI(NULL)
      output$pdfviewer2text <- renderUI(NULL) 
      output$pdfviewer2 <- renderUI(NULL) 
      output$optclustertext <- renderUI(NULL) 
      output$optcluster <- renderPlot(NULL)
      output$Files3 <- renderUI(NULL)
      
      session$sendCustomMessage('deactiveNavs', 'Taxonomic-Binning')
      session$sendCustomMessage('deactiveNavs', 'Create-Input-Tables')
      session$sendCustomMessage('deactiveNavs', 'Serial-Group-Comparisons')
      session$sendCustomMessage('deactiveNavs', 'Correlations')
      
    },
    warning=function(w) NA,
    expr = {
      
    output$Error2 <- renderUI(NULL)
    
    group_name <- input$group_name
    label_samples <- input$label_samples
    label_id <- input$label_id
    kmers_limit <- input$kmers_limit
    ###################       Read all required input files      ####################
    
    # Load the tab-delimited file containing the values to be analyzed (samples names in the first column)
    otu_file <-
      read.table (
        file = "3.Beta-Diversity/OTUs_Table-norm.tab",
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
    
    # Clean table from empty lines
    otu_file <-
      otu_file[!apply(is.na(otu_file) | otu_file == "", 1, all), ]
    
    # Load the mapping file containing individual sample information (sample names in the first column)
    meta_file <-
      read.table (
        file = input$input_meta$datapath,
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
    
    # Copy to folder 5
    file.copy(
      input$input_meta$datapath,
      "5.Serial-Group-Comparisons/mapping_file.tab",
      overwrite = T
    )
    
    # Save the column names of the mapping file
    mappingVar <- names(meta_file)
    
    # Clean table from empty lines
    meta_file <-
      data.frame(meta_file[!apply(is.na(meta_file) |
                                    meta_file == "", 1, all), ], row.names = row.names(meta_file))
    
    # Load the phylogenetic tree calculated from the OTU sequences
    tree_file <- read.tree(input$input_tree$datapath)
    
    # Create the directory where all output files are saved (is named after the target group name set above for comparisons)
    dir.create(paste0("3.Beta-Diversity/", group_name))
    
    ####################       Calculate beta-diversity          ###################
    
    # OTU-table and mapping file should have the same order and number of sample names
    # Order the OTU-table by sample names (ascending)
    otu_file <- otu_file[, order(names(otu_file))]
    
    # Transpose OTU-table and convert format to a data frame
    otu_file <- data.frame(t(otu_file))
    
    # Order the mapping file by sample names (ascending)
    meta_file <-
      data.frame(meta_file[order(row.names(meta_file)), ], row.names = row.names(meta_file))
    
    # Assign the column names to the modified mapping file
    names(meta_file) <- mappingVar
    
    # Save the position of the target group name in the mapping file
    meta_file_pos <- which(colnames(meta_file) == group_name)
    
    # Select metadata group based on the pre-set group name
    all_groups <- as.factor(meta_file[, meta_file_pos])
    
    # Root the OTU tree at midpoint
    rooted_tree <- midpoint(tree_file)
    
    # Calculate the UniFrac distance matrix for comparing microbial communities
    unifracs <-
      GUniFrac(otu_file, rooted_tree, alpha = c(0.0, 0.5, 1.0))$unifracs
    
    # Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
    unifract_dist <- unifracs[, , "d_0.5"]
    
    ################ Generate tree #######################
    
    # Save the UniFrac output as distance object
    all_dist_matrix <- as.dist(unifract_dist)
    
    # Apply a hierarchical cluster analysis on the distance matrix based on the Ward's method
    all_fit <- hclust(all_dist_matrix, method = "ward")
    
    # Generates a tree from the hierarchically generated object
    tree <- as.phylo(all_fit)
    my_tree_file_name <- paste(group_name, "/phylogram.pdf", sep = "")
    plot_color <- rainbow(length(levels(all_groups)))[all_groups]
    
    # Save the generated phylogram in a pdf file
    pdf(paste0("3.Beta-Diversity/", my_tree_file_name))
    
    # The tree is visualized as a Phylogram color-coded by the selected group name
    plot(
      tree,
      type = "phylogram",
      use.edge.length = TRUE,
      tip.color = (plot_color),
      label.offset = 0.01
    )
    print.phylo(tree)
    axisPhylo()
    tiplabels(pch = 16, col = plot_color)
    dev.off()
    
    output$bdtree <- renderPlot({
      plot(
        tree,
        type = "phylogram",
        use.edge.length = TRUE,
        tip.color = (plot_color),
        label.offset = 0.01
      )
      print.phylo(tree)
      axisPhylo()
      tiplabels(pch = 16, col = plot_color)
    })
    output$bdtreetext <- renderText({
      "<h2>Phylogenetic tree</h2>"
    })
    
    #################            Build NMDS plot           ########################
    
    # Generated figures are saved in a pdf file
    file_name <- paste(group_name, "beta-diversity.pdf", sep = "_")
    pdf(paste("3.Beta-Diversity/", group_name, "/", file_name, sep = ""))
    
    # Calculate the significance of variance to compare multivariate sample means (including two or more dependent variables)
    # Omit cases where there isn't data for the sample (NA)
    all_groups_comp <- all_groups[!is.na(all_groups)]
    unifract_dist_comp <-
      unifract_dist[!is.na(all_groups), !is.na(all_groups)]
    adonis <- adonis(as.dist(unifract_dist_comp) ~ all_groups_comp)
    all_groups_comp <-
      factor(all_groups_comp, levels(all_groups_comp)[unique(all_groups_comp)])
    
    # Calculate and display the MDS plot (Multidimensional Scaling plot)
    s.class(
      cmdscale(unifract_dist_comp, k = 2),
      col = unique(plot_color),
      cpoint =
        2,
      fac = all_groups_comp,
      sub = paste(
        "MDS plot of Microbial Profiles\n(p-value ",
        adonis[[1]][6][[1]][1],
        ")",
        sep = ""
      )
    )
    if (label_samples == 1) {
      lab_samples <- row.names(cmdscale(unifract_dist_comp, k = 2))
      ifelse (label_id != "",
              lab_samples <-
                replace(lab_samples, !(lab_samples %in% label_id), ""),
              lab_samples)
      text(
        cmdscale(unifract_dist_comp, k = 2),
        labels = lab_samples,
        cex = 0.7,
        adj = c(-.1, -.8)
      )
    }
    
    # Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
    meta <- metaMDS(unifract_dist_comp, k = 2)
    s.class(
      meta$points,
      col = unique(plot_color),
      cpoint = 2,
      fac = all_groups_comp,
      sub = paste(
        "metaNMDS plot of Microbial Profiles\n(p-value ",
        adonis[[1]][6][[1]][1],
        ")",
        sep = ""
      )
    )
    if (label_samples == 1) {
      lab_samples <- row.names(meta$points)
      ifelse (label_id != "",
              lab_samples <-
                replace(lab_samples, !(lab_samples %in% label_id), ""),
              lab_samples)
      text(
        meta$points,
        labels = lab_samples,
        cex = 0.7,
        adj = c(-.1, -.8)
      )
    }
    
    #close the pdf file
    dev.off()
    
    output$bdmds <- renderPlot({
      s.class(
        cmdscale(unifract_dist_comp, k = 2),
        col = unique(plot_color),
        cpoint =
          2,
        fac = all_groups_comp,
        sub = paste(
          "MDS plot of Microbial Profiles\n(p-value ",
          adonis[[1]][6][[1]][1],
          ")",
          sep = ""
        )
      )
      if (label_samples == 1) {
        lab_samples <- row.names(cmdscale(unifract_dist_comp, k = 2))
        ifelse (label_id != "",
                lab_samples <-
                  replace(lab_samples, !(lab_samples %in% label_id), ""),
                lab_samples)
        text(
          cmdscale(unifract_dist_comp, k = 2),
          labels = lab_samples,
          cex = 0.7,
          adj = c(-.1, -.8)
        )
      }
    })
    output$bdmdstext <-
      renderText({
        "<h2>Multidimensional scaling plot</h2>"
      })
    
    output$bdnmds <- renderPlot({
      meta <- metaMDS(unifract_dist_comp, k = 2)
      s.class(
        meta$points,
        col = unique(plot_color),
        cpoint = 2,
        fac = all_groups_comp,
        sub = paste(
          "metaNMDS plot of Microbial Profiles\n(p-value ",
          adonis[[1]][6][[1]][1],
          ")",
          sep = ""
        )
      )
      if (label_samples == 1) {
        lab_samples <- row.names(meta$points)
        ifelse (label_id != "",
                lab_samples <-
                  replace(lab_samples, !(lab_samples %in% label_id), ""),
                lab_samples)
        text(
          meta$points,
          labels = lab_samples,
          cex = 0.7,
          adj = c(-.1, -.8)
        )
      }
    })
    output$bdnmdstext <-
      renderText({
        "<h2>Non-metric multidimensional scaling plot</h2>"
      })
    
    ###############          NMDS for pairwise analysis        ###################
    
    # This plot is only generated if there are more than two groups included in the comparison
    # Calculate the pairwise significance of variance for group pairs
    # Get all groups contained in the mapping file
    unique_groups <- levels(all_groups_comp)
    if (dim(table(unique_groups)) > 2) {
      # Initialise vector and lists
      pVal = NULL
      pairedMatrixList <- list(NULL)
      pair_1_list <- NULL
      pair_2_list <- NULL
      
      for (i in 1:length(combn(unique_groups, 2)[1, ])) {
        # Combine all possible pairs of groups
        pair_1 <- combn(unique_groups, 2)[1, i]
        pair_2 <- combn(unique_groups, 2)[2, i]
        
        # Save pairs information in a vector
        pair_1_list[i] <- pair_1
        pair_2_list[i] <- pair_2
        
        # Generate a subset of all samples within the mapping file related to one of the two groups
        inc_groups <-
          rownames(subset(meta_file, meta_file[, meta_file_pos] == pair_1
                          |
                            meta_file[, meta_file_pos] == pair_2))
        
        # Convert UniFrac distance matrix to data frame
        paired_dist <- as.data.frame(unifract_dist_comp)
        
        # Save all row names of the mapping file
        row_names <- rownames(paired_dist)
        
        # Add row names to the distance matrix
        paired_dist <- cbind(row_names, paired_dist)
        
        # Generate distance matrix with samples of the compared groups (column-wise)
        paired_dist <-
          paired_dist[sapply(paired_dist[, 1], function(x)
            all(x %in% inc_groups)), ]
        
        # Remove first column with unnecessary group information
        paired_dist[, 1] <- NULL
        paired_dist <- rbind(row_names, paired_dist)
        
        # Generate distance matrix with samples of the compared group (row-wise)
        paired_dist <-
          paired_dist[, sapply(paired_dist[1, ], function(x)
            all(x %in% inc_groups))]
        
        # Remove first row with unnecessary group information
        paired_dist <- paired_dist[-1, ]
        
        # Convert generated distance matrix to data type matrix (needed by multivariate analysis)
        paired_matrix <- as.matrix(paired_dist)
        class(paired_matrix) <- "numeric"
        
        # Save paired matrix in list
        pairedMatrixList[[i]] <- paired_matrix
        
        # Applies multivariate analysis to a pair out of the selected groups
        adonis <-
          adonis(paired_matrix ~ all_groups_comp[all_groups_comp == pair_1 |
                                                   all_groups_comp == pair_2])
        
        # List p-values
        pVal[i] <- adonis[[1]][6][[1]][1]
        
      }
      
      # Adjust p-values for multiple testing according to Benjamini-Hochberg method
      pVal_BH <- p.adjust(pVal, method = "BH", n = length(pVal))
      
      # Generated NMDS plots are stored in one pdf file called "pairwise-beta-diversity-nMDS.pdf"
      file_name <-
        paste(group_name, "pairwise-beta-diversity.pdf", sep = "_")
      pdf(paste("3.Beta-Diversity/", group_name, "/", file_name, sep = ""))
      
      for (i in 1:length(combn(unique_groups, 2)[1, ])) {
        meta <- metaMDS(pairedMatrixList[[i]], k = 2)
        s.class(
          meta$points,
          col = rainbow(length(levels(
            all_groups_comp
          ))),
          cpoint = 2,
          fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                            all_groups_comp == pair_2_list[i]]),
          sub = paste(
            "NMDS plot of Microbial Profiles\n ",
            pair_1_list[i],
            " - ",
            pair_2_list[i],
            "\n(p-value ",
            pVal[i],
            ",",
            " corr. p-value ",
            pVal_BH[i],
            ")",
            sep = ""
          )
        )
      }
      dev.off()
      file.copy(
        paste(
          "3.Beta-Diversity/",
          group_name,
          "/",
          paste(group_name, "pairwise-beta-diversity.pdf", sep = "_"),
          sep = ""
        ),
        paste(
          "../www/",
          wd,
          "/",
          paste(group_name, "pairwise-beta-diversity.pdf", sep = "_"),
          sep = ""
        ),
        overwrite = T
      )
      output$pdfviewer1 <- renderUI({
        tags$iframe(
          src = paste0(
            wd,
            "/",
            paste(group_name, "pairwise-beta-diversity.pdf", sep = "_")
          ),
          height = "600",
          width = "100%"
        )
      })
      output$pdfviewer1text <-
        renderText({
          "<h2>Pairwise beta-diversity non-metric MDS</h2>"
        })
      # Generated MDS plots are stored in one pdf file called "pairwise-beta-diversity-MDS.pdf"
      file_name <-
        paste(group_name, "pairwise-beta-diversity-MDS.pdf", sep = "_")
      pdf(paste("3.Beta-Diversity/", group_name, "/", file_name, sep = ""))
      
      for (i in 1:length(combn(unique_groups, 2)[1, ])) {
        # Calculate and display the MDS plot (Multidimensional Scaling plot)
        s.class(
          cmdscale(pairedMatrixList[[i]], k = 2),
          col = rainbow(length(levels(
            all_groups_comp
          ))),
          cpoint =
            2,
          fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                            all_groups_comp == pair_2_list[i]]),
          sub = paste(
            "MDS plot of Microbial Profiles\n ",
            pair_1_list[i],
            " - ",
            pair_2_list[i],
            "\n(p-value ",
            pVal[i],
            ",",
            " corr. p-value ",
            pVal_BH[i],
            ")",
            sep = ""
          )
        )
      }
      dev.off()
      file.copy(
        paste(
          "3.Beta-Diversity/",
          group_name,
          "/",
          paste(group_name, "pairwise-beta-diversity-MDS.pdf", sep = "_"),
          sep = ""
        ),
        paste(
          "../www/",
          wd,
          "/",
          paste(group_name, "pairwise-beta-diversity-MDS.pdf", sep = "_"),
          sep = ""
        ),
        overwrite = T
      )
      output$pdfviewer2 <- renderUI({
        tags$iframe(
          src = paste0(
            wd,
            "/",
            paste(group_name, "pairwise-beta-diversity-MDS.pdf", sep = "_")
          ),
          height = "600",
          width = "100%"
        )
      })
      output$pdfviewer2text <-
        renderText({
          "<h2>Pairwise beta-diversity MDS</h2>"
        })
      
    }
    
    ######                        Determine number of clusters                           ######
    
    nclusters = NULL
    if (dim(otu_file)[1] - 1 <= kmers_limit) {
      kmers_limit = dim(otu_file)[1] - 1
    }
    for (k in 1:kmers_limit) {
      if (k == 1) {
        nclusters[k] = NA
      } else {
        # Partitioning the data into k clusters (max k is number of samples within the dataset)
        data_cluster = as.vector(pam(as.dist(unifract_dist_comp), k, diss =
                                       TRUE)$clustering)
        
        # Calculate Calinski-Harabasz Index
        nclusters[k] = calinhara(otu_file, data_cluster, k)
      }
    }
    
    # Generated plot showing the optimal number of clusters
    pdf("3.Beta-Diversity/de-novo-clustering.pdf")
    
    plot(
      nclusters,
      type = "h",
      xlab = "k clusters",
      ylab = "CH index",
      main = "Optimal number of clusters"
    )
    
    dev.off()
    
    output$optcluster <- renderPlot({
      plot(
        nclusters,
        type = "h",
        xlab = "k clusters",
        ylab = "CH index",
        main = "Optimal number of clusters"
      )
    })
    output$optclustertext <-
      renderText({
        "<h2>Optimal cluster plot</h2>"
      })
    
    #################################################################################
    ######                        Write Output Files                           ######
    #################################################################################
    
    # Write the distance matrix table in a file
    file_name <-
      paste(group_name, "distance-matrix-gunif.tab", sep = "_")
    write.table(
      unifract_dist_comp,
      paste("3.Beta-Diversity/", group_name, "/", file_name, sep = ""),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    write.table(
      unifract_dist_comp,
      "3.Beta-Diversity/distance-matrix-gunif.tab",
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    write.tree(tree, "3.Beta-Diversity/samples-Tree.nwk", tree.names = FALSE)
    
    # Graphical output files are generated in the main part of the script
    if (!flag) {
      stop(
        "
        It was not possible to install all required R libraries properly.
        Please check the installation of all required libraries manually.\n
        Required libaries:ade4, GUniFrac, phangorn"
      )
    }
    output$Files3 <-
      renderUI({
        HTML(
          "<h2>Output files:</h2>
          <pre>The script generates three graphical outputs (pdf), one text file and a newick tree
          1. A phylogram with colour-coded group clustering
          2. MDS and NMDS plots showing information about beta-diversity across all sample groups
          3. MDS and NMDS plots of all pairwise comparisons
          4. The distance matrix
          5. Plot showing the optimal number of clusters
          6. Dendogram for all samples in a newick tree file</pre>"
        )
      })
    session$sendCustomMessage('activeNavs', 'Taxonomic-Binning')
    
    }
  )
    setwd("../")
    
    #################################################################################
    ######                           End of Script                             ######
    #################################################################################
    
      })
  
  observeEvent(input$Submit4, {
    ##################################################################################
    ######                             Main Script                              ######
    ##################################################################################
    
    ###################            Read input table              ####################
    setwd(wd)
    
    tryCatch(error=function(e){
      
      output$Error3 <- renderUI(
        tags$pre(
          paste0("Error received:",
                 as.character(e),
                 "Something went wrong in previous steps... Probably an error in the taxonomy.")
                 ,style="color:red")
      )
      output$kingdom <- renderPlot(NULL)
      output$phyla <- renderPlot(NULL)
      output$class <- renderPlot(NULL)
      output$order <- renderPlot(NULL)
      output$family <- renderPlot(NULL)
      output$genus <- renderPlot(NULL)
      output$Files4 <- renderUI(NULL)
      
      session$sendCustomMessage('deactiveNavs', 'Create-Input-Tables')
      session$sendCustomMessage('deactiveNavs', 'Serial-Group-Comparisons')
      session$sendCustomMessage('deactiveNavs', 'Correlations')
      
    },
    warning=function(w) NA,
    expr = {
    # Load the tab-delimited file containing the abundances and taxonomic information to be checked (rownames in the first column)
    otu_table <-
      read.table (
        "4.Taxonomic-Binning/OTUs_Table-norm-rel-tax.tab",
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
    
    # Clean table from empty lines
    otu_table <-
      otu_table[!apply(is.na(otu_table) | otu_table == "", 1, all), ]
    
    # Create a dataframe with a number of rows identical to the number of OTUs in the dataset
    taxonomy <- otu_table[, dim(otu_table)[2]]
    
    # Test if the taxonomy column is in the correct format (delimited by semicolon)
    if (all(grepl("(?:[^;]*;){6}", taxonomy)) == FALSE) {
      #Send error message if taxonomy is not in the right format
      stop(
        "Wrong number of taxonomic classes\n
        
        Taxonomic levels have to be separated by semicolons (six in total).
        IMPORTANT: if taxonomic information at any level is missing, the semicolons are still needed:\n
        
        e.g.Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotella;
        e.g.Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales;Prevotellaceae;;"
      )
    } else {
      # Delete the taxonomy row from the OTU-table
      otuFile <- otu_table[, c(1:dim(otu_table)[2] - 1)]
      
      # Initialize empty dataframe
      taxonomy_new <- NULL
      
      for (i in 1:dim(otu_table)[1]) {
        # Split taxonomic information in its taxonomic classes
        # Kingdom - Phylum - Class - Family - Order - Genus
        splitTax <- strsplit(x = as.character(taxonomy[i]), ";")
        
        # Save the position where the first empty string (sequence of characters) occurs
        value <- which(splitTax[[1]] == "")[1]
        
        # Save the last known taxa information
        lastTaxa = splitTax[[1]][value - 1]
        
        # Replace all empty values by the last taxa information and the prefix "unkown_"
        splitTax <-
          replace(splitTax[[1]], splitTax[[1]] == "", paste("unknown_", lastTaxa))
        
        # Write new taxonomic information in the dataframe
        taxonomy_new[i] <- list(splitTax)
      }
      
      # Adjust dataframe with modified taxonomic information
      taxonomy_new <- t(as.data.frame(taxonomy_new))
      row.names(taxonomy_new) <- row.names(otuFile)
      
      # Add level information to all taxonomies
      # For taxonomies related to kingdom level
      taxonomy_new[, 1] <- sub("^", "k__", taxonomy_new[, 1])
      
      # For taxonomies related to phylum level
      taxonomy_new[, 2] <- sub("^", "p__", taxonomy_new[, 2])
      
      # For taxonomies related to class level
      taxonomy_new[, 3] <- sub("^", "c__", taxonomy_new[, 3])
      
      # For taxonomies related to order level
      taxonomy_new[, 4] <- sub("^", "o__", taxonomy_new[, 4])
      
      # For taxonomies related to family level
      taxonomy_new[, 5] <- sub("^", "f__", taxonomy_new[, 5])
      
      # For taxonomies related to genus level
      taxonomy_new[, 6] <- sub("^", "g__", taxonomy_new[, 6])
      
      #################################################################################
      
      # Create list with taxonomic information for each taxonomy level
      class_list <-
        list(
          unique(taxonomy_new[, 1]),
          unique(taxonomy_new[, 2]),
          unique(taxonomy_new[, 3]),
          unique(taxonomy_new[, 4]),
          unique(taxonomy_new[, 5]),
          unique(taxonomy_new[, 6])
        )
      
      # Clone the created list for further processing
      sample_list <- class_list
      list_length <- NULL
      
      # Iterate through all six taxonomy levels
      for (a in 1:6) {
        lis <- lapply(class_list[a], lapply, length)
        names(lis) <- lapply(class_list[a], length)
        
        # Individual number of taxonomies for each taxonomic level
        num_taxa <- as.integer(names(lis))
        list_length[a] <- num_taxa
        
        # Iterate through taxonomic class specific taxonomies
        for (b  in 1:num_taxa) {
          # Initialize list with the value zero for all taxonomies
          sample_list[[a]][[b]] <- list(rep.int(0, dim(otuFile)[2]))
          
        }
      }
      
      #################################################################################
      #################################################################################
      # Save relative abundances of all samples for each taxonomy
      
      # Iterate through all OTUs
      for (i in 1:dim(otu_table)[1]) {
        # Iterate through all taxonomic levels
        for (m in 1:6) {
          # List of m-th taxonomies of i-th taxonomic levels (e.g. m = Kingdom, i = 4th OTU -> Clostridiales)
          taxa_in_list <- list(taxonomy_new[i, ])[[1]][m]
          
          # Record the current position in a list
          position <- which(class_list[[m]] == taxa_in_list)
          
          # All rows with taxonomic information of n-th sample
          matrix <- data.matrix(otuFile)
          sub_sample_tax <-
            (subset(matrix, taxonomy_new[, m] == taxa_in_list))
          
          # Get the actual value out of the list (initialized with zero)
          temp <- unlist(sample_list[[m]][[position]])
          
          # Calculate the summed up relative abundances for the particular taxonomic class for n-th sample
          temp <- colSums(sub_sample_tax)
          
          # Replace values by new summed values
          sample_list[[m]][[position]] <- list(temp)
          
        }
      }
      
      #################################################################################
      ######                         Write output                                ######
      #################################################################################
      
      # Generate tables for each taxonomic class
      
      ##Kingdom table
      # Create table with taxonomic information (kingdom level)
      kingdom <-
        matrix(
          unlist(sample_list[[1]]),
          nrow = dim(otuFile)[2],
          ncol = list_length[1],
          dimnames = list(names(otuFile), unlist(class_list[[1]]))
        )
      kingdom <- (t(kingdom))
      
      ##Phylum table
      # Create table with taxonomic information (phylum level)
      phyla <-
        matrix(
          unlist(sample_list[[2]]),
          nrow = dim(otuFile)[2],
          ncol = list_length[2],
          dimnames = list(names(otuFile), unlist(class_list[[2]]))
        )
      phyla <- (t(phyla))
      
      # Order table according to taxonomic name (descending)
      phyla <- phyla[order(row.names(phyla)), ]
      
      ## Class table
      # Create table with taxonomic information (class level)
      classes <-
        matrix(
          unlist(sample_list[[3]]),
          nrow = dim(otuFile)[2],
          ncol = list_length[3],
          dimnames = list(names(otuFile), unlist(class_list[[3]]))
        )
      classes <- (t(classes))
      
      # Order dataframe according to taxonomic name (descending)
      classes <- classes[order(row.names(classes)), ]
      
      ## Orders
      # create table with taxonomic information (Order)
      orders <-
        matrix(
          unlist(sample_list[[4]]),
          nrow = dim(otuFile)[2],
          ncol = list_length[4],
          dimnames = list(names(otuFile), unlist(class_list[[4]]))
        )
      orders <- (t(orders))
      
      # Order dataframe according to taxonomic name (descending)
      orders <- orders[order(row.names(orders)), ]
      
      ## Family table
      # Create table with taxonomic information (family level)
      families <-
        matrix(
          unlist(sample_list[[5]]),
          nrow = dim(otuFile)[2],
          ncol = list_length[5],
          dimnames = list(names(otuFile), unlist(class_list[[5]]))
        )
      families <- (t(families))
      
      # Order dataframe according to taxonomic name (descending)
      families <- families[order(row.names(families)), ]
      
      ## Genus level
      # Create table with taxonomic information (generum level)
      genera <-
        matrix(
          unlist(sample_list[[6]]),
          nrow = dim(otuFile)[2],
          ncol = list_length[6],
          dimnames = list(names(otuFile), unlist(class_list[[6]]))
        )
      genera <- (t(genera))
      
      # Order dataframe according to taxonomic name (descending)
      genera <- genera[order(row.names(genera)), ]
      
      # Merge all dataframes
      tax_summary <-
        rbind.data.frame(kingdom, phyla, classes, orders, families, genera)
      
      # Identify duplicates and remove them
      tax_summary <-
        tax_summary[!duplicated(row.names(tax_summary)), ]
      
      ################################################################################
      ######                        Write Output Files                           ######
      #################################################################################
      
      
      # Write output files for taxonomic composition of every sample
      write.table(
        kingdom,
        "4.Taxonomic-Binning/0.Kingdom.all.tab",
        sep = "\t",
        col.names = NA
      )
      write.table(
        phyla,
        "4.Taxonomic-Binning/1.Phyla.all.tab",
        sep = "\t",
        col.names = NA
      )
      write.table(
        classes,
        "4.Taxonomic-Binning/2.Classes.all.tab",
        sep = "\t",
        col.names = NA
      )
      write.table(
        orders,
        "4.Taxonomic-Binning/3.Orders.all.tab",
        sep = "\t",
        col.names = NA
      )
      write.table(
        families,
        "4.Taxonomic-Binning/4.Families.all.tab",
        sep = "\t",
        col.names = NA
      )
      write.table(
        genera,
        "4.Taxonomic-Binning/5.Genera.all.tab",
        sep = "\t",
        col.names = NA
      )
      write.table(
        tax_summary,
        "4.Taxonomic-Binning/tax.summary.all.tab",
        sep = "\t",
        col.names = NA
      )
      suppressWarnings (try(write.table(
        tax_summary,
        "5.Serial-Group-Comparisons/tax.summary.all.tab",
        sep = "\t",
        col.names = NA,
        quote = FALSE
      ),
      silent = TRUE)
      )
      
      #################################################################################
      ######                        Write Graphical Output                       ######
      #################################################################################
      
      pdf("4.Taxonomic-Binning/taxonomic-overview.pdf")
      par(xpd  =  T, mar  =  par()$mar  +  c(0, 0, 0, 9))
      
      #Kingdom
      #k_col=distinctColorPalette(dim(kingdom)[1])
      k_col  =  rainbow(dim(kingdom)[1])
      k_col  =  sample(k_col)
      barplot(
        kingdom,
        col  =  k_col,
        cex.names  =  0.5,
        ylab  =  "cumulative relative abundance (%)",
        las  =  2,
        main  =  "Taxonomic binning at Kingdom level"
      )
      legend(
        par('usr')[2],
        par('usr')[4],
        bty = 'n',
        rev(row.names(kingdom)),
        cex = 0.7,
        col = rev(k_col),
        pch = 16,
        pt.cex = 1.2
      )
      
      #Phyla
      #p_col=distinctColorPalette(dim(phyla)[1])
      p_col = rainbow(dim(phyla)[1])
      p_col = sample(p_col)
      barplot(
        phyla,
        col = p_col,
        cex.names = 0.5,
        ylab = "cumulative relative abundance (%)",
        las = 2,
        main = "Taxonomic binning at Phyla level"
      )
      legend(
        par('usr')[2],
        par('usr')[4],
        bty = 'n',
        rev(row.names(phyla)),
        cex = 0.7,
        col = rev(p_col),
        pch = 16,
        pt.cex = 1.2
      )
      
      #Classes
      c_col = rainbow(dim(classes)[1])
      c_col = sample(c_col)
      barplot(
        classes,
        col = c_col,
        cex.names = 0.5,
        ylab = "cumulative relative abundance (%)",
        las = 2,
        main = "Taxonomic binning at Class level"
      )
      legend(
        par('usr')[2],
        par('usr')[4],
        bty = 'n',
        rev(row.names(classes)),
        cex = 0.7,
        col = rev(c_col),
        pch = 16,
        pt.cex = 1.2
      )
      
      #Orders
      o_col = rainbow(dim(orders)[1])
      o_col = sample(o_col)
      barplot(
        orders,
        col = o_col,
        cex.names = 0.5,
        ylab = "cumulative relative abundance (%)",
        las = 2,
        main = "Taxonomic binning at Order level"
      )
      legend(
        par('usr')[2],
        par('usr')[4],
        bty = 'n',
        rev(row.names(orders)),
        cex = 0.7,
        col = rev(o_col),
        pch = 16,
        pt.cex = 1.2
      )
      
      #Families
      f_col = rainbow(dim(families)[1])
      f_col = sample(f_col)
      barplot(
        families,
        col = f_col,
        cex.names = 0.5,
        ylab = "cumulative relative abundance (%)",
        las = 2,
        main = "Taxonomic binning at Family level"
      )
      legend(
        par('usr')[2],
        par('usr')[4],
        bty = 'n',
        rev(row.names(families)),
        cex = 0.7,
        col = rev(f_col),
        pch = 16,
        pt.cex = 1.2
      )
      
      #Genera
      g_col = rainbow(dim(genera)[1])
      g_col = sample(g_col)
      barplot(
        genera,
        col = g_col,
        cex.names = 0.5,
        ylab = "cumulative relative abundance (%)",
        las = 2,
        main = "Taxonomic binning at Genus level"
      )
      legend(
        par('usr')[2],
        par('usr')[4],
        bty = 'n',
        rev(row.names(genera)),
        cex = 0.7,
        col = rev(g_col),
        pch = 16,
        pt.cex = 1.2
      )
      
      dev.off()
      
      ###### Shiny plot outputs
      
      output$kingdom <- renderPlot({
        par(xpd = T, mar = par()$mar + c(0, 0, 0, 9))
        barplot(
          kingdom,
          col = k_col,
          cex.names = 0.5,
          ylab = "cumulative relative abundance (%)",
          las = 2,
          main = "Taxonomic binning at Kingdom level"
        )
        legend(
          par('usr')[2],
          par('usr')[4],
          bty = 'n',
          rev(row.names(kingdom)),
          cex = 0.7,
          col = rev(k_col),
          pch = 16,
          pt.cex = 1.2
        )
      })
      
      output$phyla <- renderPlot({
        par(xpd = T, mar = par()$mar + c(0, 0, 0, 9))
        barplot(
          phyla,
          col = p_col,
          cex.names = 0.5,
          ylab = "cumulative relative abundance (%)",
          las = 2,
          main = "Taxonomic binning at Phyla level"
        )
        legend(
          par('usr')[2],
          par('usr')[4],
          bty = 'n',
          rev(row.names(phyla)),
          cex = 0.7,
          col = rev(p_col),
          pch = 16,
          pt.cex = 1.2
        )
      })
      
      output$class <- renderPlot({
        par(xpd = T, mar = par()$mar + c(0, 0, 0, 9))
        barplot(
          classes,
          col = c_col,
          cex.names = 0.5,
          ylab = "cumulative relative abundance (%)",
          las = 2,
          main = "Taxonomic binning at Class level"
        )
        legend(
          par('usr')[2],
          par('usr')[4],
          bty = 'n',
          rev(row.names(classes)),
          cex = 0.7,
          col = rev(c_col),
          pch = 16,
          pt.cex = 1.2
        )
      })
      
      output$order <- renderPlot({
        par(xpd = T, mar = par()$mar + c(0, 0, 0, 9))
        barplot(
          orders,
          col = o_col,
          cex.names = 0.5,
          ylab = "cumulative relative abundance (%)",
          las = 2,
          main = "Taxonomic binning at Order level"
        )
        legend(
          par('usr')[2],
          par('usr')[4],
          bty = 'n',
          rev(row.names(orders)),
          cex = 0.7,
          col = rev(o_col),
          pch = 16,
          pt.cex = 1.2
        )
      })
      
      output$family <- renderPlot({
        par(xpd = T, mar = par()$mar + c(0, 0, 0, 9))
        barplot(
          families,
          col = f_col,
          cex.names = 0.5,
          ylab = "cumulative relative abundance (%)",
          las = 2,
          main = "Taxonomic binning at Family level"
        )
        legend(
          par('usr')[2],
          par('usr')[4],
          bty = 'n',
          rev(row.names(families)),
          cex = 0.7,
          col = rev(f_col),
          pch = 16,
          pt.cex = 1.2
        )
      })
      
      output$genus <- renderPlot({
        par(xpd = T, mar = par()$mar + c(0, 0, 0, 9))
        barplot(
          genera,
          col = g_col,
          cex.names = 0.5,
          ylab = "cumulative relative abundance (%)",
          las = 2,
          main = "Taxonomic binning at Genus level"
        )
        legend(
          par('usr')[2],
          par('usr')[4],
          bty = 'n',
          rev(row.names(genera)),
          cex = 0.7,
          col = rev(g_col),
          pch = 16,
          pt.cex = 1.2
        )
      })
      
      output$Files4 <-
        renderText({
          "<pre>The script is generating seven tab-delimited files
          1. Relative taxonomic abundance at the kingdom level for each sample
          2. Relative taxonomic abundance at the phyla level for each sample
          3. Relative taxonomic abundance at the class level for each sample
          4. Relative taxonomic abundance at the order level for each sample
          5. Relative taxonomic abundance at the family level for each sample
          6. Relative taxonomic abundance at the genera level for each sample
          7. Relative taxonomic abundance at all taxonomic levels for each sample</pre>"
        })
    
      
      session$sendCustomMessage('activeNavs', 'Create-Input-Tables')
    
    }
    })
    setwd("../")
    #################################################################################
    ######                           End of Script                             ######
    #################################################################################
    })
  
  observeEvent(input$Submit5, {
    setwd(wd)
    tryCatch(error=function(e){
      
      output$Error4 <- renderUI(
        tags$pre(
          paste0("Error received:",
                 as.character(e),
                 "Something went wrong in previous steps...")
          ,style="color:red")
      )
      output$Files5 <- renderUI(NULL)
      
      session$sendCustomMessage('deactiveNavs', 'Serial-Group-Comparisons')
      session$sendCustomMessage('deactiveNavs', 'Correlations')
      
    },
    warning=function(w) NA,
    expr = {
    ##################            Input Files                   #####################
    #' Please give the name of the file with alpha-diversity measures
    alpha <- "5.Serial-Group-Comparisons/alpha-diversity.tab"
    
    #' Please give the name of the file with OTUs relative abundance
    RelativeAbundanceOTUs <-
      "5.Serial-Group-Comparisons/OTUs_Table-norm-rel.tab"
    
    
    #' Please give the name of the file with relative abundances of different taxonomic levels
    TaxanomyAll <- "5.Serial-Group-Comparisons/tax.summary.all.tab"
    
    
    #' Please give the name of the meta file with sample groups and additional metadata variables if available
    MetaFile <- "5.Serial-Group-Comparisons/mapping_file.tab"
    
    
    ###################            Read input table              ####################
    
    # Reading Alpha diversity file
    alpha <-
      read.table(
        file = alpha,
        header = TRUE,
        sep = "\t",
        row.names = 1,
        check.names = F
      )
    
    # Clean table from empty lines
    alpha <- alpha[!apply(is.na(alpha) | alpha == "", 1, all), ]
    
    # Read OTU table
    RelativeAbundanceOTUs <-
      as.data.frame(t(
        read.table(
          file = RelativeAbundanceOTUs,
          header = TRUE,
          sep = "\t",
          row.names = 1,
          check.names = F
        )
      ))
    
    # Clean table from empty lines
    RelativeAbundanceOTUs <-
      RelativeAbundanceOTUs[!apply(is.na(RelativeAbundanceOTUs) |
                                     RelativeAbundanceOTUs == "",
                                   1,
                                   all), ]
    
    # Read Mapping file
    MetaFile <-
      read.table(
        file = MetaFile,
        header = TRUE,
        sep = "\t",
        comment.char = "",
        row.names = 1,
        check.names = F
      )
    
    # Clean table from empty lines
    MetaFile <-
      MetaFile[!apply(is.na(MetaFile) | MetaFile == "", 1, all), ]
    
    # Read taxonomy file
    TaxanomyAll <-
      read.table(
        file = TaxanomyAll,
        header = TRUE,
        sep = "\t",
        row.names = NULL,
        check.names = F
      )
    
    # Clean table from empty lines
    TaxanomyAll <-
      TaxanomyAll[!apply(is.na(TaxanomyAll) | TaxanomyAll == "", 1, all), ]
    
    ColnameTo_assign <- TaxanomyAll[, 1]
    TaxanomyAll[, 1] <- NULL
    TaxanomyAll <- as.data.frame(t(TaxanomyAll))
    colnames(TaxanomyAll) <- ColnameTo_assign
    
    ######################          MAIN PROGRAM                #####################
    
    # Preparing TAXA table:
    combine_taxa <-
      cbind.data.frame(MetaFile[rownames(TaxanomyAll), ], alpha[rownames(TaxanomyAll), ], TaxanomyAll) # merging Meta+Alpha+Taxa based on same order row.
    combine_taxa$SampleID <- row.names(combine_taxa)
    combine_taxa <-
      combine_taxa[, c(ncol(combine_taxa), 1:(ncol(combine_taxa) - 1))]# replacing columnn at the beginging.
    
    # Prepare OTU table:
    combine_OTUs <-
      cbind.data.frame(MetaFile[rownames(RelativeAbundanceOTUs), ], alpha[rownames(RelativeAbundanceOTUs), ], RelativeAbundanceOTUs) # merging Meta+Alpha+OTUs based on same order row.
    combine_OTUs$SampleID <- row.names(combine_OTUs)
    colnames(combine_OTUs)
    combine_OTUs <-
      combine_OTUs[, c(ncol(combine_OTUs), 1:(ncol(combine_OTUs) - 1))]# replacing columnn at the beginging.
    
    #################################################################################
    ######                        Write Output Files                           ######
    #################################################################################
    
    # Writing tables
    if (compareIgnoreOrder(row.names(TaxanomyAll), row.names((MetaFile)))$result &
        compareIgnoreOrder(row.names(alpha), row.names((MetaFile)))$result &
        compareIgnoreOrder(row.names(RelativeAbundanceOTUs), row.names((MetaFile)))$result) {
      write.table(
        combine_taxa,
        file = "5.Serial-Group-Comparisons/TaxaCombined.tab",
        sep = "\t",
        row.names = FALSE
      )
      write.table(
        combine_OTUs,
        file = "5.Serial-Group-Comparisons/OTUsCombined.tab",
        sep = "\t",
        row.names = FALSE
      )
      message("Files are combined successfully")
    } else{
      stop(
        "ATTENTION !!!! Sample names differ across files. Script aborted.",
        "Please ensure that identical sample names are used."
      )
    }
    
    if (!flag) {
      stop(
        "
        It was not possible to install all required R libraries properly.
        Please check the installation of all required libraries manually.\n
        Required libaries:compare"
      )
    }
    
    output$Files5 <- renderText({
      "<pre>2 Files have been created:
      1. TaxaCombined.tab
      2. OTUsCombined.tab</pre>"
    })
    session$sendCustomMessage('activeNavs', 'Serial-Group-Comparisons')
    })
    setwd("../")
    
    #################################################################################
    ######                           End of Script                             ######
    #################################################################################
    
    })
  
  output$grorder <- renderUI({
    setwd(wd)
    
    req(input$sgcgroup)
    
    input_filename = "5.Serial-Group-Comparisons/OTUsCombined.tab"
    original_table <-
      read.table (
        file = input_filename,
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
    
    setwd("../")
    orderInput(
      inputId = 'grorder',
      label = "Order of groups to display (Drag and drop)",
      items = levels(original_table[, input$sgcgroup])
    )
    
  })
  
  observeEvent(input$Submit6, {
    setwd(wd)
    tryCatch(error=function(e){
      
      output$Error5 <- renderUI(
        tags$pre(
          paste0("Error received:",
                 as.character(e),
                 "Something went wrong...")
          ,style="color:red")
      )
      # output$Files6 <- renderUI(NULL)
      # output$box <- renderText({NULL})
      # output$boxpoint <- renderText({NULL})
      # output$violin <- renderText({NULL})
      
      session$sendCustomMessage('deactiveNavs', 'Correlations')
      tryCatch(expr = {dev.off()}, error=function(e) NA)
    },
    expr = {
    #####################################################################################################################
    ####                                            Input handling                                                   ####
    #####################################################################################################################
    
    
    input_filename <- "5.Serial-Group-Comparisons/OTUsCombined.tab"
    
    independant_variable_name <- input$sgcgroup
    
    group_order <- input$grorder_order
    
    abundance_cutoff <- input$abundance_cutoff
    
    prevalence_cutoff <- input$prevalence_cutoff
    
    max_median_cutoff <- input$max_median_cutoff
    
    ReplaceZero <- ifelse(input$repzero == T, "YES", "NO")
    
    PlotOption <- input$plotOption
    
    sig.cutoff <- input$sig.cutoff
    
    output$box <- renderText({
      return(paste('', sep = ""))
    })
    
    output$boxpoint <- renderText({
      return(paste('', sep = ""))
    })
    
    output$violin <- renderText({
      return(paste('', sep = ""))
    })
    
    
    #####################################################################################################################
    ####                                        Functions to be used  in main Script.                            ########
    #####################################################################################################################
    
    # Function to calculate relative abundance
    rel.abundance <- function(data)
    {
      total = sum(data)
      rel.data <- 100 * data / total
      return(rel.data)
    }
    
    # Replace abundance with zero if value is below given cutoff
    abundance.fix <- function(data)
    {
      data[data < abundance_cutoff] <- 0
      return(data)
    }
    
    # Replace zero value with NA
    fill_zero.NA <- function(data, ReplaceZero)
    {
      if (ReplaceZero == "NO") {
        return(data)
      } else if (ReplaceZero == "YES") {
        data[data == 0] <- NA
        return(data)
      } else {
        return(data)
      }
    }
    
    # Return maxima and minima values of the given input value
    max.fun <- function(data)
    {
      data.max <- max(as.numeric(as.character(data)), na.rm = TRUE)
      return (data.max)
    }
    
    # Calculate the prevalence of the given input table
    pre.fun.na <- function(data)
    {
      prevalence <- nnzero(data, na.counted = FALSE)
      return(prevalence)
    }
    
    # Return the maximum median value for each group
    max.med <- function(data)
    {
      max.median <-
        max(aggregate (
          data ~ independent_variable,
          FUN = median,
          simplify = TRUE
        )[, 2])
      return(max.median)
    }
    
    # Return the prevalence ratio
    max.pre <- function(data)
    {
      # Return the number of samples (excluding NA) for each group separately
      found <-
        aggregate (
          data ~ independent_variable,
          FUN = pre.fun.na,
          simplify = TRUE,
          na.action = na.pass
        )[, 2]
      
      # Return the total number of samples (including missing values) for each group separately
      all <-
        aggregate (
          data ~ independent_variable,
          FUN = length,
          simplify = TRUE,
          na.action = na.pass
        )[, 2]
      
      # Calculate the ratio for each group and return the maximum ratio out of the groups
      max.ratio <- max(found / all)
      return(max.ratio)
    }
    
    # Set the theme to change text for plotting (ggplot - Gtable)
    mytheme <- gridExtra::ttheme_default(
      # Adjust settings for the text inside table
      core = list(fg_params = list(cex = 0.8)),
      
      # Adjust the test for column and row header
      colhead = list(fg_params = list(cex = 0.9)),
      rowhead = list(fg_params = list(cex = 1.0))
    )
    
    ###################       Read all required input files      ####################
    
    # Load the tab-delimited file containing the values to be analyzed (samples names in the first column)
    original_table <-
      read.table (
        file = input_filename,
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
    
    # Determine positions in table
    dependant_variables_start <-
      min(which(sapply(original_table, is.factor) == FALSE))
    
    taxonomic_variables_start <-
      min(which(grepl("OTU", colnames(original_table))))
    
    
    #####################################################################################################################
    ####                                      Pre-processing of OTUs Table                                       ########
    #####################################################################################################################
    
    # Convert independent variable into factor to avoid errors
    original_table[, independant_variable_name] <-
      as.factor(original_table[, independant_variable_name])
    
    # Store independent variable columns from original table
    independent_variable <-
      original_table[[independant_variable_name]]
    ifelse(
      group_order != "",
      independent_variable <-
        factor(independent_variable, levels = group_order),
      independent_variable
    )
    
    # Store metadata variable columns from original table
    my_meta_data <- original_table[1:taxonomic_variables_start - 1]
    
    # Store relative abundance values of all OTUs
    my_otu_data <-
      original_table[taxonomic_variables_start:dim(original_table)[2]]
    
    # Transform data by zeroing very low abundances (based on given abundance cutoff - abundance_cutoff)
    my_otu_mod =  as.data.frame(apply(my_otu_data, 2, abundance.fix))
    
    # Transform data by replacing all zero values with missing values
    # Column consisting of "NA" or "0" only are removed (see below)
    my_otu_mod_noz = as.data.frame(apply(my_otu_mod, 2, fill_zero.NA, ReplaceZero))
    
    # Remove column if entire OTU column contain zeros or missing values
    my_otu_mod_noz <-
      my_otu_mod_noz[, !apply(my_otu_mod_noz , 2 , function(x)
        all(is.na(x) | (x == 0))), drop = FALSE]
    
    # Transform data by removing any OTU with median relative abundance below cutoff
    t_otu_mod_noz = as.data.frame(t(my_otu_mod_noz))
    
    # Calculate median for each OTUs
    t_otu_mod_noz$max.median <- apply(t_otu_mod_noz, 1, max.med)
    
    # Select OTUs above median cutoff (med.cutoff)
    selected_max <-
      t_otu_mod_noz[t_otu_mod_noz$max.median > max_median_cutoff, ]
    
    # Remove calculated median column "max.median"
    selected_max$max.median <- NULL
    
    # Make a separate object as data frame (columns are OTUs and rows are samples)
    otu_mod_noz_max <- as.data.frame(t(selected_max))
    
    # Transpose data (columns are samples and rows are OTUs)
    t_otu_mod_noz_max = as.data.frame(t(otu_mod_noz_max))
    
    # Transform data by removing all OTUs with prevalence below the given cutoff
    t_otu_mod_noz_max$pre <- apply(t_otu_mod_noz_max, 1, max.pre)
    selected_pre <-
      t_otu_mod_noz_max[t_otu_mod_noz_max$pre > prevalence_cutoff, ]
    selected_pre$pre <- NULL
    
    # Transform and filter OTU table
    otu_mod_noz_max_pre <- as.data.frame(t(selected_pre))
    
    # Merge the metadata and OTU data in one table
    # This table will be used as input for the analysis
    input_table <- cbind(my_meta_data, otu_mod_noz_max_pre)
    
    #####################################################################################################################
    ####                              Differential Statistical Analysis                                        ########
    #####################################################################################################################
    
    
    # The number of total observations per category (independant variable) to be used in the analysis
    total <- summary(independent_variable)
    
    # Create vector with group information
    prevalence_list <- as.numeric(independent_variable)
    
    # Create vector with prevalence values
    # Iterate through all groups
    for (i in 1:nlevels(independent_variable)) {
      for (j in 1:length(prevalence_list)) {
        # Assign the prevalence value of a group to each sample
        if (as.character(independent_variable[j]) == names(total[i])) {
          prevalence_list[j] <- as.numeric(total[i])
        }
      }
    }
    
    # Create an empty dataframe for Kruskal-Wallis Rank Sum Test
    df <-
      data.frame(name = character(0),
                 pvalue = numeric(0),
                 sign = character(0))
    
    # Create an empty dataframe for Fisher's Exact Test
    Fdf <-
      data.frame(name = character(0),
                 pvalue = numeric(0),
                 sign = character(0))
    
    # Create an empty dataframe for Wilcoxon Rank Sum and Signed Rank Test
    all_pair_pval_table <-
      data.frame(
        measure = character(0),
        pair = character(0),
        Group1 = character(0),
        Group2 = character(0),
        pvalue = numeric(0),
        corrected = numeric(0)
      )
    
    # Create an empty dataframe to the results of the paired Fisher's Exact Test
    all_pair_fpval_table <-
      data.frame(
        measure = character(0),
        pair = character(0),
        Group1 = character(0),
        Group2 = character(0),
        pvalue = numeric(0),
        corrected = numeric(0)
      )
    
    ###################              Create empty lists to store information about signficant differences              ###################
    
    # Making an object for lists of boxplot
    list_box <- list()
    
    # Making an object for lists of point boxplot
    list_point <- list()
    
    # Making an object for lists of violin boxplot
    list_violin <- list()
    
    # Making a list for text table for pairwise tests
    texttable <- list()
    
    # Making a list for pairwise p-value table of Wilcoxon test
    pvaltable <- list()
    
    # Making a list for pairwise p-value table of Fisher's test
    fpvaltable <- list() #
    
    # Making a list for overall p-value table for Fisher's test
    allfpvaltable <- list()
    
    # Making a list for overall p-value table for Kruskal test
    pvaltableAll <- list()
    
    #######################               Kruskal Wallis Test for each OTU for all Groups                #######################
    
    fail <- FALSE
    
    # Start calculation with the dependant variable (e.g. Richness)
    for (i in dependant_variables_start:dim(input_table)[2])
    {
      # Take the values for all samples of the dependant variable/OTU
      my_test_vector <- input_table[, i]
      
      # Save the name of the observed variable/OTU
      my_name <- colnames(input_table)[i]
      test_table <- xtabs(my_test_vector ~ independent_variable)
      
      # Number of groups to be compared
      num_of_represented_groups <- sum(as.vector(test_table) > 0)
      
      # Test whether any group is missing
      num_of_missing_groups <-
        nlevels(independent_variable) - num_of_represented_groups
      
      # Function to return "TRUE" or "FALSE" corresponding to "not missing" and "missing"
      if (num_of_missing_groups > 0) {
        isgroupmissing <- TRUE
      } else {
        isgroupmissing <- FALSE
      }
      
      # Performs a Kruskal-Wallis rank sum test
      fit <-
        tryCatch (
          kruskal.test(my_test_vector ~ independent_variable),
          error = function(i) {
            fail <<- TRUE
          }
        )
      
      # Function to assign corrected and not-corrected pvalues
      if (fail) {
        my_pvalue <- NaN
        #my_corrected_pvalue <- 0
        fail <- FALSE
      } else {
        # Round the p-value down to four decimals
        my_pvalue <- round(fit$p.value, 4)
      }
      
      # Add the p-values to the table
      newRow <-
        data.frame(name = my_name,
                   pvalue = my_pvalue,
                   missing = isgroupmissing)
      df <- rbind(df, newRow)
    }
    
    # Applying Benjamini-Hochberg (1995) correction
    df$corrected <- round(p.adjust(df$pvalue, method = "BH"), 4)
    
    #######################   Paired-wilcoxon Test, (paired) Fisher's Exact Test    ########
    
    # Kruskal-Wallis test for analysis of variances
    # Wilcoxon signed-rank test for pairwise comparisons
    count <- 1
    x <- 0
    
    # Vector with all possible group combinations
    idx <- combn(nlevels(independent_variable), 2)
    
    for (i in dependant_variables_start:dim(input_table)[2])
    {
      flag = TRUE
      
      # The vector of a dependant variable/OTU
      my_test_vector <- input_table[, i]
      
      # The name in the header of the dependant variable/OTU
      my_name <- colnames(input_table)[i]
      
      # Save the p-value
      pvalue <- df[count, 2]
      
      # Save information about missing group
      missing_group <- df[count, 3]
      
      # Save the corrected p-value
      cpvalue <- df[count, 4]
      count <- count + 1
      signif_pairs = data.frame(
        measure = as.character(),
        name = as.character(),
        Group1 = as.character(),
        Group2 = as.character(),
        pvalue = as.numeric()
      )
      
      # Get the names of all group combinations
      idx_name <- combn(levels(independent_variable), 2)
      
      # If a pvalue is calculated out of the Kruskal Walis Test, the pairwise Wilcos Rank Test will be computed as well
      if (!is.nan(pvalue)) {
        if (pvalue <= sig.cutoff) {
          # Compute p-values from Wilcoxon test for all comparisons
          ppval_res <- numeric(ncol(idx))
          
          # Create an empty dataframe to hold the results of pairwise comparison
          pair_pval_table <-
            data.frame(
              measure = character(0),
              pair = character(0),
              Group1 = character(0),
              Group2 = character(0),
              pvalue = numeric(0),
              corrected = numeric(0)
            )
          pair_pval_table <-
            rbind(
              pair_pval_table,
              data.frame(
                measure = my_name,
                name = "All",
                Group1 = " -",
                Group2 = "- ",
                pvalue = pvalue,
                corrected = cpvalue
              )
            )
          # Compute p-values of Wilcoxon test for all comparisons
          for (i in 1:ncol(idx)) {
            # Performs a Kruskal-Wallis rank sum test
            fit <-
              tryCatch (
                wilcox.test(my_test_vector[as.numeric(independent_variable) == idx[1, i]], my_test_vector[as.numeric(independent_variable) == idx[2, i]]),
                error = function(i) {
                  fail <<- TRUE
                }
              )
            # Function to assign corrected and not-corrected pvalues
            if (fail) {
              ppval_res[i] <- NaN
              #my_corrected_pvalue <- 0
              fail <- FALSE
            } else {
              # Round the p-value down to four decimals
              ppval_res[i] <- round(fit$p.value, 4)
            }
            # Calculate p-value
            #ppval_res[i] <- wilcox.test(my_test_vector[as.numeric(independent_variable) == idx[1,i]],my_test_vector[as.numeric(independent_variable) == idx[2,i]])$p.value
            
            # Set the values of the pair and the corresponding p-value
            pair_name <-
              paste (idx_name[1, i], "-", idx_name[2, i],  sep = "")
            pair_num <- paste (idx[1, i], "-", idx[2, i],  sep = "")
            ppval <- round(ppval_res[i], 4)
            
            # Create and add a new column to the plot table and the overall pairwise comparison table
            newRow <-
              data.frame(
                measure = my_name,
                name = pair_num,
                Group1 = idx_name[1, i],
                Group2 = idx_name[2, i],
                pvalue = ppval,
                corrected = 0
              )
            pair_pval_table <- rbind(pair_pval_table, newRow)
          }
          
          # Add the corrected p-values column to the dataframe
          pair_pval_table$corrected[-1] <-
            round(p.adjust(pair_pval_table$pvalue[-1], method = "BH"), 4)
          Pforplot_table <- pair_pval_table[, c(-3, -4)]
          
          # Add the table with the corrected p-values to the complete list of pairwise p-values for all tests
          all_pair_pval_table <-
            rbind(all_pair_pval_table, pair_pval_table)
          
          # Determine which groups are significantly different based on the results of the Wilcoxon test
          signif_pairs <-
            Pforplot_table[(Pforplot_table$pvalue < sig.cutoff) &
                             !(is.na(Pforplot_table$pvalue)), ]
        }
      }
      
      # If zeros are replaced by missing values, the prevalence and the points to be plotted are the same (number of samples with a value >0)
      # If zeros are not replaced and considered as true values, the prevalence is as above, but zeros are going to be plotted in the graphs
      if (ReplaceZero == "YES") {
        plot_df <-
          cbind.data.frame(abundance = my_test_vector, variable = independent_variable)
        plot_df$samplekaname <- row.names(input_table)
        plot_df_prevalence <- plot_df
      } else {
        plot_df <-
          cbind.data.frame(abundance = my_test_vector, variable = independent_variable)
        plot_df$samplekaname <- row.names(input_table)
        my_test_vector <-
          fill_zero.NA(my_test_vector, ReplaceZero = "YES")
        plot_df_prevalence <-
          cbind.data.frame(abundance = my_test_vector, variable = independent_variable)
        plot_df_prevalence$samplekaname <- row.names(input_table)
      }
      
      
      # Calculate prevalence by counting presence or absence of the given variable in each group
      prevalence <-
        table(plot_df_prevalence[!is.na(plot_df_prevalence[, 1]), 2])
      
      # Count how many samples are absent from total number of counts
      not_found <- total - prevalence
      pre_table <- cbind(prevalence, not_found, total)
      
      # Calculate a two-sided Fisher's test
      Fishtest <-
        fisher.test(pre_table[, -3], alternative = "two.sided")
      
      # Save the p-value of the Fisher's test
      fish_pvalue <- round(Fishtest$p.value, 4)
      FnewRow <- data.frame(name = my_name, pvalue = fish_pvalue)
      Fdf <- rbind(Fdf, FnewRow)
      fppval_res <- numeric(ncol(idx))
      
      # Create an empty dataframe to hold the results of pairwise comparison for pairwise Fisher's Test
      pair_fpval_table <-
        data.frame(
          measure = character(0),
          name = character(0),
          Group1 = character(0),
          Group2 = character(0),
          pvalue = numeric(0),
          corrected = numeric(0)
        )
      
      # Compute p-values from Wilcoxon test for all comparisons
      for (i in 1:ncol(idx))
      {
        pret <- as.data.frame(pre_table)
        pretrowbind <-
          rbind(pret[idx_name[1, i], -3], pret[idx_name[2, i], -3])
        
        # Compute two-sided Fisher's test
        fppval_res[i] <-
          fisher.test(pretrowbind,
                      alternative = "two.sided",
                      workspace = 2e8)$p.value
        
        # Set values of the pair and corresponding p-value
        # Take variable name as "A-B"
        pair_name <-
          paste (idx_name[1, i], "-", idx_name[2, i],  sep = "")
        
        # Take variable name as "1-2"
        pair_num <- paste (idx[1, i], "-", idx[2, i],  sep = "")
        
        # Round the p-value to four decimals
        fppval <- round(fppval_res[i], 4)
        
        # Create and add a new line to the plot table and the overall pairwise comparison table
        newRow <-
          data.frame(
            measure = my_name,
            name = pair_num,
            Group1 = idx_name[1, i],
            Group2 = idx_name[2, i],
            pvalue = fppval
          )
        pair_fpval_table <- rbind(pair_fpval_table, newRow)
      }
      
      # Applying Benjamini-Hochberg (1995) correction
      # Add a column with the corrected p-values
      pair_fpval_table$corrected <-
        round(p.adjust(pair_fpval_table$pvalue, method = "BH"), 4)
      
      # Add the table with corrected p-values to the complete list of pairwise p-values for all tests
      all_pair_fpval_table <-
        rbind(all_pair_fpval_table, pair_fpval_table)
      
      # Make an object of "measure", "name", "p-value", "corrected" to print in PDF
      Fforplot_table <- pair_fpval_table[, c(-3, -4)]
      
      # Test whether the test is significant or not
      # If significant, then get value in "signif_fpairs"
      signif_fpairs <-
        Fforplot_table[(Fforplot_table$pvalue <= sig.cutoff) &
                         !(is.na(Fforplot_table$pvalue)), ]
      
      
      # Check if at least one of the tests (Kruskal Walis or Fisher) is signficant
      # If at least one of them is signficant a plot will be generated
      if (fish_pvalue <= sig.cutoff ||
          (pvalue <= sig.cutoff & !is.nan(pvalue))) {
        # Generate label for the X-axis
        labelsx <-
          as.data.frame(table(plot_df[!is.na(plot_df_prevalence[, 1]), 2]))
        labeling <-
          paste(labelsx$Var1, "(", labelsx$Freq, "/", total, ")", sep = "")
        plot_df$xlabeltext <-
          factor(
            paste(
              plot_df$variable,
              "(",
              labelsx[match(plot_df$variable, labelsx$Var1), "Freq"],
              "/",
              prevalence_list,
              ")",
              sep = ""
            ),
            level = labeling
          )
        x = x + 1
        
        # Determine label of the y-axis
        if ((dependant_variables_start + (count - 2)) < taxonomic_variables_start) {
          labelsy <- my_name
        } else{
          labelsy <- "Relative abundance (%)"
        }
        # Create a ggplot object of the plotted layout (including axis labels and scaling)
        g <- ggplot(plot_df, aes(x = xlabeltext, y = abundance))
        # Generate plot dependent on given PlotOption
        if (PlotOption == 1) {
          my_boxplot <-
            g + stat_boxplot(geom = "errorbar", width = 0.25) + geom_boxplot(varwidth = FALSE, width = 0.7)
          my_violinplot <-
            g + geom_violin(width = 0.7) + geom_boxplot(width = 0.1)
          my_point_boxplot <-
            g + geom_dotplot(binaxis = 'y',
                             stackdir = 'center',
                             dotsize = 0.7)
        } else if (PlotOption == 2) {
          my_boxplot <-
            g + stat_boxplot(geom = "errorbar", width = 0.25) + geom_boxplot(varwidth = FALSE, width = 0.7) + geom_dotplot(binaxis = 'y',
                                                                                                                           stackdir = 'center',
                                                                                                                           dotsize = 0.7)
          my_violinplot <-
            g + geom_violin(width = 0.7) + geom_boxplot(width = 0.1) + geom_dotplot(binaxis = 'y',
                                                                                    stackdir = 'center',
                                                                                    dotsize = 0.5)
          my_point_boxplot <-
            g + geom_dotplot(binaxis = 'y',
                             stackdir = 'center',
                             dotsize = 0.5)
        } else if (PlotOption == 3) {
          my_boxplot <-
            g + stat_boxplot(geom = "errorbar", width = 0.25) + geom_boxplot(varwidth = FALSE, width = 0.7) +
            geom_dotplot(binaxis = 'y',
                         stackdir = 'center',
                         dotsize = 0.7) +
            geom_label_repel(
              aes(
                xlabeltext,
                abundance,
                fill = NULL,
                label = samplekaname
              ),
              fontface = 'bold',
              color = 'black',
              box.padding = unit(0.25, "lines"),
              point.padding = unit(0.5, "lines"),
              size = 2
            )
          my_violinplot <-
            g + geom_violin(width = 0.7) + geom_boxplot(width = 0.1) +
            geom_dotplot(binaxis = 'y',
                         stackdir = 'center',
                         dotsize = 0.7) + geom_label_repel(
                           aes(
                             xlabeltext,
                             abundance,
                             fill = NULL,
                             label = samplekaname
                           ),
                           fontface = 'bold',
                           color = 'black',
                           box.padding = unit(0.25, "lines"),
                           point.padding = unit(0.5, "lines"),
                           size = 2
                         )
          my_point_boxplot <-
            g + geom_dotplot(binaxis = 'y',
                             stackdir = 'center',
                             dotsize = 0.7) + geom_label_repel(
                               aes(
                                 xlabeltext,
                                 abundance,
                                 fill = NULL,
                                 label = samplekaname
                               ),
                               fontface = 'bold',
                               color = 'black',
                               box.padding = unit(0.25, "lines"),
                               point.padding = unit(0.5, "lines"),
                               size = 2
                             )
        }
        
        my_boxplot <-
          my_boxplot + ggtitle(my_name) + guides(fill = FALSE) + ylab(labelsy) + xlab("") + theme_bw() + theme(
            axis.text.x = element_text(
              colour = "grey20",
              size = 12,
              angle = 45,
              hjust = 1,
              vjust = 1,
              face = "plain"
            ),
            axis.text.y = element_text(
              colour = "grey20",
              size = 12,
              angle = 0,
              hjust = 1,
              vjust = 0,
              face = "plain"
            ),
            axis.title.y = element_text(
              colour = "grey20",
              size = 14,
              angle = 90,
              hjust = .5,
              vjust = .5,
              face = "plain"
            ),
            plot.title = element_text(
              colour = "grey22",
              size = 18,
              hjust = .5,
              vjust = .5,
              face = "bold"
            )
          ) + theme(plot.margin = unit(c(1, 1, 1.5, 1.2), "cm"))
        
        my_violinplot <-
          my_violinplot + ggtitle(my_name) + ylab(labelsy) + guides(fill = FALSE) + xlab("") + theme_bw() +
          theme(
            axis.text.x = element_text(
              colour = "grey20",
              size = 12,
              angle = 45,
              hjust = 1,
              vjust = 1,
              face = "plain"
            ),
            axis.text.y = element_text(
              colour = "grey20",
              size = 12,
              angle = 0,
              hjust = 1,
              vjust = 0,
              face = "plain"
            ),
            axis.title.y = element_text(
              colour = "grey20",
              size = 14,
              angle = 90,
              hjust = .5,
              vjust = .5,
              face = "plain"
            ),
            plot.title = element_text(
              colour = "grey22",
              size = 18,
              hjust = .5,
              vjust = .5,
              face = "bold"
            )
          ) + theme(plot.margin = unit(c(1, 1, 1.5, 1.2), "cm"))
        
        my_point_boxplot <-
          my_point_boxplot + ggtitle(my_name) + theme_bw() + ylab(labelsy) + stat_summary(
            fun.y = median,
            fun.ymin = median,
            fun.ymax = median,
            geom = "crossbar",
            width = 0.5
          ) +
          guides(fill = FALSE) + xlab("") + theme(legend.position = "none") + guides(colour = FALSE) + theme(
            axis.text.x = element_text(
              colour = "grey20",
              size = 12,
              angle = 45,
              hjust = 1,
              vjust = 1,
              face = "plain"
            ),
            axis.text.y = element_text(
              colour = "grey20",
              size = 12,
              angle = 0,
              hjust = 1,
              vjust = 0,
              face = "plain"
            ),
            axis.title.y = element_text(
              colour = "grey20",
              size = 14,
              angle = 90,
              hjust = .5,
              vjust = .5,
              face = "plain"
            ),
            plot.title = element_text(
              colour = "grey22",
              size = 18,
              hjust = .5,
              vjust = .5,
              face = "bold"
            )
          ) + theme(plot.margin = unit(c(1, 1, 1.5, 1.2), "cm"))
        
        # Save the boxplot object in the list
        list_box[[x]] <- list()
        list_box[[x]] <- my_boxplot
        
        # Save the boxplot with points object in the list
        list_point[[x]] <- list()
        list_point[[x]] <- my_point_boxplot
        
        # Save the violin object in the list
        list_violin[[x]] <- list()
        list_violin[[x]] <- my_violinplot
        
        
        
        # Make a list object to store all "prevalence table" to print in PDF
        texttable[[x]] <- list()
        
        # Creat a gtable containing text grobs representing a character matrix
        # "mytheme" is a text theme defined at the begining
        texttable[[x]] <- tableGrob(pre_table, theme = mytheme)
        
        # Define the title of the tables
        title <-
          textGrob("Prevalence table", gp = gpar(fontsize = 7))
        padding <- unit(3, "mm")
        texttable[[x]] <-
          gtable_add_rows(texttable[[x]],
                          heights = grobHeight(title) + padding,
                          pos = 0)
        texttable[[x]] <-
          gtable_add_grob(texttable[[x]], title, 1, 1, 1, ncol(texttable[[x]]))
        
        
        # Test for the case that there are no significant pairs
        # Check whether Fisher test is significant or not
        if (fish_pvalue <= sig.cutoff) {
          if (!is.na(signif_fpairs[1, 1])) {
            signif_fpairs$measure <- as.character(signif_fpairs$measure)
            signif_fpairs$name <- as.character(signif_fpairs$name)
            colnames(signif_fpairs) <-
              c("Species", "Groups", "p-value", "Adj. p-value")
            fpvaltable[[x]] <- list()
            
            # Pavlue table for significant pairs
            fpvaltable[[x]] <-
              tableGrob(signif_fpairs[2:4],
                        rows = NULL,
                        theme = mytheme)
            title <-
              textGrob("Fisher's Exact Test - pairwise", gp = gpar(fontsize = 9))
            padding <- unit(3, "mm")
            fpvaltable[[x]] <-
              gtable_add_rows(fpvaltable[[x]],
                              heights = grobHeight(title) + padding,
                              pos = 0)
            fpvaltable[[x]] <-
              gtable_add_grob(fpvaltable[[x]], title, 1, 1, 1, ncol(fpvaltable[[x]]))
          } else {
            FnewRow$name <- as.character(FnewRow$name)
            FnewRow$pvalue <- as.character(FnewRow$pvalue)
            FnewRow <- cbind(FnewRow$name, "-", FnewRow$pvalue, 0)
            colnames(FnewRow) <-
              c("Species", "Groups", "p-value", "Adj. p-value")
            allfpvaltable[[x]] <- list()
            allfpvaltable[[x]] <-
              tableGrob(FnewRow, rows = NULL, theme = mytheme)
            #fpvaltable[[x]] <- NULL
          }
        }
        # Check whether Kruskal-Wallis test is significant or not
        if (!is.na(signif_pairs[1, 1])) {
          signif_pairs$measure <- as.character(signif_pairs$measure)
          signif_pairs$name <- as.character(signif_pairs$name)
          colnames(signif_pairs) <-
            c("Species", "Groups", "p-value", "Adj. p-value")
          pvaltable[[x]] <- list()
          pvaltableAll[[x]] <- list()
          
          # Pvalue table for significant pairs
          signif_all <- signif_pairs[1, c(1, 3, 4)]
          #if(dim(signif_pairs)[1] > 1) {
          signif_pairs <- signif_pairs[2:dim(signif_pairs)[1], ]
          pvaltable[[x]] <-
            tableGrob(signif_pairs[2:4], rows = NULL, theme = mytheme)
          
          # Title of tables in the PDF
          title <-
            textGrob("Wilcoxon Rank Sum Test - pairwise", gp = gpar(fontsize = 9))
          padding <- unit(2, "mm")
          pvaltable[[x]] <-
            gtable_add_rows(pvaltable[[x]],
                            heights = grobHeight(title) + padding,
                            pos = 0)
          pvaltable[[x]] <-
            gtable_add_grob(pvaltable[[x]], title, 1, 1, 1, ncol(pvaltable[[x]]))
          colnames(signif_all) <- c(" ", "p-value", "Adj. p-value")
          pvaltableAll[[x]] <-
            tableGrob(signif_all, rows = NULL, theme = mytheme)
          title <-
            textGrob("Kruskal-Wallis Rank Sum Test - all groups ",
                     gp = gpar(fontsize = 9))
          padding <- unit(2, "mm")
          pvaltableAll[[x]] <-
            gtable_add_rows(pvaltableAll[[x]],
                            heights = grobHeight(title) + padding,
                            pos = 0)
          pvaltableAll[[x]] <-
            gtable_add_grob(pvaltableAll[[x]], title, 1, 1, 1, ncol(pvaltableAll[[x]]))
        }
      }
      
    }
    
    
    # Apply Benjamini-Hochberg (1995) correction
    # Add a column with the corrected p-values
    Fdf$corrected <- round(p.adjust(Fdf$pvalue, method = "BH"), 4)
    sig_Fdf <- subset(Fdf, Fdf$pvalue <= sig.cutoff)
    counter = 1
    
    if (length(allfpvaltable) != length(fpvaltable)) {
      if (length(allfpvaltable) > length(fpvaltable)) {
        for (i in (length(fpvaltable) + 1):(length(allfpvaltable) + 1)) {
          fpvaltable[[i]] <- list()
          fpvaltable[[i]] <- NULL
        }
      }
      else{
        for (i in (length(allfpvaltable) + 1):(length(fpvaltable) + 1)) {
          allfpvaltable[[i]] <- list()
          allfpvaltable[[i]] <- NULL
        }
      }
    }
    
    if (length(fpvaltable) != 0) {
      for (i in 1:length(fpvaltable)) {
        if (!is.null(fpvaltable[[i]]) || !is.null(allfpvaltable[[i]])) {
          colnames(sig_Fdf) <- c(" ", "p-value", "Adj. p-value")
          allfpvaltable[[i]] <- list()
          title_all <-
            textGrob("Fisher's Exact Test - all groups", gp = gpar(fontsize = 9))
          padding <- unit(3, "mm")
          allfpvaltable[[i]] <-
            tableGrob(sig_Fdf[counter, ], rows = NULL, theme = mytheme)
          allfpvaltable[[i]] <-
            gtable_add_rows(allfpvaltable[[i]],
                            heights = grobHeight(title_all) + padding,
                            pos = 0)
          allfpvaltable[[i]] <-
            gtable_add_grob(allfpvaltable[[i]],
                            title_all,
                            1,
                            1,
                            1,
                            ncol(allfpvaltable[[i]]))
          counter = counter + 1
        }
      }
    }
    
    # If the last pvalues are only signficant for fisher OR Wilcox additional empty lists has to be added to the array of lists
    if (length(pvaltableAll) < length(allfpvaltable)) {
      # Add empty list entries in the Wilcox and Kruskal Walis tables
      for (i in (length(pvaltableAll) + 1):(length(allfpvaltable) + 1)) {
        pvaltable[[i]] <- list()
        pvaltable[[i]] <- NULL
        pvaltableAll[[i]] <- list()
        pvaltableAll[[i]] <- NULL
      }
    }
    
    if (length(pvaltableAll) > length(allfpvaltable)) {
      # Add empty list entries in the Fisher tables
      for (i in (length(allfpvaltable) + 1):(length(pvaltableAll) + 1)) {
        fpvaltable[[i]] <- list()
        fpvaltable[[i]] <- NULL
        allfpvaltable[[i]] <- list()
        allfpvaltable[[i]] <- NULL
      }
    }
    
    ####################                 Print plots in PDF               ############################
    
    # Take current path in one variable to store results in seperate folders in further steps
    OriginalPath <- getwd()
    
    # Take the name of the independent variable to name the folder
    prefix = paste(independant_variable_name, "OTUsCombined", sep = "_")
    
    # Make a directory name with independent variable name and date
    newdir <- paste(prefix, Sys.Date(), sep = "_")
    
    # Create a directory
    if (!dir.exists(paste("5.Serial-Group-Comparisons/", newdir, sep = ""))){
      dir.create(paste("5.Serial-Group-Comparisons/", newdir, sep = ""))
    }
    else {
      system(paste("rm -Rf", paste("5.Serial-Group-Comparisons/", newdir, sep = "")))
      dir.create(paste("5.Serial-Group-Comparisons/", newdir, sep = ""))
    }
    
    # Set path for all outputs to the new directory
    setwd(paste("5.Serial-Group-Comparisons/", newdir, sep = ""))
    if (!dim(all_pair_pval_table)[1] == 0) {
      # Passing type of plots
      for (plotType in c("box", "boxpoint", "violin"))
      {
        # Switching to plots type
        plotb <-
          switch(
            plotType,
            box = list_box,
            boxpoint = list_point,
            violin = list_violin
          )
        
        # Open a PDF to print all outputs
        pdf(paste(input$sgcgroup, "plot", plotType, ".pdf", sep = ""),
            onefile = T)
        for (pos in 1:x) {
          # No significant results for Kruskal Walis Test
          if (is.null(pvaltableAll[[pos]])) {
            # Fisher Test is not signficant
            if (is.null(allfpvaltable[[pos]])) {
              plotb[[pos]]
            }
            # Signfiicant Fisher Test
            else {
              # Check if paired Fisher is signficant as well
              if (is.null(fpvaltable[[pos]])) {
                # Only the overall Fisher Test is significant
                grid.arrange(
                  plotb[[pos]],
                  arrangeGrob(
                    allfpvaltable[[pos]],
                    nrow = 1,
                    ncol = 1
                  ),
                  nrow = 2,
                  ncol = 2,
                  heights = c(3, 1),
                  as.table = T
                )
              }
              else {
                # Both Fisher Tests are significant
                if (nlevels(independent_variable) > 5) {
                  print(plotb[[pos]])
                  print(
                    plot_grid(
                      allfpvaltable[[pos]],
                      fpvaltable[[pos]],
                      ncol = 1,
                      nrow = 2,
                      align = "v",
                      rel_heights = c(1, 10)
                    )
                  )
                }
                else{
                  grid.arrange(
                    plotb[[pos]],
                    arrangeGrob(
                      allfpvaltable[[pos]],
                      fpvaltable[[pos]],
                      nrow = 2,
                      ncol = 1
                    ),
                    nrow = 2,
                    ncol = 2,
                    heights = c(3, 1),
                    as.table = T
                  )
                }
              }
            }
          }
          # Significant Kruskal Walis Test
          else {
            # Fisher test is not signficant
            if (is.null(allfpvaltable[[pos]])) {
              # Add table only for Kruskal Walis Test
              if (is.null(pvaltable[[pos]])) {
                grid.arrange(
                  plotb[[pos]],
                  arrangeGrob(
                    pvaltableAll[[pos]],
                    nrow = 1,
                    ncol = 1
                  ),
                  nrow = 2,
                  ncol = 2,
                  heights = c(3, 1),
                  as.table = T
                )
              }
              # Add table for Kruskal Walis and Wilcoxon test
              else{
                if (nlevels(independent_variable) > 5) {
                  print(plotb[[pos]])
                  print(
                    plot_grid(
                      pvaltableAll[[pos]],
                      pvaltable[[pos]],
                      ncol = 1,
                      nrow = 2,
                      align = "v",
                      rel_heights = c(1, 10)
                    )
                  )
                }
                else{
                  grid.arrange(
                    plotb[[pos]],
                    arrangeGrob(
                      pvaltableAll[[pos]],
                      pvaltable[[pos]],
                      nrow = 2,
                      ncol = 1
                    ),
                    nrow = 2,
                    ncol = 2,
                    heights = c(3, 1),
                    as.table = T
                  )
                }
              }
            }
            # Significant Kruskal Walis and Fisher Test
            else {
              # The paired Fisher test is not significant
              if (is.null(fpvaltable[[pos]])) {
                if (nlevels(independent_variable) > 5) {
                  print(plotb[[pos]])
                  print(
                    plot_grid(
                      allfpvaltable[[pos]],
                      pvaltableAll[[pos]],
                      pvaltable[[pos]],
                      ncol = 2,
                      nrow = 2,
                      align = "v",
                      rel_heights = c(1, 10)
                    )
                  )
                } else{
                  # Add table for both tests
                  grid.arrange(
                    plotb[[pos]],
                    arrangeGrob(
                      pvaltableAll[[pos]],
                      pvaltable[[pos]],
                      allfpvaltable[[pos]],
                      nrow = 4,
                      ncol = 1
                    ),
                    nrow = 2,
                    ncol = 2,
                    heights = c(3, 1),
                    as.table = T
                  )
                }
              }
              else {
                # The paired Fisher Test is significant
                if (nlevels(independent_variable) > 5) {
                  print(plotb[[pos]])
                  print(
                    plot_grid(
                      allfpvaltable[[pos]],
                      pvaltableAll[[pos]],
                      fpvaltable[[pos]],
                      pvaltable[[pos]],
                      ncol = 2,
                      nrow = 2,
                      align = "v",
                      rel_heights = c(1, 10)
                    )
                  )
                } else{
                  # Add table for both tests
                  grid.arrange(
                    plotb[[pos]],
                    arrangeGrob(
                      pvaltableAll[[pos]],
                      pvaltable[[pos]],
                      allfpvaltable[[pos]],
                      fpvaltable[[pos]],
                      nrow = 4,
                      ncol = 1
                    ),
                    nrow = 2,
                    ncol = 2,
                    heights = c(3, 1),
                    as.table = T
                  )
                }
              }
            }
          }
        }
        dev.off()
        file.copy(
          paste(input$sgcgroup, "plot", plotType, ".pdf", sep = ""),
          paste(
            "../../../www/",
            wd,
            "/",
            input$sgcgroup,
            "plot",
            plotType,
            ".pdf",
            sep = ""
          ),
          overwrite = T
        )
      }
    }
    ################################################################################
    ######                        Write Output Files in Seperate folder       ######
    #################################################################################
    
    # Main file generated after all pre-processing and used for statistical analysis
    # Serial-Group-Comparisons input Table
    write.table(
      input_table,
      paste(
        "OTUsCombined",
        "-",
        independant_variable_name,
        "-modified.txt",
        sep = ""
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Table with the results of Kruskal-Wallis test
    write.table(
      df[, c(1, 2, 4)],
      paste(
        "OTUsCombined",
        "-",
        independant_variable_name,
        "-pvalues.tab",
        sep = ""
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Table with the results of Wilcoxon Rank Sum Test
    write.table(
      all_pair_pval_table[, -2],
      paste(
        "OTUsCombined",
        "-",
        independant_variable_name,
        "-sign_pairs.tab",
        sep = ""
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Table with the results of Fisher's Exact Test
    write.table(
      Fdf,
      paste(
        "OTUsCombined",
        "-",
        independant_variable_name,
        "-FisherTestAll.tab",
        sep = ""
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Table with the results of pairwise Fisher's Exact Test
    write.table(
      all_pair_fpval_table[, -2],
      paste(
        "OTUsCombined",
        "-",
        independant_variable_name,
        "-FisherTestPairWise.tab",
        sep = ""
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Input file for correlation script (6.Correlation)
    Corr_input_table <-
      input_table[, -c(1:(dependant_variables_start - 1))]
    suppressWarnings (try(write.table(
      Corr_input_table,
      paste(
        "../../6.Correlations/",
        "OTUsCombined",
        "_",
        "Corr_input_table.tab",
        sep = ""
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    ),
    silent = TRUE)
    )
    
    # Adding log file in analysis
    sink(file = "my_analysis_log.txt")
    cat ("***************************************", "\n")
    cat ("Parameters Used for Analysis", "\n")
    cat ("***************************************", "\n", "\n")
    cat ("study.name:OTUsCombined", "\n", "\n")
    cat (
      "independant_variable_name:",
      independant_variable_name,
      "\n",
      "\n"
    )
    cat (
      "dependant_variables_start:",
      dependant_variables_start,
      "\n",
      "\n"
    )
    cat (
      "taxonomic_variables_start:",
      taxonomic_variables_start,
      "\n",
      "\n"
    )
    cat ("abundance_cutoff:", abundance_cutoff, "\n", "\n")
    cat ("prevalence_cutoff:", prevalence_cutoff, "\n", "\n")
    cat ("max_median_cutoff:", max_median_cutoff, "\n", "\n")
    cat ("PlotOption:", PlotOption, "\n", "\n")
    cat ("ReplaceZero:", ReplaceZero, "\n", "\n")
    cat ("sig.cutoff:", sig.cutoff, "\n", "\n")
    sink()
    setwd(OriginalPath)
    if (!flag) {
      stop(
        "
        It was not possible to install all required R libraries properly.
        Please check the installation of all required libraries manually.\n
        Required libaries:plotrix,PerformanceAnalytics,reshape,ggplot2,gridExtra,grid,ggrepel,gtable,Matrix,cowplot"
      )
    }
    
    ####### Render pdfs that were created #########
    output$box <- renderUI({
      tags$iframe(
        src = paste0(wd, "/", input$sgcgroup, "plotbox.pdf"),
        height = "600",
        width = "100%"
      )
    })
    
    output$boxpoint <- renderUI({
      tags$iframe(
        src = paste0(wd, "/", input$sgcgroup, "plotboxpoint.pdf"),
        height = "600",
        width = "100%"
      )
    })
    
    output$violin <- renderUI({
      tags$iframe(
        src = paste0(wd, "/", input$sgcgroup, "plotviolin.pdf"),
        height = "600",
        width = "100%"
      )
    })
    
    output$Files6 <- renderText({
      "<pre>Files created:
      Data as tab-delimited table modified from the input table according to selected abundance and prevalence thresholds.
      Output tables containing taxa names, p-values, and corrected p-values (for Wilcoxon and Kruskal-Wallis Test separately).
      Box-, dot-, or violin-plot (as PDF) with corresponding p-value tables.</pre>"
    })
    
    session$sendCustomMessage('activeNavs', 'Correlations')
  })
    
    setwd("../")
    
    ########################################################
    ##    Script Ended !!!!
    #########################################################
    
    })
  
  observeEvent(input$Submit8, {
    setwd(wd)
    tryCatch(error=function(e){
      if (as.character(e) == "Error received:Error in 1:otu_variables_start: result would be too long a vector") {
        errm="Most probable cause of error: Serial-Group-Comparisons step didn't return enough data. Please be less strict and try again"
      }
      else {
        errm=""
      }
      output$Error6 <- renderUI(
        tags$pre(
          
          
          paste("Error received:",
                 as.character(e),
                 "Something went wrong...",errm,sep='\n')
          ,style="color:red")
      )
      output$corrplot <- renderPlot(NULL)
      output$Finalpdf <- renderPlot(NULL)
      tryCatch(expr = {dev.off()}, error=function(e) NA)
    },
    expr = {
    ##### Input Handling #####
    
    input_file <- "6.Correlations/OTUsCombined_Corr_input_table.tab"
    
    signf_cutoff <- input$signf_cutoff
    
    includeTax <- input$includeTax
    
    includeMeta <- input$includeMeta
    
    fill_NA <- input$fill_NA
    
    replace_zeros <- input$replace_zeros
    
    prevalence_exclusion <- input$prevalence_exclusion
    
    min_pair_support <- input$min_pair_support
    
    plot_pval_cutoff <- input$plot_pval_cutoff
    
    plot_corr_cutoff <- input$plot_corr_cutoff
    
    
    ###################            Read input table              ####################
    # Load the tab-delimited file containing the values to be checked (rownames in the first column)
    my_data <-
      read.table (
        file = input_file,
        check.names = FALSE,
        header = TRUE,
        dec = ".",
        sep = "\t",
        row.names = 1,
        comment.char = ""
      )
    print(my_data)
    
    # Clean table from empty lines
    my_data <- my_data[!apply(is.na(my_data) | my_data == "", 1, all), ]
    
    otu_variables_start <-
      min(which(grepl("OTU", colnames(my_data))))
    
    ####################            Functions                  #####################
    
    # Function for filling missing values with the mean of the column
    fill_NA.mean <- function(vec)
    {
      # Calculate mean value of each column (excluding missing values)
      m <- mean(vec, na.rm = TRUE)
      # Replace missing values with mean
      vec[is.na(vec)] <- m
      # Return the new input data frame
      return(vec)
    }
    
    # Function to logarithmically normalized OTU values
    log_ratio <- function(data)
    {
      # Compute the logarithmus
      log_data <- log(data)
      # Calculate exponential function of column-wise mean values for finite log transformed data
      gm <- exp(mean(log_data[is.finite(log_data)]))
      # Compute the logarithmus
      log_gm <- log(gm)
      # Take the difference of both log-transformed datasets
      data <- log_data - log_gm
      # Return the new OTU table
      return(data)
    }
    
    ##################### END of FUNCTIONS ###########################
    
    ######################  MAIN PROGRAM #####################
    my_data <- as.data.frame(apply(my_data, 2, as.numeric))
    
    first_OTU <- colnames(my_data)[otu_variables_start]
    
    # Split the meta and taxonomic parts of the table
    # Choose the continuous scaled variables
    my_meta_data <<- my_data[1:otu_variables_start - 1]
    
    # Choose the taxonomic variables
    my_otu_data <- my_data[otu_variables_start:dim(my_data)[2]]
    
    # Process the meta measurements according to selection
    if (fill_NA == 0) {
      # Do not do anything, just rename the file
      my_meta_fixed <-  my_meta_data
    }
    
    if (fill_NA == 1) {
      # Fill non-zero missing meta-values with the mean of the column (optional)
      # Apply the previously implemented function 'fill_Na.mean' to the meta-data subset
      my_meta_fixed <- apply(my_meta_data, 2, fill_NA.mean)
      my_meta_fixed <- as.data.frame(my_meta_fixed)
      
    }
    
    # The maximal number of absence
    prevalence_cutoff <-
      dim(my_otu_data)[1] - (prevalence_exclusion * dim(my_otu_data)[1])
    
    # Count how many missing values are found for each OTU
    na_count <-
      sapply(my_otu_data, function(y)
        sum(length(which(is.na(
          y
        )))))
    
    # Count how many zeros are found for each OTU
    zero_count <-
      sapply(my_otu_data, function(y)
        sum(length(which(y == 0))))
    prevalence_count <- na_count + zero_count
    
    # A new OTU-table is generated, where the number of missing values is below the set cutoff
    my_otu_data <-
      my_otu_data[, prevalence_count <= prevalence_cutoff]
    
    # If the parameter is set, zeros are replaced with missing values
    if (replace_zeros == 1) {
      my_otu_data[my_otu_data == 0] <- NA
    }
    
    # Replace zeros with 0.0001 to avoid infinite number when calculating logarithmus (log(0)=-INF)
    my_otu_data[my_otu_data == 0] <- 0.0001
    
    # Transform compositional data by log ratio transformation
    my_otu_fixed = apply(my_otu_data, 2, log_ratio)
    
    # Merge the meta- and OTU data in one table
    transformed_data <- cbind(my_meta_fixed, my_otu_fixed)
    
    # Centre and scale the values
    my_scaled_data <-
      scale(transformed_data, center = TRUE, scale = TRUE)
    
    # Calculate all pairwise correlations using Pearson correlation method
    my_rcorr <- rcorr(as.matrix(my_scaled_data, type = "pearson"))
    
    # Generate vector with variable names
    var_names <- row.names(my_rcorr$r)
    
    # Depending on which parameters were set at the beginning, one query type is selected
    # In each query type, three matrices are generated: p-value matrix, correlation matrix, support matrix
    # All possibles pairs are saved in a vector (pairs)
    if (includeTax == 1 & includeMeta == 0) {
      # Correlation among OTUs and NO correlation among meta-variables
      row_names <-
        var_names[c(otu_variables_start:dim(my_rcorr$r)[1])]
      col_names <- var_names
      pairs <- expand.grid(row_names, col_names)
      my_cor_matrix <-
        my_rcorr$r[c(otu_variables_start:dim(my_rcorr$r)[1]), ]
      my_pvl_matrix <-
        my_rcorr$P[c(otu_variables_start:dim(my_rcorr$P)[1]), ]
      my_num_matrix <-
        my_rcorr$n[c(otu_variables_start:dim(my_rcorr$n)[1]), ]
      
      # Set variable for plotting
      diagonale = 0
      
    } else if (includeTax == 1 & includeMeta == 1) {
      # Correlation among OTUs and correlation among meta-variables
      row_names <- var_names
      col_names <- var_names
      pairs <- expand.grid(row_names, col_names)
      my_cor_matrix <- my_rcorr$r
      my_pvl_matrix <- my_rcorr$P
      my_num_matrix <- my_rcorr$n
      # Set variable for plotting
      diagonale = 0
      
    } else if (includeTax == 0 & includeMeta == 1) {
      # NO correlation among OTUs and correlation among meta-variables
      
      row_names <- var_names[c(1:(otu_variables_start - 1))]
      col_names <- var_names
      pairs <- expand.grid(row_names, col_names)
      my_cor_matrix <- my_rcorr$r[c(1:(otu_variables_start - 1)), ]
      my_pvl_matrix <- my_rcorr$P[c(1:(otu_variables_start - 1)), ]
      my_num_matrix <- my_rcorr$n[c(1:(otu_variables_start - 1)), ]
      # Set variable for plotting
      diagonale = 1
      
    } else {
      # NO correlation among OTUs and NO correlation among meta-variables
      
      row_names <- var_names[c(1:(otu_variables_start - 1))]
      col_names <- var_names[otu_variables_start:dim(my_rcorr$r)[1]]
      pairs <- expand.grid(row_names, col_names)
      my_cor_matrix <-
        my_rcorr$r[c(1:(otu_variables_start - 1)), c(otu_variables_start:dim(my_rcorr$r)[1])]
      my_pvl_matrix <-
        my_rcorr$P[c(1:(otu_variables_start - 1)), c(otu_variables_start:dim(my_rcorr$P)[1])]
      my_num_matrix <-
        my_rcorr$n[c(1:(otu_variables_start - 1)), c(otu_variables_start:dim(my_rcorr$n)[1])]
      # Set variable for plotting
      diagonale = 1
      
    }
    
    # Select the corresponding p-value for each pair
    p_vector <- as.vector(my_pvl_matrix)
    
    # Select the corresponding correlation coefficient for each pair
    c_vector <- as.vector(my_cor_matrix)
    
    # Select the corresponding number of observations for each pair
    n_vector <- as.vector(my_num_matrix)
    
    # Generate matrix with the pairwise comparisons
    my_pairs <-
      matrix(ncol = 5,
             c(
               as.character(pairs[, 2]),
               as.character(pairs[, 1]),
               c_vector,
               p_vector,
               n_vector
             ))
    
    # Delete all pairs with insufficient number of pairs
    my_pairs <-
      subset(my_pairs, as.numeric(my_pairs[, 5]) > min_pair_support)
    
    # Adjust p-value for multiple testing using the Benjamin-Hochberg method
    pVal_BH <- round(p.adjust(my_pairs[, 4], method = "BH"), 4)
    
    # Add the corrected p-value in the table
    my_pairs <- cbind(my_pairs, as.numeric(pVal_BH))
    
    # Remove similar pairs (values along the diagonal)
    my_pairs <-
      my_pairs[!as.character(my_pairs[, 1]) == as.character(my_pairs[, 2]), ]
    
    # Remove duplicate pairs
    my_pairs <- my_pairs[!duplicated(my_pairs[, 3]),]
    
    # Created matrix columns represent correlation coefficients, p-values, number of observations, and corrected p-values
    # Rows represent the pairs
    matrix_names <- list(
      c(rep("", times = dim(my_pairs)[1])),
      c(
        "variable1",
        "variable2",
        "correlation",
        "pValue",
        "support",
        "Corrected"
      )
    )
    
    dimnames(my_pairs) <- matrix_names
    
    # Create subset of pairs with significant p-values
    my_pairs_cutoff <-
      my_pairs[as.numeric(my_pairs[, 4]) <= signf_cutoff,]
    
    # Convert to matrix
    my_pairs_cutoff <-
      matrix(my_pairs_cutoff,
             ncol = 6,
             dimnames = list(
               c(rep("", times = dim(
                 my_pairs_cutoff
               )[1])),
               c(
                 "variable1",
                 "variable2",
                 "correlation",
                 "pvalue",
                 "support",
                 "corrected pvalue"
               )
             ))
    
    # Create subset of significant pairs with strong correlation (above 0.5)
    my_pairs_cutoff_corr <-
      my_pairs_cutoff[abs(as.numeric(my_pairs_cutoff[, 3])) >= 0.5,]
    
    # Remove columns containing no information
    my_cor_matrix <-
      my_cor_matrix[, colSums(is.na(my_cor_matrix)) != nrow(my_cor_matrix)]
    
    # Missing values in the correlation matrix are set to zero
    my_cor_matrix[is.na(my_cor_matrix)] <- 0
    #################################################################################
    ######                        Generate Graphs                              ######
    #################################################################################
    
    # Take the name of the inputfile to name the folder
    prefix = paste("OTUsCombined_Corr_input_table", sep = "_")
    
    # Make a directory name with inputfile name and date
    newdir <- paste(prefix, Sys.Date(), sep = "_")
    
    # Create a directory
    dir.create(paste0("6.Correlations/", newdir))
    
    # Save visualized correlation matrix in "corrplot.pdf"
    pdf(file = paste0("6.Correlations/", newdir, "/corrplot.pdf"))
    
    # Visualization of all correlations between meta-variables and OTUs
    corrplot(
      matrix(
        data = na.omit(my_cor_matrix),
        nrow = dim(as.data.frame(row_names))[1],
        ncol = dim(as.data.frame(col_names))[1],
        dimnames = list(row_names, col_names)
      ),
      tl.col = "black",
      tl.srt = 65,
      tl.cex = 0.6,
      cl.cex = 0.5,
      diag = diagonale
    )
    
    dev.off()
    
    output$corrplot <- renderPlot({
      corrplot(
        matrix(
          data = na.omit(my_cor_matrix),
          nrow = dim(as.data.frame(row_names))[1],
          ncol = dim(as.data.frame(col_names))[1],
          dimnames = list(row_names, col_names)
        ),
        tl.col = "black",
        tl.srt = 65,
        tl.cex = 0.6,
        cl.cex = 0.5,
        diag = diagonale
      )
    })
    
    # Check if significance value for graphical output was modified by the user
    if (plot_pval_cutoff != signf_cutoff |
        plot_corr_cutoff != 0.5) {
      # Generate a new matrix with the signficance cutoff
      my_pairs_cutoff <-
        my_pairs[as.numeric(my_pairs[, 4]) <= plot_pval_cutoff,]
      
      # Extract all significant pairs with the set correlation cutoff
      corr_pval_cutoff <-
        my_pairs_cutoff[abs(as.numeric(my_pairs_cutoff[, 3])) >= plot_corr_cutoff,]
      corr_pval_cutoff <-
        matrix(corr_pval_cutoff,
               ncol = 6,
               dimnames = list(
                 c(rep("", times = dim(
                   corr_pval_cutoff
                 )[1])),
                 c(
                   "variable1",
                   "variable2",
                   "correlation",
                   "pvalue",
                   "support",
                   "corrected pvalue"
                 )
               ))
    } else  {
      # If the significance cutoff is 0.05 an the correlation cutoff is 0.5 (for plotting)
      # Take the previously generated matrix
      corr_pval_cutoff <- my_pairs_cutoff_corr
      corr_pval_cutoff <-
        matrix(corr_pval_cutoff,
               ncol = 6,
               dimnames = list(
                 c(rep("", times = dim(
                   corr_pval_cutoff
                 )[1])),
                 c(
                   "variable1",
                   "variable2",
                   "correlation",
                   "pvalue",
                   "support",
                   "corrected pvalue"
                 )
               ))
    }
    
    # Save linearized transformed correlations of significant pairs in "linear_sign_pairs.pdf"
    pdf(
      file = paste0("6.Correlations/", newdir, "/linear_sign_pairs.pdf"),
      family = "sans",
      fillOddEven = TRUE
    )
    
    # Iterate through all significant pairs
    for (i in 1:dim(corr_pval_cutoff)[1]) {
      # Save log-scaled transformed values of the first variable of the pair
      tryCatch({
        x_df <-
          transformed_data[names(transformed_data) %in% corr_pval_cutoff[i, 2]]
      },
      error = function(e) {
        print(e)
        return(NULL)
      })
      
      # Save as numerical vector
      x <- x_df[, 1]
      
      # Save log-scaled transformed values of the second variable of the pair
      y_df <-
        transformed_data[names(transformed_data) %in% corr_pval_cutoff[i, 1]]
      
      # Save as numerical vector
      y <- y_df[, 1]
      
      # Create a linear model for the pair (excluding missing values)
      clm <- lm(y ~ x, na.action = na.exclude)
      
      # Determine the number of steps for the generation of the confidence intervals
      steps <- (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) / 1000
      
      # Sequence of more finely/evenly spaced data than original one for the calculation of the confidence interval
      newx <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), steps)
      
      # Predict confidence interval based on the generated linear model
      # Prediction intervals are calculated based on the residuals of the regression equation
      # Prediction intervals account for the variability around the mean response inherent in any prediction
      # It represents the range where a single new observation is likely to fall
      a <-
        predict(clm, newdata = data.frame(x = newx), interval = "confidence")
      
      # Graphic display of all log-scaled transformed values of the ith-pair
      plot(
        x,
        y,
        ylab = "",
        xlab = "",
        cex.axis = 0.75 ,
        xaxt = 'n',
        yaxt = 'n',
        xlim = c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
        xaxs = "i"
      )
      title(
        xlab = names(x_df),
        line = 0.5,
        font.lab = 2,
        cex.lab = 1.4
      )
      title(
        ylab = names(y_df),
        line = 0.5,
        font.lab = 2,
        cex.lab = 1.4
      )
      
      # Draw the confidence interval around the fitted line
      polygon(
        c(newx, rev(newx)),
        c(a[, 2], rev(a[, 3])),
        col = "grey91",
        border = TRUE,
        lty = "dashed"
      )
      
      # Samples are shown as dots
      points(x, y)
      
      # Draw linear regression line
      abline(clm, lwd = 2)
      
      # Take calculated pairwise p-value from rcorr
      pvalue_text <-
        paste("P-value:", round(as.numeric(corr_pval_cutoff[i, 4]), 4), sep = "")
      
      # Take corrected pairwise p-value
      pvalue_corr_text <-
        paste("Adj. p-value:", as.numeric(corr_pval_cutoff[i, 6]), sep = "")
      
      # Take calculated pairwise correlation coefficient from rcorr
      corr_text <-
        paste("Pearson's r:", round(as.numeric(corr_pval_cutoff[i, 3]), 4), sep = "")
      
      # Take calculated pairwise number of observations from rcorr
      support_text <-
        paste("supported by ",
              round(as.numeric(corr_pval_cutoff[i, 5]), 4),
              " observations",
              sep = "")
      
      # Show correlation coefficient and p-value in the plot
      mtext(corr_text, side = 3, line = 2)
      mtext(pvalue_text, side = 3, line = 1)
      mtext(pvalue_corr_text, side = 3, line = 0)
      mtext(support_text, side = 1, line = 2)
      
      abline(par("usr")[3], 0)
      segments(par("usr")[1], a[1, 2], par("usr")[1], a[1, 3])
      abline(par("usr")[4], 0)
      
    }
    
    dev.off()
    file.copy(
      paste0("6.Correlations/", newdir, "/linear_sign_pairs.pdf"),
      paste0("../www/", wd, "/linear_sign_pairs.pdf"),
      overwrite = T
    )
    
    output$Finalpdf <- renderUI({
      tags$iframe(
        src = paste0(wd, "/", "linear_sign_pairs.pdf"),
        height = "600",
        width = "100%"
      )
    })
    #################################################################################
    ######                        Write Output Files                           ######
    #################################################################################
    
    # Take the name of the inputfile to name the folder
    prefix = paste("OTUsCombined_Corr_input_table", sep = "_")
    
    # Make a directory name with inputfile name and date
    newdir2 <- paste(prefix, Sys.Date(), sep = "_")
    
    # Create a directory
    dir.create(paste0("6.Correlations/", newdir, "/", newdir2))
    
    # Write the log-scale transformed table
    write.table(
      my_scaled_data,
      paste0(
        "6.Correlations/",
        newdir,
        "/",
        newdir2,
        "/transformed.tab"
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Write the correlation table
    write.table(
      my_cor_matrix,
      paste0(
        "6.Correlations/",
        newdir,
        "/",
        newdir2,
        "/correlation-table.tab"
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Write the pvalue table
    write.table(
      my_pvl_matrix,
      paste0("6.Correlations/", newdir, "/", newdir2, "/pval-table.tab"),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Write the number of effective samples table
    write.table(
      my_num_matrix,
      paste0(
        "6.Correlations/",
        newdir,
        "/",
        newdir2,
        "/support-table.tab"
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Write the significant correlations
    write.table(
      my_pairs_cutoff,
      paste0(
        "6.Correlations/",
        newdir,
        "/",
        newdir2,
        "/cutoff-pairs-corr-sign.tab"
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    # Write plotted pairs
    write.table(
      corr_pval_cutoff,
      paste0(
        "6.Correlations/",
        newdir,
        "/",
        newdir2,
        "/plotted-pairs-stat.tab"
      ),
      sep = "\t",
      col.names = NA,
      quote = FALSE
    )
    
    if (!flag) {
      stop(
        "
        It was not possible to install all required R libraries properly.
        Please check the installation of all required libraries manually.\n
        Required libaries:ade4, GUniFrac, phangorn, randomcoloR, Rcpp"
      )
    }
    
    
    #################################################################################
    ######                           End of Script                             ######
    #################################################################################
    })
    
    setwd("../")
    })
  
  createZip <- function() {
    setwd(wd)
    fs <-
      c(
        dir("1.Normalization", full.names = T),
        dir("2.Alpha-Diversity", full.names = T),
        dir("3.Beta-Diversity", full.names = T),
        dir("4.Taxonomic-Binning", full.names = T),
        dir("5.Serial-Group-Comparisons", full.names = T),
        dir("6.Correlations", full.names = T)
      )
    zip(zipfile = "Rhea_data_download.zip", files = fs)
    setwd("../")
  }
  
  output$Download <-
    downloadHandler(
      filename = 'Rhea_data_download.zip',
      content = function(file) {
        createZip()
        file.copy(paste0(wd, "/Rhea_data_download.zip"), file, overwrite = T)
      },
      contentType = "application/zip"
    )
  
  }

# Run the application
shinyApp(ui = ui, server = server)