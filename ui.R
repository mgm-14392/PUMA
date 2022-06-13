
library(shiny)

shinyUI(fluidPage(
  titlePanel("PUMA: Platform for Unified Molecular Analysis, Version 1.0"),
  tabsetPanel(
    tabPanel("About",
             sidebarLayout(position="left",
                           sidebarPanel(
                             p("If you use this app, please cite it as follows:"),
                             a("Gonzalez-Medina M, Medina-Franco JL. Platform for Unified Molecular Analysis: PUMA.
                               J. Chem. Inf. Model. (2017), 2017;57(8):1735-40.",
                               style="text-align:justify",
                               href= "http://pubs.acs.org/doi/10.1021/acs.jcim.7b00253"),
                             img(src="Logo.png", height = 150, width = 170, style="display: block; margin-left: auto; margin-right: auto;"),
                             a("Visit our website", href= "http://www.difacquim.com/", target="_blank"),
                             h6("This Shiny app was developped by", a("Mariana Gonzalez-Medina", 
                                                                      href="https://www.researchgate.net/profile/Mariana_Gonzalez-Medina2",
                                                                      target="_blank"), "and",
                                a("Jose L. Medina-Franco", 
                                  href="https://www.researchgate.net/profile/Jose_Medina-Franco", target="_blank"), style="text-align:justify"),
                             p("Last update July, 2017")
                             ),
                           mainPanel(
                             h3("Welcome to PUMA!"),
                             p("This is an online server to visualize 
                               the chemical space and compute the molecular diversity 
                               of your data sets. You only need a comma (,) delimited file with
                             the SMILES of all your compounds, the name of your data sets and an ID for each compound.", style="text-align:justify"),
                             p("You can download the User Guide, an example of an input file and template", 
                               a("here.", 
                                 href= "https://www.difacquim.com/d-tools/", target="_blank")),
                             h3("Contact"),
                             h4("We would appreciate your feedback!"),
                             p("If you have questions or suggestions,
                               please send an email to:"),
                             p(a("Mariana Gonzalez-Medina", 
                                 href="https://www.researchgate.net/profile/Mariana_Gonzalez-Medina2",
                                 target="_blank"),
                               strong("mgm_14392@comunidad.unam.mx")),
                             p(a("Jose L. Medina-Franco", 
                                 href="https://www.researchgate.net/profile/Jose_Medina-Franco", 
                                 target="_blank"),
                               strong("medinajl@unam.mx")),
                             strong("The non-shiny scripts are available upon request.",
                                    style="text-align:justify"),
                             h3("Acknowledgements"),
                             p("UNAM: PAPIME PE200116; PAIP 5000-9163")
                           )
             )),
    tabPanel("Chemical Space",
             sidebarLayout(position = "left",
                           sidebarPanel(
                             fileInput('file1', 'Choose a csv file',
                                       accept=c('text/csv', 'csv', 'comma-separated-values','.csv',
                                                'text/comma-separated-values/plain')),
                             selectInput("TOTD","Select 2D or 3D visualization",
                                                      choices = c("2D","3D")),
                             selectInput("PC1","Select a PCa",
                                         choices = c("PC1"=4,"PC2"=5,"PC3"=6,
                                                     "PC4"=7,"PC5"=8)),
                             selectInput("PC2","Select a PCb",
                                         choices = c("PC2"=5,"PC1"=4,"PC3"=6,
                                                     "PC4"=7,"PC5"=8)),
                             selectInput("PC3","Select a PCc",
                                         choices = c("PC3"=6,"PC1"=4,"PC2"=5,
                                                     "PC4"=7,"PC5"=8)),
                             downloadButton("down1","Download PCs and properties"),
                             downloadButton("down2","Download PCA loadings"),
                             downloadButton("down13","Download PCA summary")
                             
                           ),
                           mainPanel(plotlyOutput("plot",width = 700, height = 700),
                                     downloadButton("down","Download image")
                                     )
                           )
             ),
    tabPanel("Properties - Statistics",
             sidebarLayout(position = "left",
                           sidebarPanel(
                             fileInput('file2', 'Choose a csv file', 
                                       accept=c('text/csv', 'csv', 'comma-separated-values','.csv',
                                                'text/comma-separated-values/plain')),
                             selectInput("GRAF","Select a plot type",
                                         choices = c("Density","Boxplot","Histogram")),
                             selectInput("PROP","Select a molecular property",
                                         choices = c("MW"=4,"TPSA"=5,"nRotB"=6,
                                                     "nHBDon"=7,"nHBAcc"=8,"ALogP"=9)),
                             downloadButton("down4","Download property statistics"),
                             p("This subset includes compounds that satisfy all Lipinski's and
                               Veber's rules."),
                             downloadButton("downDL","Download drug-like subset"),
                             p("This subset includes compounds that violate no more than one of
                               Lipinski's rules."),
                             downloadButton("downDL2","Download subset")
                             ),
                           mainPanel(plotlyOutput("plot2", width = 700, height = 700),
                                     downloadButton("down5","Download image")
                                     )
                           )
             ),
    tabPanel("Properties - Similarity and Distance",
             sidebarLayout(position = "left",
                           sidebarPanel(
                             fileInput('file5', 'Choose a csv file', 
                                       accept=c('text/csv', 'csv', 'comma-separated-values','.csv',
                                                'text/comma-separated-values/plain')),
                             selectInput("SODM","Select a similarity or distance metric",
                                         choices = c("Euclidean","Tanimoto")),
                             p("File with pairwise similarity or distance.", 
                               style="text-align:justify"),
                             downloadButton("down10","Download the pairwise results"),
                             p("From the pairwise results, these compounds have the largest
                               distances (or lowest similarity) to the other compounds in the
                               data set.",style="text-align:justify"),
                             downloadButton("downdivs","Download diverse subset"),
                             p("File with inter and intra-data set distances.",
                               style="text-align:justify"),
                             downloadButton("down11","Download data sets results")
                           ),
                           mainPanel(plotlyOutput("plot5", width = 700, height = 700),
                                     downloadButton("down12","Download image")
                           )
             )
    ),
    tabPanel("Scaffolds - CSR Curves",
             sidebarLayout(
               sidebarPanel(
                 p("Cyclic System Recovery (CSR) curves"),
                 fileInput('file3',"Choose a csv file",
                           accept = c('text/csv', 'csv', 'comma-separated-values','.csv',
                                      'text/comma-separated-values/plain')
                           ),
                           downloadButton(
                             "down6","Download CSR curves data"),
                 downloadButton("downscaf","Download scaffolds")
                            ),
               mainPanel(
                 plotlyOutput("plot3", width = 700, height = 700),
                 downloadButton("down7","Download image")
                        )
             )
             ),
    tabPanel("Scaffolds - Shannon Entropy",
             sidebarLayout(
               sidebarPanel(
                 p("Scaled Shannon Entropy (SSE)"),
                 fileInput('file4', "Choose a csv file",
                           accept = c('text/csv', 'csv', 'comma-separated-values','.csv',
                                      'text/comma-separated-values/plain')
                           ),
                 p("Select the entropy of the n most populated scaffolds to generate the frequency plot.",
                   style="text-align:justify"),
                 selectInput("EI","SSEn",
                             choices = c("SSE10","SSE20","SSE30",
                                         "SSE40","SSE50","SSE60")),
                 p("Write the number of data set you want to plot.",
                   style="text-align:justify",strong("See the number in the SSE data file.",
                                                     style="text-align:justify")),
                 numericInput("EIN",
                              "Data set number",
                              min= 1, max = 30, value = 1, step = 1),
                 p("You will obtain the Scaled Shannon Entropy (SSE)
                   of the n most populated scaffolds (n up to 60).",
                   style="text-align:justify"),
                 downloadButton(
                   "down8","Download SSE data"),
                 downloadButton("dScaffs","Download scaffolds and ID"),
                 p("This subset includes all the unique scaffolds in the data set.",
                   style="text-align:justify"),
                 downloadButton("dunscafs","Download unique scaffolds")
               ),
               mainPanel(
                 plotlyOutput("plot4", width = 700, height = 700),
                 downloadButton("down9","Download image")
               )
             )
             ),
    tabPanel("Fingerprint Diversity",
             sidebarLayout(position = "left",
                           sidebarPanel(
                             fileInput("fileCDF", "Choose a csv file",
                                       accept = c('text/csv', 'csv', 'comma-separated-values','.csv',
                                                  'text/comma-separated-values/plain')),
                             selectInput("FP3",
                                         "Select a fingerprint",
                                         c("ECFP","Pubchem","MACCS")),
                             radioButtons("dep3",
                                          "Select the diameter for topological fingerprints",
                                          choices=c("4"=4, "6"=6)),
                             downloadButton("downCDF", "Download the similarity statistics")
               ),
               mainPanel(
                 plotlyOutput("CDFplot", width = 700, height = 700),
                 downloadButton("CDF","Download CDP plot")
               )
             )),
    tabPanel("D-Tools",
             h3("DIFACQUIM tools for Cheminformatics"),
             fluidRow(
               column(6,
                      a(img(src="CDPs.png", height = 300, width = 450),
                        href="http://132.248.103.152:3838/CDPlots/",
                        target="_blank")),
               column(6,
                      a(img(src="ALP.png", height = 300, width = 400),
                        href="http://132.248.103.152:3838/ActLSmaps/",
                        target="_blank"))
               ),
             fluidRow(
               column(6,
                      a("Consensus Diversity Plots (CDPs)",
                        href="http://132.248.103.152:3838/CDPlots/",
                        target="_blank")),
               column(6,
                      a("Activity Landscape Plotter (ALPloter)",
                        herf="http://132.248.103.152:3838/ActLSmaps/",
                        target="_blank")))
             )
    )
  ))

