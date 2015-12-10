shinyUI(pageWithSidebar(
    headerPanel(h4("Rename tips in a phylogenetic tree:")),
    sidebarPanel(
        img(src = "DarwinsTree.jpg", height = 50, width = 50),
        fileInput('file1', label="1) Select the phylo tree:"),
        radioButtons("tree_informat", 
                           label = "Input phylo tree format", 
                           choices = c(newick = "newick", 
                                          nexus = "nexus", BEAST = "beast"), selected = "newick", inline = TRUE),
        hr(),
        helpText("One column of the metadata file should contain the tree tip-labels"),
        fileInput('file2', label="2) Select the metadata file:",
                  accept=c('text/csv', 'text/comma-separated-values,text/plain',
                           '.csv')),
        helpText("Select your metadata file header, delimiter and quote style (click the 'Metadata table' tab to preview):"),          
        checkboxInput('header', 'Header row?', TRUE),
                  radioButtons('sep', 'Delimiter',
                      c(comma=',',
                       semicolon=';',
                       tab='\t'),
                       selected = ',', inline=TRUE),
                  radioButtons('quote', 'Quote style',
                     c(None='',
                       'Double Quote'='"',
                       'Single Quote'="'"),
                     selected = '"', inline=TRUE),
        hr(),
        uiOutput("choose_original_label_column"),
        uiOutput("choose_new_label_columns"),
        helpText("Click 'download' only after viewing tree in the 'Processed tree' tab"),
        downloadButton("downloadData")
    ),
    mainPanel(
        textOutput("text1"),
        textOutput("text2"),
        tabsetPanel(
            tabPanel("Input tree", plotOutput('tree_in', height=700)),
            tabPanel("Metadata table", tableOutput('metadata')),
            tabPanel("Processed tree", plotOutput('tree_out', height=700))
        )
    )
))
