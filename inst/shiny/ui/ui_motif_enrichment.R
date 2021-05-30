
tabPanel("Motif Search",

         # Application title
         titlePanel("Motif Search"),
         
         # Place for uploading data
         sidebarLayout(
           sidebarPanel(
             h3("Upload meme files"),
             tags$div(tags$p(
               'Currently meme files are accepted'
             )),
             fileInput('meme_input','Upload meme files',multiple = F),
             actionButton(inputId = 'meme_upload',label = 'Upload')
             
           ),
           
           # Show a table of the inputted data
           mainPanel(
             tabsetPanel(
               tabPanel('Motif Visualization',
                selectizeInput('Select_Motif','Select a motif to Visualize',multiple=F,choices=NULL,selected = NULL),
                plotOutput('motif_visualization',width = '75%')
               ),
               tabPanel('Motif Enrichment',

               )

             )
           )
         )
)