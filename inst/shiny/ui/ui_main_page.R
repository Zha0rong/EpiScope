
tabPanel("EpiScope",
         #img(src = "header.png", height="50%", width="50%", align="right"),

         # Application title
         titlePanel("Upload Bed files"),


         # Place for uploading data
         sidebarLayout(
           sidebarPanel(
             h3("Bed file Upload"),
             fileInput(
               "Bed",
               "Bedfiles",
               multiple = FALSE
             ),
            selectizeInput(inputId = 'Species',label = 'Species',choices = NULL,multiple = F,selected = NULL,
                        options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))),
            selectizeInput(inputId = 'GenomeVersion',label = 'GenomeVersion',choices = NULL,multiple = F,selected = NULL,
                        options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))),
            selectizeInput(inputId = 'AnnotationTitle',label = 'AnnotationTitle',choices = NULL,multiple = F,selected = NULL,
                        options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }'))),

            actionButton(inputId = 'submit',label = 'Submit')

           ),

           # Show a table of the inputted data
           mainPanel(
             tabsetPanel(
               tabPanel('Peak Summmary',
                        h3('Distribution of Domain in Genome'),
                        DT::dataTableOutput('chromosomal_peak_distribution',width = '75%')
               ),

               tabPanel('Peak Length Distribution',
                        plotOutput('length_distribution',width = '75%'),
                        DT::dataTableOutput(outputId = 'length_distribution_table',width='75%')

               ),
               tabPanel('TSS Heatmap',
                        h3('Heatmap visualization of TSS coverage'),
                        numericInput('TSS_range','Range of TSS site',value = 3000,min = 0),
                        actionButton(inputId = 'TSS_heatmap_submit',label = 'Submit'),
                        plotOutput('TSS_Heatmap',width = '50%')
               ),
               tabPanel('Functional Annotation of peaks',
                        plotOutput('upsetandvenn')
               ),
               tabPanel('Detail Annotation of peaks',
                        DTOutput(outputId = 'annotation_table',width='75%')
               )

             )
           )
         )
)
