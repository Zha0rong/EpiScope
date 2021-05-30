observe( if (!is.null(reactivevalue$species)) {
  updateSelectizeInput(session,'Species','Select Species',choices=reactivevalue$species,
                       selected=NULL,options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))

})

observe( if (!is.null(input$Species)) {
  reactivevalue$select_species=input$Species
  msg <- sprintf('Retrieving Genome Version...')
  withProgress(message=msg, {
    setProgress(0.25)
    available_genome_version=unique(query(reactivevalue$AnnotationHub,c('GRanges',reactivevalue$select_species))$genome)

    setProgress(1, 'Completed')

  })

  updateSelectizeInput(session,'GenomeVersion','Select Genome Version',choices=available_genome_version[!is.na(available_genome_version)],
                       selected=NULL,options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))

})

observe( if (!is.null(input$GenomeVersion)) {
  reactivevalue$GenomeVersion=input$GenomeVersion


  msg <- sprintf('Retrieving Annotation Files...')
  withProgress(message=msg, {
    setProgress(0.25)
    available_annotation=unique(query(reactivevalue$AnnotationHub,c('GRanges',reactivevalue$GenomeVersion,reactivevalue$select_species))$title)

    setProgress(1, 'Completed')

  })



  updateSelectizeInput(session,'AnnotationTitle','Select AnnotationTitle',choices=available_annotation[!is.na(available_annotation)],
                    selected=NULL,options=list(placeholder = 'Please select an option below',onInitialize = I('function() { this.setValue(""); }')))

})

observe( if (!is.null(input$AnnotationTitle)) {
  reactivevalue$SelectedAnnotation=input$AnnotationTitle
})

observeEvent( input$submit, {
  if (!is.null(input$Bed)) {
    msg <- sprintf('Uploading...')

    withProgress(message=msg, {
      reactivevalue$Bed=input$Bed$datapath
      reactivevalue$peak=readPeakFile(reactivevalue$Bed)
      seqlevelsStyle(reactivevalue$peak) = 'NCBI'
      setProgress(0.25, 'Finished Reading in Peaks.')
      reactivevalue$txdb=query(reactivevalue$AnnotationHub,c('GRanges',reactivevalue$select_species,
                                                             reactivevalue$GenomeVersion,
                                                             reactivevalue$AnnotationTitle))[[1]]

      reactivevalue$txdb=makeTxDbFromGRanges(reactivevalue$txdb,drop.stop.codons = T)
      reactivevalue$annodb=query(reactivevalue$AnnotationHub,c('OrgDb',reactivevalue$select_species))[[1]]

      setProgress(0.5, 'Finished Building Annotation.')
      reactivevalue$peakAnno <- annotatePeak_local_AnnoDb(reactivevalue$peak, tssRegion=c(-3000, 3000),level = 'gene',
                                             TxDb=reactivevalue$txdb,verbose = F,annoDb = reactivevalue$annodb)
      reactivevalue$geneannotation=ExtractGeneInfo(reactivevalue$annodb)
      setProgress(0.75, 'Finished Annotating Peaks.')
      output$upsetandvenn = renderPlot(upsetplot(reactivevalue$peakAnno, vennpie=TRUE))
      setProgress(1, 'Completed')

    })
    reactivevalue$peakAnnodataframe=data.frame(reactivevalue$peakAnno)
  }
})


observe( if (!is.null(reactivevalue$peak)) {
  chromosomal_peak_distribution=data.frame(table(reactivevalue$peak@seqnames),stringsAsFactors = F)
  colnames(chromosomal_peak_distribution)=c('Chromosome','Number of Peaks')
  chromosomal_peak_distribution=chromosomal_peak_distribution[order(as.character(chromosomal_peak_distribution$Chromosome)),]
  output$chromosomal_peak_distribution=renderDataTable(datatable(chromosomal_peak_distribution))
  length_distribution=as.numeric(reactivevalue$peak@ranges@width)
  length_distribution=data.frame(length_distribution)
  colnames(length_distribution)='width'
  length_distribution_table=data.frame(as.numeric(summary(length_distribution$width)))
  rownames(length_distribution_table)=c( 'Minimum','First Quartile',  'Median',    'Mean', 'Third Quartile',    'Maximum')
  colnames(length_distribution_table)='Statistics'
  length_distribution_table=datatable(length_distribution_table)
  plot=ggplot(length_distribution,aes(x=width))+
    geom_histogram(aes(y=..density..),
                   binwidth=100,
                   colour="black", fill="white") +
    geom_density(alpha=.2, fill="#FF6666")
  output$length_distribution=renderPlot(plot)
  output$length_distribution_table=renderDataTable(length_distribution_table)
  selectedcolumns=c('seqnames',
                    'start',
                    'end',
                    'width',
                    'annotation',
                    'geneId',
                    'distanceToTSS',
                    'SYMBOL',
                    'GENENAME')
  output$annotation_table=DT::renderDT(DT::datatable(reactivevalue$peakAnnodataframe[,selectedcolumns],rownames = F),filter = "top")

})

#


observeEvent( input$TSS_heatmap_submit, {
  output$TSS_Heatmap=renderPlot(peakHeatmap(reactivevalue$peak, TxDb=reactivevalue$txdb, upstream=(input$TSS_range), downstream=input$TSS_range, color="blue"))
})





observeEvent( input$meme_upload, {
  if (!is.null(input$meme_input)) {
    motifs=read_meme(input$meme_input$datapath)
    reactivevalue$motifs=list()
    for (i in 1:length(motifs)) {
      reactivevalue$motifs[[paste(motifs[[i]]@altname,motifs[[i]]@name)]]=motifs[[i]]
    }
    updateSelectizeInput(session,'Select_Motif','Select a motif to Visualize',choices=names(reactivevalue$motifs),selected = NULL)


  }
})

observeEvent(input$Select_Motif, {
  if (!is.null(input$Select_Motif)) {
    tryCatch({output$motif_visualization=renderPlot(plotMotifLogo(reactivevalue$motifs[[input$Select_Motif]]@motif,
                                                        motifName = input$Select_Motif))})
  }
})




