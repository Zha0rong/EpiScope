
options(shiny.maxRequestSize=600*1024^2)
source('server/observer.R',local = T)


AnnotationHub=AnnotationHub()

Species=as.character(unique(AnnotationHub$species))

reactivevalue=reactiveValues(Bed=NULL,
                             peak=NULL,
                             peakAnno=NULL,
                             txdb=NULL,
                             peakAnnodataframe=NULL,
                             SymboltoID=NULL,
                             motifs=NULL,
                             AnnotationHub=AnnotationHub,
                             species=Species,
                             select_species=NULL,
                             GenomeVersion=NULL,
                             SelectedAnnotation=NULL,
                             annodb=NULL
                             )



