library(data.table)
library(ggplot2)
library(viridis)
library(VennDiagram)

annots <- fread("cricket_emapper.emapper.annotations", skip=4, fill=T, na.string="-")
table_annots <- subset(annots, !(grepl("## ", annots$`#query`)))

##all have description - not all have GO, KEGG, PFAM

length(!is.na(table_annots$PFAMs))

##mean evalue
mean(table_annots$evalue)
min(table_annots$evalue)
max(table_annots$evalue)
##evalue results look horrible if you plot it

total.number.genes <- 191455
##how many unique IDs in query column
number.egnogg.annot.genes <- length(unique(table_annots$`#query`))

##genes with GO, KEGG, PFAM
##find rows where value isn't NA (to remove genes with no annotation)
##and then make list of querys/genes that had that an annotation of that category
go <- table_annots[!is.na(GOs), unique(`#query`)]
kegg <- table_annots[!is.na(KEGG_Pathway), unique(`#query`)]
pfam <- table_annots[!is.na(PFAMs), unique(`#query`)]

##venn diagram overlap
vd <- venn.diagram(x=list("Pfam"=pfam, "KEGG"=kegg, "GO"=go), filename=NULL,
                   fill=c("#440154FF", "#21908CFF", "#FDE725FF"), alpha=0.7,
                   main=paste("Total number of EggNOG annotated transcripts:", number.egnogg.annot.genes))
grid.newpage()
grid.draw(vd)

##you can save the file from the plot window (as a jpeg/pdf) or with this as an svg file
##you can open and edit svg files in inkscape which can be handy sometimes
ggsave(file="venn_diagram.svg", plot=vd, width=10, height=8)


##what percentage of genes had an eggnogg annotation?
(number.egnogg.annot.genes/total.number.genes)*100

