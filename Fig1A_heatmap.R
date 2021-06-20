
######################################
#### Read in the files
######################################

options(stringsAsFactors = FALSE)
library(ComplexHeatmap)
library(svglite)
library(RColorBrewer)

read.in <- function(path.in, filename.in){
  df <- read.delim(paste0(path.in,filename.in))
  samples <- df[1,]
  samples <- colnames(samples[,samples[1,] =="E" | samples[1,] =="#!{Type}E"])
  df.annotation <- df[grep("#",df[,1]),samples]
  df.annotation[,1] <- apply(as.matrix(df.annotation[,1]), 1, function(x){
    y <- strsplit(x, split="}")[[1]][2]
  })
  row.names(df.annotation) <- apply(as.matrix(df[grep("#",df[,1]),1]), 1, function(x){
    y <- strsplit(x, split=":|}")[[1]][2]
  })
  df.annotation <- df.annotation[-1,]
  df.annotation <- apply(df.annotation, c(1,2), function(x){
    if(grepl(",",x)){
      as.numeric(gsub(",", ".", gsub("\\.", "", x)))
    } else {
      x
    }
  })
  
  df.expression <- df[-grep("#",df[,1]),samples]
  row.names(df.expression) <- df$Name[-grep("#",df[,1])]
  df.expression <- apply(df.expression, c(1,2),as.numeric)
  df.expression[is.nan(df.expression)] <- NA
  
  sample_annotation <- as.data.frame(t(df.annotation))
  return(list(Expr = df.expression, Annot = sample_annotation))
}

path.in <- "C:/Users/User/PhD/MM_segundo"

rawProteomics <- read.delim(paste0(path.in, 
                                   "20200310_Global Proteomics_Top 500 Anova FDR 0,005 with annotations.txt"))
rawPhospho <- read.delim(paste0(path.in, 
                                "20200310_Phosphoproteomics_Top 500 Anova FDR 0,05.txt"))
rawTransc <- read.delim(paste0(path.in, 
                               "Top 500 Anova FDR 0,05 Transcriptomics_5 subtypes.txt"))


Proteomics <- read.in(path.in,
                      "20200310_Global Proteomics_Top 500 Anova FDR 0,005 with annotations.txt")
Phospho <- read.in(path.in,
                   "20200310_Phosphoproteomics_Top 500 Anova FDR 0,05.txt")
Transc <- read.in(path.in,
                  "Top 500 Anova FDR 0,05 Transcriptomics_5 subtypes.txt")



######################################
#### Clean the annotation tables
######################################

PClusters.Prot <- Proteomics$Annot$Cluster_
names(PClusters.Prot) <- row.names(Proteomics$Annot)
PClusters.Phospho <- Phospho$Annot$Cluster_
names(PClusters.Phospho) <- row.names(Phospho$Annot)
names(PClusters.Phospho) %in% names(PClusters.Prot)
PClusters.Transc <- Transc$Annot$Cluster_
names(PClusters.Transc) <- row.names(Transc$Annot)
names(PClusters.Transc) %in% names(PClusters.Prot)
names(PClusters.Transc) %in% names(PClusters.Prot)

Common.Annot <- Proteomics$Annot

Common.Annot <- as.data.frame(apply(Common.Annot, 2, function(x){
  y <- ifelse(x =="", NA, x)
  y <- gsub("#N/A", NA, y, fixed = T)
  y <- ifelse(y =="-", NA, y)
  y <- ifelse(y =="??", NA, y)
}), stringsAsFactors = F)

Common.Annot$Subtype <- Common.Annot$Cluster_
Common.Annot$`Age at diagnosis` <- cut(as.numeric(Common.Annot$age.diag), breaks = quantile(as.numeric(Common.Annot$age.diag), na.rm=T),
                                       include.lowest = T)
Common.Annot$survival <- as.numeric(Common.Annot$survival)
Common.Annot$`Disease specific survival` <- cut(as.numeric(Common.Annot$survival)/365.242199, 
                                                breaks = c(0,1,3,5,26), include.lowest = T)

Common.Annot$Gender <- Common.Annot$gender
Common.Annot$`Disease stage` <- Common.Annot$dis.stage
Common.Annot$`Status` <- Common.Annot$`Remarks `
Common.Annot$`Adjacent lymphnode (%)` <- as.numeric(Common.Annot$`adjacent lymphnode`)
Common.Annot$`Average tumor content (%)` <- as.numeric(Common.Annot$`Average tumor`)
Common.Annot$`Lymphatic score` <- as.numeric(Common.Annot$`Average lymphatic`)

Common.Annot$Lund <- ifelse(Common.Annot$Lund =="High.imm", "High immune", Common.Annot$Lund)
colnames(Common.Annot)[colnames(Common.Annot)=="Lund"] <- "Cirenajwis et al. 2015"

Common.Annot$`Adjacent lymphnode (%)` <- ifelse(is.na(Common.Annot$`Adjacent lymphnode (%)`), 
                                                0, Common.Annot$`Adjacent lymphnode (%)`)
Common.Annot$`Adjacent lymphnode2` <- 0
Common.Annot$`Average tumor content (%)` <- ifelse(is.na(Common.Annot$`Average tumor content (%)`), 
                                                   0, Common.Annot$`Average tumor content (%)`)

for (i in 1:nrow(Common.Annot)){
  if (is.na(Common.Annot$`Lymphatic score`[i])){
    Common.Annot$`Lymphatic score`[i] <- 0
    Common.Annot$`Lymphatic score2`[i] <- -0.5
  } else {
    Common.Annot$`Lymphatic score2`[i] <- 0
  }
}
Common.Annot$TCGA <- ifelse(Common.Annot$TCGA =="keratin", "Keratin",Common.Annot$TCGA)
Common.Annot$TCGA <- ifelse(Common.Annot$TCGA =="MITF.low", "MITF low",Common.Annot$TCGA)
Common.Annot$TCGA <- ifelse(Common.Annot$TCGA =="immune", "Immune",Common.Annot$TCGA)



######################################
#### Clean the expression tables
######################################

Expr.tr <- scale(t(Transc$Expr))
Expr.pr <- t(Proteomics$Expr)
Expr.ph <- scale(t(Phospho$Expr))


tr.missingpat <- row.names(Expr.pr)[!row.names(Expr.pr) %in% row.names(Expr.tr)]
ph.missingpat <- row.names(Expr.pr)[!row.names(Expr.pr) %in% row.names(Expr.ph)]

missing <- c("Mm1","Mm2", "Mm3","Mm4","Mm5",
             "Mm6","Mm7","Mm8", "Mm9",
             "Mm10", "Mm11", "Mm12", "Mm13",
             "Mm14","Mm15", "Mm16","Mm17","Mm18","Mm19","Mm20",  "Mm21", "Mm22", "Mm23",
             "Mm24", "Mm25",  
             "Mm26", "Mm27", "Mm28", "Mm29", "Mm30")
not.missing <- row.names(Expr.pr)[!row.names(Expr.pr) %in% missing]

#length(missing) + length(not.missing)



EC.order <- c(row.names(Common.Annot[Common.Annot$Subtype=="EC" &
                                       row.names(Common.Annot) %in% not.missing, ]),
              c("Mm1","Mm2", "Mm3","Mm4","Mm5"))

EC.I.order <- c(row.names(Common.Annot[Common.Annot$Subtype=="EC-Immune" &
                                         row.names(Common.Annot) %in% not.missing, ]),
                c("Mm6","Mm7","Mm8", "Mm9"))


EC.M.order <- c(row.names(Common.Annot[Common.Annot$Subtype=="EC-Mit" &
                                         row.names(Common.Annot) %in% not.missing, ]),
                c( "Mm10", "Mm11","Mm12", "Mm13"))

Mit.order <- c(row.names(Common.Annot[Common.Annot$Subtype=="Mit" &
                                        row.names(Common.Annot) %in% not.missing, ]),
               c("Mm14", "Mm15","Mm16", "Mm17","Mm18","Mm19","Mm20",
                 "Mm21", "Mm22", "M23","Mm24", "Mm25"))

Mit.I.order <- c(row.names(Common.Annot[Common.Annot$Subtype=="Mit-Immune" &
                                          row.names(Common.Annot) %in% not.missing, ]),
                 c("Mm26", "Mm28", "Mm29", "Mm30", "Mm31"))


Expr.pr <- Expr.pr[c(EC.order,EC.I.order,
                     EC.M.order, Mit.order, Mit.I.order),]


add.decoy <- function(df, missing, order){
  Decoy.table <- df[1:length(missing),]
  Decoy.table[!is.na(Decoy.table)] <- NA
  Decoy.table[is.na(Decoy.table)] <- -20
  row.names(Decoy.table) <- missing
  df.suppl <- rbind(df, Decoy.table)
  df.suppl <- df.suppl[order,]
  return(df.suppl)
}
Expr.ph <- add.decoy(Expr.ph, missing=ph.missingpat, order=row.names(Expr.pr))
Expr.tr <- add.decoy(Expr.tr, missing=tr.missingpat, order=row.names(Expr.pr))


all(row.names(Expr.ph) == row.names(Expr.pr))
all(row.names(Expr.tr) == row.names(Expr.pr))

Common.Expr.table <- rbind(t(Expr.pr),t(Expr.ph),t(Expr.tr))

data.types <- c(rep("Proteomics", ncol(Expr.pr)), rep("Phosphoproteomics", ncol(Expr.ph)),
                rep("Transcriptomics", ncol(Expr.tr)))



Common.Annot <- Common.Annot[row.names(Expr.pr),]


######################################
#### Create the column annotation and the heatmaps
######################################

column_ha = HeatmapAnnotation(df =Common.Annot[,c("Subtype", "Cirenajwis et al. 2015", "TCGA", 
                                                  "BRAF", "NRAS", "WT",
                                                  "Age at diagnosis", "Gender",
                                                  "Disease specific survival", "Disease stage", "Status")],
                              which="col",
                              col=list(
                                "Subtype" =
                                  c("EC"="#FBB4AE",
                                    "Mit"="#B3CDE3",
                                    "EC-Immune" ="#CCEBC5",
                                    "EC-Mit"="#DECBE4",
                                    "Mit-Immune"="#FED9A6"),
                                "Cirenajwis et al. 2015" = 
                                  c("High immune"="#d8407d",
                                    "Normal"="#e4fb65",
                                    "Pigmentation"="#5dd1a8",
                                    "Proliferative"="#6f7ebd"),
                                "TCGA" = 
                                  c("Immune"="#fd8f66",
                                    "Keratin"="#e4fb65",
                                    "MITF low"="#6a408d"),
                                "Disease specific survival" = 
                                  c("[0,1]"="white", 
                                    "(1,3]"="#7FCDBB",
                                    "(3,5]" ="#1D91C0",
                                    "(5,26]"="#253494"),
                                "Status"=
                                  c("dead" = "#33A02C",
                                    "dead (likely melanoma)" = "#B2DF8A",
                                    "dead other reason" = "#FFFF99",
                                    "dead unknown reason" = "#6A3D9A",
                                    "alive" = "white"),
                                "Disease stage"=
                                  c("1" = "white",
                                    "3" = "#FB9A99",
                                    "4" = "#E31A1C"),
                                "BRAF"= 
                                  c("0" = "White",
                                    "1" = "#CE1256"),
                                "NRAS"= 
                                  c("0" = "White",
                                    "1" = "#CE1256",
                                    "NA" ="grey"),  
                                "WT"= 
                                  c("0" = "white",
                                    "1" = "#666666",
                                    "NA" ="grey"),
                                "Age at diagnosis"= 
                                  c("[24,54.5]" = "white",
                                    "(54.5,64]" ="#FDBF6F",
                                    "(64,71.5]"="#FF7F00",
                                    "(71.5,86]" ="#B15928"),
                                "Gender"= 
                                  c("Female" = "#666666",
                                    "Male" = "white")),
                              
                              "Average tumor and\nadjacent lymphnode content" = anno_barplot(Common.Annot[,c("Average tumor content (%)",
                                                                                                             "Adjacent lymphnode (%)",
                                                                                                             "Adjacent lymphnode2")],
                                                                                             gp = gpar(fill = c("#5E4FA2", "#FF7F00", "white")),
                                                                                             height = unit(1.5, "cm"),
                                                                                             axis_param = list(
                                                                                               side = "left",
                                                                                               at = c(30, 70), 
                                                                                               labels = c("30%", "70%") )),
                              "Lymphatic score" = anno_barplot(Common.Annot[,c("Lymphatic score", #Set3
                                                                               "Lymphatic score2")],
                                                               gp = gpar(fill = c("#FDBF6F", "grey")),
                                                               height = unit(1.5, "cm"),
                                                               axis_param = list(
                                                                 side = "left",
                                                                 at = c(0, 2, 4, 6), 
                                                                 labels = c("0","2", "4", "6") )),
                              annotation_name_side = "right",
                              gap = unit(c(3.5, 0.5, 2.5, 0.5, 0.5, 2.5, 0.5,0.5,0.5,0.5, 2.5, 0.5), 'mm'),
                              gp = gpar(col = "black"),
                              show_legend = FALSE,
                              show_annotation_name = TRUE)



ht1 <- Heatmap(t(Expr.pr), name = "Z-score", 
               top_annotation = column_ha, 
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               show_row_names = FALSE, 
               height = 5,
               column_split = Common.Annot$Subtype,
               cluster_column_slices = FALSE,
               cluster_row_slices  = FALSE,
               row_gap = unit(2, "mm"),
               column_gap = unit(2, "mm"),
               clustering_distance_rows = "pearson",
               clustering_method_rows = "average",
               clustering_distance_columns = "pearson",
               clustering_method_columns = "average",
               col = circlize::colorRamp2(c(-20,-10,-1, 0, 1), 
                                          c("white","#0000FF","#0000FF", "#FFFFFF", "#FF0000")),
               row_title_rot = 0, 
               show_column_dend = FALSE, show_row_dend = FALSE,
               row_title_side = "right", show_heatmap_legend = FALSE)


ht2 <- Heatmap(t(Expr.ph), name = "Z-score", 
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               show_row_names = FALSE, 
               height = 2,
               column_split = Common.Annot$Subtype,
               cluster_column_slices = FALSE,
               cluster_row_slices  = FALSE,
               row_gap = unit(2, "mm"),
               column_gap = unit(2, "mm"),
               clustering_distance_rows = "pearson",
               clustering_method_rows = "average",
               clustering_distance_columns = "pearson",
               clustering_method_columns = "average",
               col = circlize::colorRamp2(c(-20,-10,-1, 0, 1), 
                                          c("white","#0000FF","#0000FF", "#FFFFFF", "#FF0000")),
               row_title_rot = 0,
               show_column_dend = FALSE, show_row_dend = FALSE,
               row_title_side = "right", show_heatmap_legend = FALSE)

ht3 <- Heatmap(t(Expr.tr), name = "Z-score", 
               cluster_rows = TRUE,
               height = 2,
               cluster_columns = FALSE,
               show_row_names = FALSE, 
               show_column_names = FALSE,
               column_split = Common.Annot$Subtype,
               cluster_column_slices = FALSE,
               cluster_row_slices  = FALSE,
               row_gap = unit(2, "mm"),
               column_gap = unit(2, "mm"),
               clustering_distance_rows = "pearson",
               clustering_method_rows = "average",
               clustering_distance_columns = "pearson",
               clustering_method_columns = "average",
               col = circlize::colorRamp2(c(-20,-10,-1, 0, 1), 
                                          c("white","#0000FF","#0000FF", "#FFFFFF", "#FF0000")),
               row_title_rot = 0, 
               show_column_dend = FALSE, show_row_dend = FALSE,
               row_title_side = "right", show_heatmap_legend = FALSE)



######################################
#### Draw (and save) the heatmap list
######################################

ht_list = ht1 %v% ht2 %v% ht3
draw(ht_list)


#svg(filename="C:/Users/User/PhD/MM_segundo/Fig1heatmap.svg",width = 20, height = 15)
#ht_list
#dev.off()


######################################
#### Draw (and save) the legends
######################################

lgd <- Legend(col_fun = circlize::colorRamp2(c(-1, 0, 1), c("#0000FF", "#FFFFFF", "#FF0000")), 
              title = "Scaled Expression", at = c(-1, 0, 1),border = "black")

lgd.Subtype <- Legend(at = c("EC", "Mit","EC-Immune","EC-Mit", "Mit-Immune"), 
                      legend_gp = gpar(fill = c("#FBB4AE","#B3CDE3","#CCEBC5", "#DECBE4", "#FED9A6")), 
                      title = "Subtype\n ",border = "black",nr = 2)

lgd.Lund <- Legend(at = c("High immune", "Normal","Pigmentation","Proliferative"), 
                   legend_gp = gpar(fill = c("#D8407D","#E4FB65","#5DD1a8", "#6F7EBD")), 
                   title = "Cirenajwis et al. 2015\n ",border = "black",nr = 2)

lgd.TCGA <- Legend(at = c("Immune", "Keratin","MITF low"), 
                   legend_gp = gpar(fill = c("#FD8F66","#eFFB65","#6A408D")), 
                   title = "TCGA\n ",border = "black",nr = 2)

lgd.DSS <- Legend(at = c("[0,1]", "(1,3]","(3,5]","(5,26]"),
                  legend_gp = gpar(fill = c("white","#7FCDBB","#1D91C0","#253494")), 
                  title = "Disease specific survival (yrs)\n ",border = "black",nr = 2)

lgd.status <- Legend(at = c("dead", "dead (likely melanoma)",
                            "dead other reason","dead unknown reason",
                            "alive"),
                     legend_gp = gpar(fill = c("#33A02C","#B2DF8A","#FFFF99","#6A3D9A","white")), 
                     title = "Status\n ",border = "black",nr = 2)

lgd.stage <- Legend(at = c("1", "3","4"), 
                    legend_gp = gpar(fill = c("white","#FB9A99","#E31A1C")), 
                    title = "Disease stage\n ",border = "black",nr = 2)

lgd.BRAF <- Legend(at = c("0", "1"), 
                   legend_gp = gpar(fill = c("white","#CE1256")), 
                   title = "BRAF\n ",border = "black",nr = 2)

lgd.NRAS <- Legend(at = c("0", "1"), 
                   legend_gp = gpar(fill = c("white","#CE1256")), 
                   title = "NRAS\n ",border = "black",nr = 2)

lgd.WT <- Legend(at = c("0", "1"), 
                 legend_gp = gpar(fill = c("white","#666666")), 
                 title = "WT\n ",border = "black",nr = 2)

lgd.Age <- Legend(at = c("[24,54.5]", "(54.5,64]","(64,71.5]","(71.5,86]"),
                  legend_gp = gpar(fill = c("white","#FDBF6F","#FF7F00","#B15928")), 
                  title = "Age at diagnosis\n ",border = "black",nr = 2)


lgd.Gender <- Legend(at = c("Female", "Male"), 
                     legend_gp = gpar(fill = c("#666666","white")), 
                     title = "Gender\n ",border = "black",nr = 2)


lgd.barplot = Legend(at = c("Lymphatic score", "Average tumor content (%)", 
                            "Adjacent lymphnode (%)"), 
                     title = "Barplots\n", 
                     legend_gp = gpar(fill = c("#FDBF6F", "#5E4FA2", "#FF7F00")),
                     border = "black",nr = 2)


pd = packLegend(lgd, lgd.Subtype, lgd.Lund, lgd.TCGA, 
                lgd.BRAF, lgd.NRAS, lgd.WT,
                lgd.Age, lgd.Gender,lgd.DSS,lgd.stage,lgd.status,   
                lgd.barplot, direction = "horizontal",
                max_width = unit(15, "cm"),column_gap = unit(10, "mm"), row_gap = unit(10, "mm"))

draw(pd)



#svg(filename="C:/Users/User/PhD/MM_segundo/Fig1_heatmap_legends.svg",width = 10, height = 8)
#draw(pd)
#dev.off()
