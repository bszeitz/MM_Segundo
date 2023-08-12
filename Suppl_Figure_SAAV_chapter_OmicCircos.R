
#### Suppl.Figure for SAAV chapter (OmicCircos plot)

options(stringsAsFactors = F)
library(OmicCircos)
library(readxl)


setwd("C:/Users/User/PhD/MM_segundo/")

excelfile <- "SAAV_chapter_SupplTable_final.xlsx"
Verified.SAAVs <- as.data.frame(read_xlsx(excelfile, sheet=2))

Verified.SAAVs$Name <- paste(Verified.SAAVs$Master.Protein.Accession, Verified.SAAVs$Mutation.in.isoform, sep="_")

## seg.f.mydata: the segment data, which is used to draw the outmost anchor track
## column 1: segment names, 2 and 3: start and end positions, 4: additional info of the segments (optional)
seg.f.mydata <- matrix(nrow=0, ncol=4)
colnames(seg.f.mydata) <- c("No.of.Batches", "Start", "End", "Name")
#no.seq <- c(1,4,7,10,13, 16)
no.seq <- base::seq(1,15,1)
#cat.names <- c("One-Three", "Four-Seven", "Eight-Ten", "Eleven-Thirteen", "Sixteen")
Verified.SAAVs <- Verified.SAAVs[order(Verified.SAAVs$Master.Protein.Accession),]
for (i in 1:length(no.seq)) {
  table.sub <- Verified.SAAVs[which( Verified.SAAVs$No.of.Batches == no.seq[i] ),c("Name", "No.of.Batches")]
  table.sub <- table.sub[order(table.sub$No.of.Batches),]
  table.sub[, "Start"] <- seq(0, nrow(table.sub)-1, 1)
  table.sub[, "End"] <- seq(1, nrow(table.sub), 1)
  table.sub <- table.sub[,c("No.of.Batches", "Start", "End", "Name")]
  seg.f.mydata <- rbind(seg.f.mydata, table.sub)
}


## seg.v.mydata: the mapping data, a data frame, which includes the values to be drawn in the graph
## column 1: segment names, 2: position, 3 and beyond: values and names (optional)
seg.v.mydata <- seg.f.mydata[,c("No.of.Batches", "End", "Name")]
for (i in 1:nrow(Verified.SAAVs)){
  Verified.SAAVs$Cosmic.Found[i] <- ifelse(Verified.SAAVs$Cosmic_ID[i] == "Not Found", 0, 10)
  peptatlas.match <- sum(as.numeric(strsplit(Verified.SAAVs$No.PeptideAtlas.ExactMatch[i], split=" + ", fixed = T)[[1]])) +
    sum(as.numeric(strsplit(Verified.SAAVs$No.PeptideAtlas.PartialMatch[i], split=" + ", fixed = T)[[1]]))
  Verified.SAAVs$PeptideAtlas.Found[i] <- ifelse(peptatlas.match == 0, 0, 10)
  Verified.SAAVs$cs.IDs[i] <- ifelse(Verified.SAAVs$cs.IDs[i] == "Not Found", 0, 10)
}
Verified.SAAVs.short <- Verified.SAAVs[,c("Name", "SAAV.PSM.sum", "PSM.ratio", 
                                          "PeptideAtlas.Found", "cs.IDs")]
seg.v.mydata <- merge(seg.v.mydata, Verified.SAAVs.short, by="Name")
seg.v.mydata <- seg.v.mydata[,c(2,3,1,4:(ncol(seg.v.mydata)))]
seg.v.mydata$cs.IDs <- as.numeric(as.character(seg.v.mydata$cs.IDs))

sapply(seg.v.mydata, class)
colnames(seg.f.mydata) <- c("seg.name", "seg.Start", "seg.End", "the.v")
seg.f.mydata$NO <- NA

seg.num <- length(unique(seg.f.mydata[,1]))
seg.name <- paste("B", 1:seg.num, sep="")
seg.f.mydata[,1] <- paste("B", seg.f.mydata[,1], sep="")
#new
seg.v.mydata <- seg.v.mydata[order(seg.v.mydata$No.of.Batches),]

summary(seg.v.mydata$PSM.ratio)



#colors <- rainbow(seg.num, alpha=0.5)

seg.v.mydata[is.na(seg.v.mydata$PSM.ratio),"PSM.ratio"] <- -0.3
col.ratio <- vector()
for (i in 1:nrow(seg.v.mydata)) {
  if (is.na(seg.v.mydata[i,"PSM.ratio"])) { col.ratio[i] <- NA
  } else if (seg.v.mydata[i,"PSM.ratio"] >= 0.50) { col.ratio[i] <- "#E31A1C" #red
  } else if (seg.v.mydata[i,"PSM.ratio"] >= 0.30) { col.ratio[i] <- "#FF7F00" #orange
  } else if (seg.v.mydata[i,"PSM.ratio"] >= 0.002) { col.ratio[i] <- "#1F78B4" #blue
  } else { col.ratio[i] <- "grey"}
} 
names(col.ratio) <- seg.v.mydata[,3]
seg.v.mydata$No.of.Batches <- paste("B", seg.v.mydata$No.of.Batches, sep="")
seg.v.mydata$SAAV.PSM.sum <- log10(seg.v.mydata$SAAV.PSM.sum)


name.order <- seg.f.mydata$the.v
name.row.nos <- vector()
for (i in 1:nrow(seg.v.mydata)){
  name.row.nos <- c(name.row.nos, grep(paste0("\\b", seg.v.mydata$Name[i], "\\b"), name.order))
}

row.names(seg.v.mydata) <- seq(1, nrow(seg.v.mydata), 1)
seg.v.mydata <- seg.v.mydata[name.row.nos,]

all(seg.f.mydata$the.v == seg.v.mydata$Name)

seg.v.mydata.short <- seg.v.mydata[,c("No.of.Batches", "End", "Name", "SAAV.PSM.sum", "PSM.ratio")]

min(seg.v.mydata.short$PSM.ratio)
max(seg.v.mydata.short$PSM.ratio)
seg.v.mydata.short[1,6] <- -0.3
seg.v.mydata.short[1015,6] <- max(seg.v.mydata.short$PSM.ratio)
seg.v.mydata.short[2:1014,6] <- 0.5

min(seg.v.mydata.short$SAAV.PSM.sum)
max(seg.v.mydata.short$SAAV.PSM.sum)
seg.v.mydata.short[1,7] <- min(seg.v.mydata.short$SAAV.PSM.sum)
seg.v.mydata.short[1015,7] <- max(seg.v.mydata.short$SAAV.PSM.sum)
seg.v.mydata.short[2:1014,7] <- log10(2)


# color the batch segments
colors.start <- RColorBrewer::brewer.pal(5, "Greens")
colors <- vector()
colors[1:5] <- colors.start[2]
colors[6:10] <- colors.start[3]
colors[11:15] <- colors.start[5]

# color the PSM sum
colors.start.PSM.sum <- RColorBrewer::brewer.pal(11, "Spectral")

colors.PSM <- vector()
for (i in 1:nrow(seg.v.mydata)) {
  if (seg.v.mydata[i,"End"] ==1) { 
    #colors.PSM[i] <- "grey"
    next
  }
  if (seg.v.mydata[i,"SAAV.PSM.sum"] > log10(20) ) { colors.PSM[i] <- colors.start.PSM.sum[1]
  } else if (seg.v.mydata[i,"SAAV.PSM.sum"] > log10(10) ) { colors.PSM[i] <- colors.start.PSM.sum[3]
  } else if (seg.v.mydata[i,"SAAV.PSM.sum"] > log10(5)) { colors.PSM[i] <- colors.start.PSM.sum[4]
  } else if (seg.v.mydata[i,"SAAV.PSM.sum"] > log10(2)) { colors.PSM[i] <- colors.start.PSM.sum[5] 
  } else { colors.PSM[i] <- colors.start.PSM.sum[6] }
} 

# color the PSM ratios
colors.start.ratio <- RColorBrewer::brewer.pal(9, "Set1")


colors.ratio <- vector()
for (i in 1:nrow(seg.v.mydata)) {
  if (seg.v.mydata[i,"End"] ==1) { 
    #colors.ratio[i] <- "grey"
    next
  }
  if (seg.v.mydata[i,"PSM.ratio"] > 0.50 ) { colors.ratio[i] <- colors.start.ratio[1]
  } else if (seg.v.mydata[i,"PSM.ratio"] > 0.30 ) { colors.ratio[i] <- colors.start.ratio[5]
  } else if (seg.v.mydata[i,"PSM.ratio"] > 0.002 ) { colors.ratio[i] <- colors.start.ratio[2]
  } else { colors.ratio[i] <- "grey" }
} 

db <- segAnglePo(seg.f.mydata, seg=seg.name, angle.start = 90, angle.end = 360)
par(mar=c(2,2,2,2))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")

circos(R=400, cir=db, type="chr", col=colors, print.chr.lab=F, W=75, scale = F)

circos(R=300, cir=db, W=100, mapping=seg.v.mydata, col.v=grep("SAAV.PSM.sum", colnames(seg.v.mydata)), 
       type="h", B=T , col=na.omit(colors.PSM), lwd=0.1, scale=F)

circos(R=180, cir=db, W=100, mapping=seg.v.mydata, col.v=grep("PSM.ratio", colnames(seg.v.mydata)), 
       type="s", B=T , col=colors.ratio, lwd=1, scale=F, cex=0.6)

circos(R=180, cir=db, W=100, mapping=seg.v.mydata.short, col.v=grep("V6", colnames(seg.v.mydata.short)), 
       type="lh", B=F , col="black", lwd=1, scale=F)

circos(R=140, cir=db, W=20, mapping=seg.v.mydata, col.v=grep("PeptideAtlas.Found", colnames(seg.v.mydata)), 
       type="h", col = "#A6761D",lwd=0.05, col.bar=TRUE, cluster=F)

circos(R=120, cir=db, W=20, mapping=seg.v.mydata, col.v=grep("cs.IDs", colnames(seg.v.mydata)), 
       type="h", col = "#6A3D9A",lwd=0.05, col.bar=TRUE, cluster=F)


all(seg.f.mydata$the.v == seg.v.mydata$Name)
all(seg.f.mydata$seg.End == seg.v.mydata$End)



#### Suppl.Figure for correlation analyis

cor.test(as.numeric(Verified.SAAVs$PSM.ratio), 
         as.numeric(Verified.SAAVs$Alt.AF.EU), 
         method = "spearman", use="pairwise.complete.obs")

library(ggplot2)

ggplot(Verified.SAAVs, aes(x=as.numeric(Alt.AF.EU), y=as.numeric(PSM.ratio))) +
  geom_point() + theme_bw() + ylab("PSMr") + xlab("AAF")#+ 
  #annotate("text", x = 0.3, y = 0.98,
           #label = "Spearman's rank correlation coefficient = 0.75\np-value < 2.2e-16") 



