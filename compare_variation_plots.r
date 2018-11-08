library(data.table)

proportions_1331 <- read.csv("p1331_variant_comparison_oct30.txt", sep=",", as.is=TRUE, header=TRUE)
proportions_1332 <- read.csv("p1332_variant_comparison_oct30.txt", sep=",", as.is=TRUE, header=TRUE)
proportions_1348 <- read.csv("p1348_variant_comparison_oct30.txt", sep=",", as.is=TRUE, header=TRUE)
p1331<-data.table(proportions_1331)
p1332<-data.table(proportions_1332)
p1348<-data.table(proportions_1348)


t_col <- function(color, percent = 50, name = NULL) {
#	  color = color name
#	percent = % transparency
#	   name = an optional name for the color
## Get RGB values for named color
  rgb.val <- col2rgb(color)
## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
## Save the color
  invisible(t.col)
}

# Compare Illumina RCA and nanopore:
pdf("Fig6A_diversity_IlluminaRCA_Nanopore_oct30.pdf", height=6, width=10)
par(mfrow=c(1,2))
plot(x=p1348$Illumina1, y=p1348$Nanopore,xlim=c(0,1), ylim=c(0,1), xlab="Illumina RCA proportion non-consensus", ylab="Nanopore proportion non-consensus", pch=19, col=t_col("cornflowerblue",70), cex=1.4)
points(c(0.0, 1.0), c(0.0,1.0), type="l", lty="dashed", col=t_col("red3", 70), lwd=2)
points(x=p1332$Illumina1,y=p1332$Nanopore, pch=19, col=t_col("gray4",70), cex=1.4)
points(x=p1331$Illumina1,y=p1331$Nanopore, pch=19, col=t_col("orange3",70), cex=1.4)

# calculate number of reads in the 1% box
in_box_1331<-sum(p1331$Nanopore < 0.01 & p1331$Illumina1 < 0.01, na.rm=T)
in_box_1332<-sum(p1332$Nanopore < 0.01 & p1332$Illumina1 < 0.01, na.rm=T)
in_box_1348<-sum(p1348$Nanopore < 0.01 & p1348$Illumina1 < 0.01, na.rm=T)
total_1331<-sum(p1331$Nanopore <= 1 & p1331$Illumina1 <= 1, na.rm=T)
total_1332<-sum(p1332$Nanopore <= 1 & p1332$Illumina1 <= 1, na.rm=T)
total_1348<-sum(p1348$Nanopore <= 1 & p1348$Illumina1 <= 1, na.rm=T)
(in_box_1331 + in_box_1332 + in_box_1348)/(total_1331 + total_1332 + total_1348)

plot(x=p1348$Illumina1_nodels, y=p1348$Nanopore_nodels,xlim=c(0,1), ylim=c(0,1), xlab="Illumina RCA proportion non-consensus base", ylab="Nanopore proportion non-consensus base", pch=19, col=t_col("cornflowerblue",70), cex=1.4)
points(c(0.0, 1.0), c(0.0,1.0), type="l", lty="dashed", col=t_col("red3", 70), lwd=2)
points(x=p1332$Illumina1_nodels,y=p1332$Nanopore_nodels, pch=19, col=t_col("gray4",70), cex=1.4)
points(x=p1331$Illumina1_nodels,y=p1331$Nanopore_nodels, pch=19, col=t_col("orange3",70), cex=1.4)
dev.off()

# now plot Illumina RCA vs IlluminaCL
proportions <- read.csv("~/Dropbox/Nanopore/HBV/variation_coverage_above_100.csv", sep=",", as.is=TRUE, header=TRUE)
pdf("Fig2C_CL_vs_RCA_variation_all_samples.pdf")
plot(proportions[,5], proportions[,6], col=t_col("cornflowerblue", 70), xlim=c(0,0.5), ylim=c(0,0.5), pch=19, xlab="", ylab="", cex=1.4, cex.axis=1.3)
title(main="Per-site variation in RCA vs CL only samples", xlab="Illumina CL only", ylab="Illumina CL + RCA", cex.lab=1.3, cex.main=1.3)
points(proportions[,3], proportions[,4], col=t_col("gray4",70), pch=19, cex=1.4)
points(proportions[,1], proportions[,2], col=t_col("orange3",70), pch=19, cex=1.4)
points(c(0.0, 1.0), c(0.0,1.0), type="l", lty="dashed", col=t_col("red2"))
dev.off()


