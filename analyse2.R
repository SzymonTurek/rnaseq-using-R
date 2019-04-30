plot(log2(WT.R1.g.exp[as.character(TR.or[,1])]+0.5), log2(WT.R1.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot WT.g vs WT.t")
if (!requireNamespace("BiocManager"))install.packages("BiocManager")
BiocManager::install("Rsubread")

library(Rsubread)


buildindex(basename="genome.index", reference="/data/next-gen-sequencing/data/TAIR9.fa", indexSplit=TRUE, memory=4000)


subjunc(index="genome.index", readfile1="/data/next-gen-sequencing/data/OE_1_R1.fastq", output_file = "OE_1_R1.BAM", nthreads=4, input_format = "FASTQ", sortReadsByCoordinates=TRUE)

subjunc(index="genome.index", readfile1="/data/next-gen-sequencing/data/WT_R1.fastq", output_file = "WT_R1.BAM", nthreads=4, input_format = "FASTQ", sortReadsByCoordinates=TRUE)

fc <- featureCounts(c("/data/next-gen-sequencing/users/plabaj/OE_1_R1.BAM", "/data/next-gen-sequencing/users/plabaj/WT_R1_sorted.bam"), annot.ext="/data/next-gen-sequencing/data/TAIR9.gtf", isGTFAnnotationFile=TRUE)


colnames(fc$stat) <- c("Status","WT", "OE1")

pchs <- rep('+', dim(fc$counts)[1])
cols <- rep('black', dim(fc$counts)[1])
names(pchs) <- names(cols) <- rownames(fc$counts)
sel.gene <- "AT3G01150"
pchs[sel.gene] <- "x"
cols[sel.gene] <- "red"
 

plot(log2(fc$counts[,2]+0.5), log2(fc$counts[,1]+0.5),xlab="WT", ylab="OE_1", pch='+', col='black')
##smoothScatter(log2(fc$counts[,2]+0.5), log2(fc$counts[,1]+0.5),xlab="WT", ylab="OE_1")
dim(fc$counts)
plot(log2(fc$counts[,2]+0.5), log2(fc$counts[,1]+0.5),xlab="WT", ylab="OE_1", pch=pchs, col=cols, main = "Scatter plot WT vs OE_1")

dev.copy2pdf(file="Fig1.pdf")


TR.or <- read.table(file='/data/next-gen-sequencing/data/TAIR9.tr')
TR.gene <- TR.or
TR.gene[,1] <- substr(TR.gene[,2],1,9)

write.table(TR.gene, file="TAIR_gene.tr", quote=FALSE, row.names=FALSE, col.names=FALSE)

##BitSeq expression estimates analisis

OE.R1.t <- read.table(file='/data/next-gen-sequencing/res/OE_1_R1.counts')
OE.R2.t <- read.table(file='/data/next-gen-sequencing/res/OE_1_R2.counts')
WT.R1.t <- read.table(file='/data/next-gen-sequencing/res/WT_R1.counts')
WT.R2.t <- read.table(file='/data/next-gen-sequencing/res/WT_R2.counts')
##genes
OE.R1.g <- read.table(file='/data/next-gen-sequencing/res/OE_1_R1_gene.counts') 
OE.R2.g <- read.table(file='/data/next-gen-sequencing/res/OE_1_R2_gene.counts')
WT.R1.g <- read.table(file='/data/next-gen-sequencing/res/WT_R1_gene.counts')
WT.R2.g <- read.table(file='/data/next-gen-sequencing/res/WT_R2_gene.counts')

##mean value of expreseion for each transcript/gene

OE.R1.t.exp <- apply(as.matrix(OE.R1.t),1,mean)
OE.R2.t.exp <- apply(as.matrix(OE.R2.t),1,mean)
WT.R1.t.exp <- apply(as.matrix(WT.R1.t),1,mean)
WT.R2.t.exp <- apply(as.matrix(WT.R2.t),1,mean)

OE.R1.g.exp <- apply(as.matrix(OE.R1.g),1,mean)
OE.R2.g.exp <- apply(as.matrix(OE.R2.g),1,mean)
WT.R1.g.exp <- apply(as.matrix(WT.R1.g),1,mean)
WT.R2.g.exp <- apply(as.matrix(WT.R2.g),1,mean)

plot(density(as.matrix(OE.R1.t.exp)[1,]))

TR.or <- read.table(file='/data/next-gen-sequencing/users/ngseq04/TAIR_gene.tr')
gene.names <- unique(as.character(TR.or[,1]))
##on the gene level
pchs <- rep('+', length(gene.names))
cols <- rep('black', length(gene.names))
names(pchs) <- names(cols) <- gene.names
sel.gene <- "AT3G01150"

pchs[sel.gene] <- 'x'
cols[sel.gene] <- "cyan"

#plot(log2(fc$counts[,2]+0.5), log2(fc$counts[,1]+0.5),xlab="WT", ylab="OE_1", pch=pchs, col=cols, main = "Scatter plot WT vs OE_1")

plot(log2(WT.R1.g.exp+0.5), log2(OE.R1.g.exp+0.5), xlab = "WT", ylab="R2", pch=pchs, col=cols, main = "Scatter plot WT vs R2")

plot(log2(OE.R1.g.exp+0.5), log2(OE.R2.g.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")
##########

plot(log2(WT.R1.t.exp+0.5), log2(OE.R2.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")
# on the transcript level
pchs <- rep('.', dim(TR.or)[1])
cols <- rep('black', dim(TR.or)[1])
names(pchs) <- names(cols) <- as.character(TR.or[,2])
sel.gene <- "AT3G01150"
sel.trans <- as.character(TR.or[TR.or[,1]==sel.gene,2])

pchs[sel.trans] <- 'x'
cols[sel.trans] <- 'red'

plot(log2(OE.R1.t.exp+0.5), log2(OE.R2.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")
plot(log2(WT.R1.t.exp+0.5), log2(OE.R2.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")

plot(log2(WT.R1.t.exp+0.5), log2(OE.R1.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")
abline(a=0,b=1,col='red')

plot(log2(WT.R1.t.exp+0.5), log2(OE.R2.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")

cd plot(log2(WT.R1.t.exp+0.5), log2(OE.R1.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")

##gene vs transcript

chs <- rep('.', dim(TR.or)[1])
cols <- rep('black', dim(TR.or)[1])
names(pchs) <- names(cols) <- as.character(TR.or[,2])
sel.gene <- "AT3G01150"
sel.trans <- as.character(TR.or[TR.or[,1]==sel.gene,2])

pchs[sel.trans] <- 'x'
cols[sel.trans] <- 'red'

head(TR.or)
names(OE.R1.g.exp) <- gene.names
plot(log2(OE.R1.g.exp[as.character(TR.or[,1])]+0.5), log2(OE.R1.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot R1 vs R2")

names(WT.R1.g.exp) <- gene.names


plot(log2(WT.R1.g.exp[as.character(TR.or[,1])]+0.5), log2(WT.R1.t.exp+0.5), xlab = "R1", ylab="R2", pch=pchs, col=cols, main = "Scatter plot WT.g vs WT.t")


##LIMMA,...
if (!requireNamespace("BiocManager", quietly = TRUE)) # wczytanie pakietów, limma, edgeR, DESeq2
    install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")
BiocManager::install("edgeR", version = "3.8")
BiocManager::install("DESeq2", version = "3.8")

##BitSeq expression estimates analisis
#wczytanie macierzy dla genu i transkryptu, nadekspresja oraz typ dziki
OE.R1.t <- read.table(file='/data/next-gen-sequencing/res/OE_1_R1.counts')
OE.R2.t <- read.table(file='/data/next-gen-sequencing/res/OE_1_R2.counts')
WT.R1.t <- read.table(file='/data/next-gen-sequencing/res/WT_R1.counts')
WT.R2.t <- read.table(file='/data/next-gen-sequencing/res/WT_R2.counts')
##genes
OE.R1.g <- read.table(file='/data/next-gen-sequencing/res/OE_1_R1_gene.counts') 
OE.R2.g <- read.table(file='/data/next-gen-sequencing/res/OE_1_R2_gene.counts')
WT.R1.g <- read.table(file='/data/next-gen-sequencing/res/WT_R1_gene.counts')
WT.R2.g <- read.table(file='/data/next-gen-sequencing/res/WT_R2_gene.counts')

##mean value of expreseion for each transcript/gene
#średnia wartość ekspresji dla transkyptu/genu
OE.R1.t.exp <- apply(as.matrix(OE.R1.t),1,mean)
OE.R2.t.exp <- apply(as.matrix(OE.R2.t),1,mean)
WT.R1.t.exp <- apply(as.matrix(WT.R1.t),1,mean)
WT.R2.t.exp <- apply(as.matrix(WT.R2.t),1,mean)

OE.R1.g.exp <- apply(as.matrix(OE.R1.g),1,mean)
OE.R2.g.exp <- apply(as.matrix(OE.R2.g),1,mean)
WT.R1.g.exp <- apply(as.matrix(WT.R1.g),1,mean)
WT.R2.g.exp <- apply(as.matrix(WT.R2.g),1,mean)
#wczytanie pliku TAIR_gene.tr, kolumna 1 dla genu, 2. dla transkryptu
TR.or <- read.table(file='/data/next-gen-sequencing/users/ngseq04/TAIR_gene.tr')
gene.names <- unique(as.character(TR.or[,1]))


names(OE.R1.g.exp) <-names(OE.R2.g.exp) <- names(WT.R1.g.exp) <-names(WT.R2.g.exp) <-gene.names#nazwanie kolejno wczytanych danych nazwami z pliku TAIR_gene.tr
genes.merged <- cbind(OE.R1.g.exp,OE.R2.g.exp,WT.R1.g.exp,WT.R2.g.exp)
colnames(genes.merged) <- c("OE.R1", "OE.R2", "WT.R1", "WT.R2")#nazwanie kolumn

samples <- substr(colnames(genes.merged),0,2) #wycina 2 znaki i mam "OE" "OE" "WT" "WT"
design <- data.frame(OEs=ifelse(samples=="OE",1,0), WTs=ifelse(samples=="WT",1,0))# tworzymy macierz na podstanie poprzednich nazw, obrazującą nadekspresję oraz brak nadekspresji

rownames(design) <- colnames(genes.merged)

cm <- makeContrasts(OEvsWT=OEs-WTs, levels=design)

dge <- DGEList(counts=genes.merged)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=FALSE)

dge$samples
summary(dge$counts)
v$targets
v%design
summary(v$E)

f.t <- lmFit(v, design)

cf <- contrasts.fit(f.t, cm)
fe <- eBayes(cf, proportion=0.01)

adj.method='BH'
limma.countsTMMvoom.genes <- topTable(fe, number=Inf, adjust.method=adj.method, sort.by='none')
sum(limma.countsTMMvoom.genes$adj.P.Val <0.05)
which(limma.countsTMMvoom.genes$adj.P.Val <0.05)
limma.countsTMMvoom.genes[which(limma.countsTMMvoom.genes$adj.P.Val <0.05),]



pchs <- rep('.', length(gene.names))
cols <- rep('black', length(gene.names))
names(pchs) <- names(cols) <- gene.names
sel.gene <- "AT3G01150"

pchs[sel.gene] <- 'x'
cols[sel.gene] <- "cyan"

plot(limma.countsTMMvoom.genes$AveExpr, limma.countsTMMvoom.genes$logFC, pch=pchs, col=cols, main = "GEN")
dev.copy2pdf(file="Fig1.pdf")


#dla transkryptu

TR.or <- read.table(file='/data/next-gen-sequencing/users/ngseq04/TAIR_gene.tr')
gene.names <- unique(as.character(TR.or[,2]))

names(OE.R1.t.exp) <-names(OE.R2.t.exp) <- names(WT.R1.t.exp) <-names(WT.R2.t.exp) <-gene.names
genes.merged <- cbind(OE.R1.t.exp,OE.R2.t.exp,WT.R1.t.exp,WT.R2.t.exp)
colnames(genes.merged) <- c("OE.R1", "OE.R2", "WT.R1", "WT.R2")

samples <- substr(colnames(genes.merged),0,2) 
design <- data.frame(OEs=ifelse(samples=="OE",1,0), WTs=ifelse(samples=="WT",1,0))# czy oe czy wt

rownames(design) <- colnames(genes.merged)

cm <- makeContrasts(OEvsWT=OEs-WTs, levels=design)

dge <- DGEList(counts=genes.merged)
dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=FALSE)

dge$samples
summary(dge$counts)
v$targets
v%design
summary(v$E)

f.t <- lmFit(v, design)

cf <- contrasts.fit(f.t, cm)
fe <- eBayes(cf, proportion=0.01)

adj.method='BH'
limma.countsTMMvoom.genes <- topTable(fe, number=Inf, adjust.method=adj.method, sort.by='none')
sum(limma.countsTMMvoom.genes$adj.P.Val <0.05)
which(limma.countsTMMvoom.genes$adj.P.Val <0.05)
limma.countsTMMvoom.genes[which(limma.countsTMMvoom.genes$adj.P.Val <0.05),]



pchs <- rep('.', length(gene.names))
cols <- rep('black', length(gene.names))
names(pchs) <- names(cols) <- gene.names
sel.gene <- "AT3G01150_ID4"

pchs[sel.gene] <- 'x'
cols[sel.gene] <- "cyan"

plot(limma.countsTMMvoom.genes$AveExpr, limma.countsTMMvoom.genes$logFC, pch=pchs, col=cols, main = "TRANSKRYPT")
dev.copy2pdf(file="transkrypt.pdf")

