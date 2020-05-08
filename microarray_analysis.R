#Chapter 6 - Micro Array Analysis
#Example1
#Membangun dataset MLL.B dari package ALLMLL
library(affy)
data(MLL.B, package = "ALLMLL")
MLL.B
#Mencetak struktur objek
slotNames(MLL.B)
#Mencetak jumlah baris dan kolom nilai ekspresi MLL.B
dim(exprs(MLL.B))
#Mengekstraksi anotasi
annotation(MLL.B)
#mencetak 10 nama pertama dari probe yang digunakan
probeNames(MLL.B)[1:10]
#Mencetak nilai PM
pm(MLL.B, "200000_s_at")[1:4,1:3]
#Membuat grafik variabilitas data probe yang akan didapat - Plot nilai intensitas
matplot(pm(MLL.B,"200000_s_at"),type="l", xlab = "Probe No.",ylab="PM Probe Intensity")
#Membuat grafik kepadatan data MLL.B
hist(MLL.B)
matplot(MLL.B, pairs = TRUE, plot.method="smoothScatter") #ndabisa
image(MLL.B)


#Preprocessing Methods
#Mencetak latarbelakang yang tersedia
bgcorrect.methods()
#Mengoreksi nilai PM - menghapus nilai MM
pmcorrect.methods()
#Normalisasai Data
normalize.methods(MLL.B)
#Pengumpulan nilai intensitas probe menjadi ekspresi gen
express.summary.stat.methods()

#Example 1
#Menggabungkan koreksi latar belakang RMA dengan normalisasi konstan 
eset <- expresso(MLL.B, bgcorrect.method = "rma", normalize.method = "constant", pmcorrect.method = "pmonly", summary.method = "avgdiff")
#Example 2
#Menggabungkan koreksi latar belakang konvolusi, normalisasi kuantil, dan peringkasan berdasarkan model multi-array sesuai dengan cara yang kuat oleh algoritma polian median
library(affy)
data(MLL.B, package = "ALLMLL")
eset3 <- rma(MLL.B)
boxplot(data.frame(exprs(eset3)))

#Example 3
#Memperoleh gambaran umum tentang jumlah pasien yang berada dalam fase penyakit tertentu
#BiocManager::install("ALL")
library(ALL)
data(ALL, package = "ALL")
slotNames(ALL)
row.names(exprs(ALL))[1:10]
feno <- pData(ALL)

#menghitung kolom mad dan median - mengurangi median dari setiap entri kolom
#Membagi setiap entri kolom oleh MAD
ALL1pp <- ALL1 <- ALL[,ALL$mol == "ALL1/AF4"]
mads <- apply(exprs(ALL1), 2, mad)
meds <- apply(exprs(ALL1), 2, median)
dat <- sweep(exprs(ALL1), 2, meds)
exprs(ALL1pp) <- sweep(dat, 2, mads, FUN="/")

#Gene Filtering
#Menghitung koefisien variasi per gen untuk data ALL1pp
#Example1
cvval <- apply(exprs(ALL1pp),1,function(x){sd(x)/abs(mean(x))})
sum(cvval<0.2)

#Example 2
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("genefilter")
#Penggabungan beberapa filter
library("genefilter")
f1 <- function(x)(IQR(x)>0.5)
f2 <- pOverA(.25, log2(100))
f3 <- function(x) (median(2^x) > 300)
f4 <- function(x) (shapiro.test(x)$p.value > 0.05)
f5 <- function(x) (sd(x)/abs(mean(x))<0.1)
f6 <- function(x) (sqrt(10)* abs(mean(x))/sd(x) > qt(0.975,9))
ff <- filterfun(f1,f2,f3,f4,f5,f6)
library("ALL"); data(ALL)
selected <- genefilter(exprs(ALL[,ALL$BT=="B"]), ff)
sum(selected)

#Example 3
#Membuat faktor logis patientB
#menunjukkan pasien dengan B-cell ALL (TRUE) dan T-cell ALL (FALSE)
library("genefilter");library("ALL"); data(ALL)
patientB <- factor(ALL$BT %in% c("B","B1","B2","B3","B4"))
f1 <- function(x) (shapiro.test(x)$p.value > 0.05)
f2 <- function(x) (t.test(x ~ patientB)$p.value < 0.05)
sel1 <- genefilter(exprs(ALL[,patientB==TRUE]), filterfun(f1))
sel2 <- genefilter(exprs(ALL[,patientB==FALSE]), filterfun(f1))
sel3 <- genefilter(exprs(ALL), filterfun(f2))
selected <- sel1 & sel2 & sel3
ALLs <- ALL[selected,]
sum(ALLs, )

#Membuat diagram venn
#BiocManager::install("limma")
library(limma)
x <- matrix(as.integer(c(sel1,sel2,sel3)),ncol = 3,byrow=FALSE)
colnames(x) <- c("sel1","sel2","sel3")
vc <- vennCounts(x, include="both")
vennDiagram(vc)

#Applications of Linear Models
#Menganalisis variansi
library("ALL"); library("limma");
data(ALL, package = "ALL")
allB <- ALL[,which(ALL$BT %in% c("B","B1","B2"))] 
design.ma <- model.matrix(~ 0 + factor(allB$BT))
colnames(design.ma) <- c("B","B1","B2")
fit <- lmFit(allB, design.ma)
fit <- eBayes(fit) 
toptab <- topTable(fit, coef=2,5,adjust.method="fdr") 
print(toptab[,1:5],digits=4) 

#Mencetak hipotesis tertenti
cont.ma <- makeContrasts(B-B1,B1-B2, levels=factor(allB$BT))
cont.ma
#Membuat kontras matriks
fit1 <- contrasts.fit(fit, cont.ma)
fit1 <- eBayes(fit1)
toptabcon <- topTable(fit, coef=2,5,adjust.method="fdr")
print(toptabcon[,1:5],digits=4)
toptabcon <- topTable(fit1, coef=2,5,adjust.method="fdr") 
print(toptabcon[,1:5],digits=4)

#Example 2
#Menggabungkan output khas
#BiocManager::install("annaffy")
#BiocManager::install("hgu95av2.db")
library("annaffy");library("hgu95av2.db")
saveHTML(anntable, "ALLB123.html", title = "B-cell 012 ALL")
aafTableAnn()

#Example 3
#Merangkum hasil dalam tabel HTML berdasarkan nilai-p dari mis. analisis varians (ANOVA)
#BiocManager::install("multtest")
library("multtest"); library("annaffy"); library("hgu95av2.db") 
library("ALL"); data(ALL, package = "ALL") 
ALLB <- ALL[,which(ALL$BT %in% c("B","B1","B2"))] 
panova <- apply(exprs(ALLB), 1, function(x) anova(lm(x ~ ALLB$BT))$Pr[1]) 
genenames <- featureNames(ALLB)[panova<0.000001]
atab <- aafTableAnn(genenames, "hgu95av2.db", aaf.handler()[c(1:3,8:9,11:13)])
saveHTML(atab, file="ANOVAonB-cellGroups.html")

#Example 4
#Menganalisis public available data
#memilih gen dengan efek interaksi yang signifikan
#BiocManager::install("GEOquery")
library(GEOquery); library(limma); library(hgu95av2.db);library(annaffy)
gds <- getGEO("GDS1365")
eset <- GDS2eSet(gds,do.log2=T)
prot <- pData(eset)$protocol
time <- pData(eset)$time
pvalt<- apply(exprs(eset)[1:12625,], 1,
  function(x) anova(lm~prot*time)) $Pr[1:3])
pvalt <- data.frame(t(pvalt))
colnames(pvalt) <- c("meffprot","mefftime","interaction")
genenames <- featureNames(eset)[pvalt$meffprot< 0.01 &
  pvalt$mefftime < 0.01 & pvalt$interaction < 0.01] 
atab <- aafTableAnn(genenames,"hgu95av2.db",aaf.handler()[c(1:3,8:9,11:13)]) 
saveHTML(atab, file="Two-way ANOVA protocol by time.html")

#Mencari package anotasi
library("ALL"); data(ALL)
annotation(ALL)
# the annotation package we need is hgu95av2.db
library(hgu95av2.db) 
ls("package:hgu95av2.db") 
#Membuat konten environment
ChrNrOfProbe <- as.list(hgu95av2CHR)
ChrNrOfProbe[1]

#Mencari environment berdasarkan nama
get("1389_at", env = hgu95av2ACCNUM)
get("1389_at", env = hgu95av2ENTREZID) 
get("1389_at", env = hgu95av2SYMBOL) 
get("1389_at", env = hgu95av2GENENAME) 
get("1389_at", env = hgu95av2SUMFUNC)
get("1389_at", env = hgu95av2UNIGENE)
#Mencari nukleotida database
library(annotate)
genbank("J03779",disp="browser")
genbank(1430782,disp="data",type="uid")
get("1389_at", env = hgu95av2CHRLOC)
#Lokasi cytoband
get("1389_at", env = hgu95av2MAP)
ll1<-GOENTREZID2GO[["4121"]]

#Pencarian Literatur menggunakan anotasi
library(hgu95av2.db);library(annotate); library(ALL); data(ALL) 
pmid <- get("1389_at",env=hgu95av2PMID) 
pubmed(pmid,disp="browser")
#Untuk melihat informasi penulis
absts <- pm.getabst("1389_at", "hgu95av2")
pm.titles(absts)
#Pencarian ekspresi regular
ne <- pm.abstGrep("neutral endopeptidase",absts[[1]])
# pmAbst2HTML(absts[[1]],filename="pmon1389_at.html")

#Mencari nomor GO dan evidence - 6.7
#extract a list from the annotation ﬁles hgu95av2GO
go1389 <- get("1389_at", env = hgu95av2GO)
idl <- lapply(go1389,function(x) x$GOID) 
idl[[1]]
#memilih GO numbers yang berelasi dengan proses biologi
library(annotate)
getOntology(go1389,"BP")
getEvidence#go1389)
#Membuat list
go1389TAS <- subset(go1389,getEvidence(go1389)=="TAS")
#Mengekstrak informasi dari list
sapply(go1389TAS,function(x) x$GOID)
sapply(go1389TAS,function(x) x$Evidence)
sapply(go1389TAS,function(x) x$Ontology)

#6.8 Go Parents and children
#Example 1
#Mengumpulkan Informasi Go
GOMFPARENTS$"GO:0003700"
GOMFCHILDREN$"GO:0003700" 
#Mengumpulkan ontologi, orang tua, dan pengidentifikasi anak dalam vektor
go1389 <- get("1389_at", env = hgu95av2GO)
gonr <- getOntology(go1389, "BP")
gP <- getGOParents(gonr)
gC <- getGOChildren(gonr) 
gPC <- c(gonr,gP,gC)
pa <- sapply(gP,function(x) x$Parents) 
ch <- sapply(gC,function(x) x$Children)
gonrc <- c(gonr,unlist(pa),unlist(ch))
          
#Example 2
#Seleksi Probe oleh GO
#BiocManager::install("GO")
library(GO); library(annotate); library("ALL"); data(ALL)
go1389 <- get("1389_at", env = hgu95av2GO)
gonr <- getOntology(go1389, "BP")
gP <- getGOParents(gonr)
pa <- sapply(gP,function(x) x$Parents)
probes <- mget(pa,hgu95av2GO2ALLPROBES)
probeNames <- unlist(probes)
ALLpr <- ALL[probeNames,]
dim(exprs(ALLpr)) 
         
#6.9 Filtering berdasarkan istilah biologis
#Example 1 ”transcriptional repression
library("GO"); library("annotate"); library("hgu95av2.db")
GOTerm2Tag <- function(term) { 
  GTL <- eapply(GOTERM, function(x) {grep(term, x@Term, value=TRUE)}) 
  Gl <- sapply(GTL, length) 
  names(GTL[Gl>0])
}
GOTerm2Tag("transcriptional repressor")

#Menyelidiki data ALLs
tran1 <- hgu95av2GO2ALLPROBES$"GO:0016564"
tran2 <- hgu95av2GO2ALLPROBES$"GO:0016566"
tran3 <- hgu95av2GO2ALLPROBES$"GO:0017053"
tran <- c(tran1,tran2,tran3)
inboth <- tran %in% row.names(exprs(ALLs))
ALLtran <- ALLs[tran[inboth],]

#6.10 Signifikansi per kromosom
library("ALL"); data(ALL); library("hgu95av2.db")
rawp <- apply(exprs(ALL), 1, function(x) t.test(x ~ ALL$remission)$p.value)
xx <- as.list(hgu95av2CHR)
AffimIDChr <- names(xx[xx=="19"])
names(rawp) <- featureNames(ALL)
f <- matrix(NA,2,2)
f[1,1] <- sum(rawp[AffimIDChr]<0.05); f[1,2] <- sum(rawp[AffimIDChr]>0.05)
f[2,1] <- sum(rawp<0.05) - f[1,1] ; f[2,2] <- sum(rawp>0.05) - f[1,2]
print(f)
fisher.test(f)
chisq.test(f)
