library(FreeHiCLite)
library(pryr)

#dat1 <- hicData("/Users/pgzhou/Documents/Juicebox/data/inter.hic", "All", "All", "BP", 50000)
# dat2 <- hicData("/Users/pgzhou/Documents/Juicebox/data/inter.hic", "1", "2", "BP", 5000)


f = "DESCRIPTION"
normalizePath(f)
file1 <- "/Users/pgzhou/Documents/hicLib//data/rabbit_rep1_30.hic"
#file2 <- "/Users/pgzhou/Documents/Juicebox/data/inter.hic"
file3 <- "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic"


system.time({dat1 <- FreeHiCLite::readJuicer(file1, c("chr2", "chr2"), NULL, "FRAG", 50L)})
info <- FreeHiCLite::readJuicerInformation(file1)
frag <- FreeHiCLite::readJuicerFragmentSites(file1, c("chr1"))
str(info)

info <- FreeHiCLite::readJuicerInformation(file3)
str(info)

tmp <- test2 <- dat1$contact$`2_2`
df = convertJuicer(tmp, "chr2", "chr2")
df <- df[, c("chr1", "pos1", "chr2", "pos2", "score")]
df$chr1 <- sample(c("chr1", "chr2", "chr3"), NROW(df), replace = T)
df$chr2 <- sample(c("chr1", "chr2", "chr3"), NROW(df), replace = T)

tmp <- df[, c("pos1", "pos2", "score")]

pairs <- paste(df$chr1, df$chr2, sep="_")
res <- split(tmp, pairs)

p <- names(res)
a1 = sapply(p, FUN = function(x) {strsplit(x, '_')})
a2 =  sapply(res, NROW)
do.call('rbind', sapply(1:2, FUN = function(i) {rep(a1[[i]], a2[[i]])}))

ch <- do.call('rbind', res)

system.time({dat1 <- FreeHiCLite::readJuicer(file1, c("chr1", "chr2"), c("1_1", "2_2"), "FRAG", 50L)})
system.time({dat2 <- FreeHiCLite::readJuicer(file3, c("chr1"), "BP", c(500000L))})

#system.time({dat1 <- readJuicer(file1, c("chr1", "chr2"), "FRAG", 1L)})
test <- read.table('/Users/pgzhou/Documents/hicLib//data/out2.txt')
NROW(test)
test2 <- dat1$contact$chr1_chr2
NROW(test2)
all(test == test2)
str(dat1)
N = NROW(dat1)
spikeIn <- dat1[sort(sample(1:N, 0.1 * N)),]
spikeIn[,3] = spikeIn[,3] * 30

dat2 <- spikein(dat1, SpikeInSignal = spikeIn, 50000, smooth = TRUE)

dat <- list('1_1' = dat1)

simu <- hicDataSimuPure(dat, '1_1', 50000, sequenceDepth = 1000, countScale = 2, noiseRate = 0.5, neighborZeroRate = 0.5)


# dat2 <- hicDataHttp("https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic", "1", "1", "BP", 50000)

file <- "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic"

system.time({dat3 = hicDataSimu(fileName = "/Users/pgzhou/Documents/Juicebox/data/inter.hic",
                                chromosomes =  c("1", "2", "3"), unit = "BP", resolution = 5000, 
                                sequenceDepth = 1000, countScale = 2, noiseRate = 0.5, neighborZeroRate = 0.5)})

dat3 = hicDataSimu(fileName = "/Users/pgzhou/Documents/Juicebox/data/inter.hic",chromosomes =  c("1", "2", "3"), unit = "BP", resolution = 5000, 
                   sequenceDepth = 1000, countScale = 1, noiseRate = 0.0, neighborZeroRate = 0.5)
dat3[['1_1']]
str(dat3)


library(FreeHiCLite)

N <- 2000
maxBinX <- maxBinY <- 2000000L
maxCounts <- 10
binX <- sample(1:maxBinX, N, replace=TRUE)
binY <- sample(1:maxBinY, N, replace=TRUE)
counts <- sample(1:10, N, replace=TRUE)

seqDepth <- 20000
countScale <- 0
noiseRate <- 0.1
neighborZeroRate <- 0
resolution <- 5000

## Matrix version
## matrix layout as x, y, counts
contacts <- matrix(0, N, 3)
contacts[,1] <- binX
contacts[,2] <- binY
contacts[,3] <- counts

res <- FreeHiC(contacts, seqDepth, countScale, noiseRate, neighborZeroRate, resolution)
head(res)

## List version
N2 = 3000
binX2 <- sample(1:maxBinX, N2, replace=TRUE)
binY2 <- sample(1:maxBinY, N2, replace=TRUE)
counts2 <- sample(1:10, N2, replace=TRUE)

contacts2 <- matrix(0, N2, 3)
contacts2[,1] <- binX2
contacts2[,2] <- binY2
contacts2[,3] <- counts2

contactsMap <- list("1_1" = contacts, "2_2" = contacts2)
res <- FreeHiC(contactsMap, seqDepth, countScale, noiseRate, neighborZeroRate, resolution)
str(res)

## Dataframe version

chr1 <- c(rep("1", N), rep("2", N2))
chr2 <- c(rep("1", N), rep("2", N2))

contactsAll <- rbind(contacts, contacts2)

contactsDf <- data.frame(chr1 = chr1, x = contactsAll[,1], chr2 = chr2, y = contactsAll[,2], counts = contactsAll[,3])

res <- FreeHiC(contactsDf, seqDepth, countScale, noiseRate, neighborZeroRate, resolution)
str(res)

