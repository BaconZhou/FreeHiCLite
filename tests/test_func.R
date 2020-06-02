library(FreeHiCLite)
library(pryr)

#dat1 <- hicData("/Users/pgzhou/Documents/Juicebox/data/inter.hic", "All", "All", "BP", 50000)
# dat2 <- hicData("/Users/pgzhou/Documents/Juicebox/data/inter.hic", "1", "2", "BP", 5000)


f = "DESCRIPTION"
normalizePath(f)
file1 <- "/Users/pgzhou/Documents/hicLib//data/rabbit_rep1_30.hic"
#file2 <- "/Users/pgzhou/Documents/Juicebox/data/inter.hic"
file3 <- "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic"


system.time({dat1 <- FreeHiCLite::readJuicer(file1, c("chr200", "chr2"), NULL, "FRAG", 50L)})
info <- FreeHiCLite::readJuicerInformation(file1)
str(info)

info <- FreeHiCLite::readJuicerInformation(file3)
str(info)


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
