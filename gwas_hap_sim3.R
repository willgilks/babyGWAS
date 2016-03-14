## R code to test two traits for haploid genotype association

## libraries ########
require(broom)
require(tidyr)
require(parallel)
require(dplyr)
require(ggplot2)
# require(ggthemes)
require(cowplot)
# require(gridExtra)
# require(ggExtra)

### SIMULATION ##########
## simulation variables
nflies = 100
nsnps = 100
t1setMean = 90
t1setSd = 30
t2setMean = 80
t2setSd = 40

## simulate genotypes: snps by row, flies by column.
## freq distr could be more biologically accurate.
## need to add in NA for failed genotypes.
g = matrix(sample(0:1,nsnps*nflies, replace=TRUE),nsnps,nflies)
row.names(g) = paste("snp_", 1:nrow(g), sep="")
colnames (g) = paste("fly_", 1:ncol(g), sep="")

## assess missingness and call rates
# fly.means = colMeans(g, na.rm = TRUE)
# snp.means = rowMeans(g,na.rm = TRUE)


## simulate snp map with chromosome positions.
bp = 10 * seq(1:nsnps)
chr = c(rep("2L", (length(bp))/2),rep("2R", (length(bp))/2))
m = data.frame(chr, bp)
row.names(m)= paste("snp_", 1:nsnps , sep="")
head(m)
rm(bp, chr)


## simulate trait values
t1=rnorm (nflies, t1setMean, t1setSd)
t2=rnorm (nflies, t2setMean, t2setSd)

## standardise trait variables
svar = function (trait.vals) {
  mtv = mean (trait.vals)
  stv = (trait.vals - mtv) / (sd (trait.vals))
  return (stv) }

## standardise and bind trait variables
traitvals = cbind (svar(t1), svar(t2))
typeof(traitvals)
dim(traitvals)
head(traitvals)

## trait anatognism per genome
calant = function (t1, t2) {
  a = (t1 - t2) / ( 2*sin (pi/4) )
  return(a) }


## trait concordance per genome
calcon = function (t1, t2) {
  c = (t1 + t2) / ( 2*sin (pi/4) )
  return(c) }



## ASSOCIATION TEST
hwas = function (genos, phenos, map, test_type) {
  ## genos is a matrix with n columns.
  ## phenos is a matrix with n rowms.
  ## need to check merge map params and input.
  ## add data input warnings.
  # final.res$fdr = p.adjust(final.res$p.value, method = "fdr", n=length(final.res$p.value))
  ## (use 'parallel' for >100k snps.)
  r=apply(genos, 1, function(snp) test_type (phenos ~ snp))
          t=do.call(rbind, lapply(r, tidy))
          rm(r)
          f=merge(t, map, by=0, all=TRUE)
          rm(t)
          return(f) }


## do the tests.
asst1 = hwas(g, traitvals[,1], m, t.test)
asst2 = hwas(g, traitvals[,2], m, t.test)
assma = hwas(g, traitvals, m, manova)


traitvals = cbind (svar(t1), svar(t2))
testvals = tbl_df(data.frame(traitvals))

testvals$geno1 = g[which.min(assma$p.value),]

head (testvals)
# testvals = rename(testvals, X1="10", X2="20")
testvals$diff=testvals$X1-testvals$X2
hist(testvals$diff)

ggplot(testvals, aes(x = 0, y=X1)) +
  geom_segment(aes(xend = 1, yend = X2, 
                   colour=as.factor(geno1))) +
  theme(panel.background=element_rect(fill="black")) +
  scale_colour_manual(values=c("light blue", "green"))


ggplot(testvals, aes(x = 0, y=0)) +
  geom_segment(aes(xend = 1, yend = diff, 
                   colour=as.factor(geno1))) +
  theme(panel.background=element_rect(fill="black")) +
  scale_colour_manual(values=c("light blue", "green"))



# calculate average change per genome.
ggplot (testvals, aes(diff, fill=as.factor(geno1))) +
  geom_histogram()


## extract out mean trait values for each genotype from ttests.
## use these as line start-end coordinates on plot.
## this move is dodgy. Depends on snp sequence being equally matched.
geno_means = tbl_df(data.frame(snp=assma$Row.names, 
               pman=assma$p.value, 
               t1e1=asst1$estimate1, 
               t1e2=asst1$estimate2, 
               t2e1=asst2$estimate1, 
               t2e2=asst2$estimate2))
head(geno_means)

geno_means$log10p = -(log10(geno_means$pman))

topSnpMeans = filter(geno_means, log10p > 3)
head(topSnpMeans)


### multi-line plot showing effect size direction, and signficiance.
### longer lines, genotype greater effect size.
### brighter lines have lower p-value.
# png(filename = "~/Desktop/urchin1000.png")
manova.plot = ggplot(geno_means, aes(t1e1, t2e1)) +
  geom_segment(aes(xend = t1e2, yend = t2e2, colour=log10p)) +
  scale_colour_gradientn(colours = rainbow(7)) +
  labs(x="trait 1", y="trait 2") + theme_classic() + 
  theme(panel.background=element_rect(fill="black"))
manova.plot
dev.off()
###########
head(topSnpMeans)

png(filename = "~/Desktop/urchin1000.png")
# manova.plot2 = 
  ggplot(topSnpMeans, aes(t1e1, t2e1)) +
    geom_segment(aes(xend = t1e2, yend = t2e2)) +
    labs(x="trait 1", y="trait 2") + theme_classic()
# manova.plot2
dev.off()


## fdr ##########
# pvalue<-c(.03, .002,.002,.93) 
# sorted.pvalue = sort(pvalue)
# j.alpha = (1:length(pvalue))*(.05/length(pvalue)) 
# dif = sorted.pvalue-j.alpha
# neg.dif = dif[dif<0]
# pos.dif = neg.dif[length(neg.dif)]
# index<-dif==pos.dif
# p.cutoff<-sorted.pvalue[index]
# p.sig<-pvalue[pvalue<=p.cutoff]
# p.sig
########

## make data for qq plot try ggplot stat for this.
  obs = -log10 (sort(res$p.value))
  exp = -log10 (1:length (obs)/length (obs))
  r = tbl_df (data.frame(obs, exp))
head(r)
qplot(obs, exp)


# Make QQ DATA (maybe use list or for loop)  ########### 
p.trait1.obs =  -log10 (sort (assoc.trait1$p.value ))
p.trait1.exp =  -log10 ( 1 : length (p.trait1.obs) / length (p.trait1.obs))
p.trait2.obs =  -log10 (sort (assoc.trait2$p.value ))
p.trait2.exp =  -log10 ( 1 : length (p.trait2.obs) / length (p.trait2.obs))
p.trait.ant.obs =  -log10 (sort (assoc.antagonism$p.value ))
p.trait.ant.exp =  -log10 ( 1 : length (p.trait.ant.obs) / length (p.trait.ant.obs))
p.maxlim = 0.5 + max(p.trait1.obs, p.trait1.exp, p.trait2.obs, p.trait2.exp, p.trait)



## fdr ######
# pval.ctd2 = p.adjust(pval2, method="fdr")
# snps.sel.univari2 = which(pval.ctd2 < 0.05)
# n.snps.sel.univar2 = length(snps.sel.univar2)
# n.snps.sel.univar2
# snps.sel.univar2

## references ######
# http://adegenet.r-forge.r-project.org/files/Leuven2014/practical-GWAS.pdf

