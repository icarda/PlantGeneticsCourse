#' https://cran.r-project.org/web/packages/statgenGWAS/
#' https://biometris.github.io/statgenGWAS/articles/GWAS.html

library(statgenGWAS)


#' Markers can either be coded as character strings or as numerical values.
#' Common coding styles include the following options:
#' * AA, AB, BA, BB -> 0, 1, 1, 2 (A is reference allele, and B is alternative allele)
#' * CC, CT, TC, TT -> 0, 1, 1, 2 (e.g., ref. allele is C)
#' * C, Y, Y, T     -> 0, 1, 1, 2 (Nucleic Acid Notation: https://en.wikipedia.org/wiki/Nucleic_acid_notation)
#' * numerical code:   0, 1, 1, 2 (0 ref./major allele, 2 alt./minor allele, and 1 is heterozygous)  statgenGWAS
#' * numerical code:  -1, 0, 0, 1 (-1 ref./major allele, 1 alt./minor allele, and 0 is heterozygous) rrBLUP
#' 
#' Common Genotypic File Formats:
#' https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load
#' 
#' read the marker matrix

genotypic_data <- read.table('Raw_SNPS_80_Present_5_MAF.hmp.txt', 
                             header = TRUE, 
                             comment.char = '')


#' Marker map
#' The data.frame map is used to describe the physical positions of the markers 
#' on the chromosomes. The data consists of two columns, "chr" for the name or 
#' number of the chromosome and "pos" for the position of the marker on the 
#' chromosome. The position can be in basepair or in centimorgan. 
#' The names of the markers should be the row names of the data.frame.

marker_map <- genotypic_data[, 3:4]
colnames(marker_map) <- c('chr', 'pos')
rownames(marker_map) <- genotypic_data[, 1]


#' Marker matrix
#' 
#' The marker matrix is stored in the matrix marker within the gData object. 
#' It has the names of the markers in its column names and the genotypes in 
#' its row names. Markers can either be coded as character strings or as 
#' numerical values. 

marker_matrix <- genotypic_data[, -c(1:11)]
rownames(marker_matrix) <- genotypic_data[, 1]
marker_matrix <- t(marker_matrix)


#' Create gData object

gData <- createGData(geno = marker_matrix, map = marker_map)


#' Recoding and cleaning of markers
#' Marker data has to be numerical and without missing values in order to do 
#' GWAS analysis. This can be achieved using the codeMarkers function.


#' Important note, before performing any analysis, the marker matrix has to be 
#' converted to a numerical matrix. This can be done using this function.
#' 
#' - minor (alternative) allele numerical code is 2.
#' - major (reference) allele numerical code is 0.
#' - heterozygous allele numerical code is 1.
#' 
#' Remove duplicate SNPs from gData

gData <- codeMarkers(gData, 
                     naStrings = 'N', 
                     removeDuplicates = TRUE, 
                     impute = FALSE, 
                     verbose = TRUE)


#' * nMiss: A numerical value between 0 and 1. 
#' SNPs with a fraction of missing values higher than nMiss will be removed.

gData <- codeMarkers(gData, 
                     impute = FALSE, 
                     verbose = TRUE, 
                     nMiss = 0.1)


#' * nMissGeno: A numerical value between 0 and 1. 
#' Genotypes with a fraction of missing values higher than nMissGeno will be removed.

gData <- codeMarkers(gData, 
                     impute = FALSE, 
                     verbose = TRUE, 
                     nMissGeno = 0.2)


#' * MAF: A numerical value between 0 and 1. 
#' SNPs with a Minor Allele Frequency (MAF) below this value will be removed.

gData <- codeMarkers(gData, 
                     impute = FALSE, 
                     verbose = TRUE, 
                     MAF = 0.05)


#' Total number of remaining missing values and their ratio
#' Random imputation is OK when ratio of missing values is less than 5%

(d <- dim(gData$markers))
(n <- sum(is.na(gData$markers)))
(r <- n/(d[1]*d[2]))


#' To impute the missing values added above, codeMarkers has to be run with 
#' impute = TRUE. Then the type of imputation can be chosen by setting the 
#' parameter imputeType:
#
#' * fixed:  Impute all missing values by a single fixed value. 
#'           Use the parameter fixedValue to set this value.
#'           
#' * random: Impute missing values with a random value based on 
#'           the non-missing values for a SNP.
#'           
#' * beagle: Impute missing values using the independent beagle 
#'           software (Browning and Browning 2007).

gData <- codeMarkers(gData, 
                     impute = TRUE, 
                     imputeType = 'beagle', 
                     verbose = TRUE) 


#' Plot genetic map

plot(gData)


################################################################################


#' Phenotypic data, either directly from field trials or after summarizing can 
#' be stored in pheno in the gData object. Pheno can either be a single data.frame 
#' or a list of data.frames for storing data for different trials or different 
#' summarizations of the original data. The first column of all elements of 
#' pheno should be genotype and all the other columns should represent different 
#' traits. Storing additional variables should be done in covar. 

#' read the phenotypic data (first column should be "genotype")

phenotypic_data <- read.csv('Pheno_REML.csv')


#' Add phenotypic data to gData

gData <- createGData(gData, pheno = phenotypic_data)


#' Single trait GWAS
#'
#' Variance covariance matrix
#' There are two ways to compute the phenotypic variance covariance matrix used 
#' in the GWAS analysis. Either the EMMA algorithm (Kang et al. 2008) or the 
#' Newton-Raphson algorithm (Tunnicliffe 1989). Specify the method by setting 
#' the parameter remlAlgo to either "EMMA" or "NR". By default the EMMA 
#' algorithm is used.
#'
#' Significance thresholds
#' The threshold for selecting significant SNPs in a GWAS analysis is computed 
#' by default using Bonferroni correction, with an alpha of 0.05. The alpha can 
#' be modified setting the option alpha when calling runSingleTraitGwas. 

#' good Bonferroni correction approximation
(pvalue.thr <- 0.05 / ncol(gData$markers))
format(pvalue.thr, scientific = FALSE)

#' LOD score (logarithm of the odds)
(LOD.thr <- -log10(pvalue.thr))

#' reverse the conversion to get p-value threshold
10^(-LOD.thr)

#' exact Bonferroni correction formula
(pvalue.thr <- 1 - (1 - 0.05) ^ (1 / ncol(gData$markers)))
(LOD.thr <- -log10(pvalue.thr))


#' Two other threshold types can be used: a fixed threshold (thrType = "fixed") 
#' specifying the -log10(p) (LODThr) value of the threshold, or a threshold that 
#' defines the n SNPs with the highest log10(p) scores as significant SNPs. 
#' Set thrType = "small" together with nSnpLOD = n to do this.
#'
#' Kinship matrix
#' The kinship matrix describes the genetic relatedness between the different 
#' genotypes. This should be a square matrix with genotypes in both row and 
#' column names and a measure for the genetic relatedness in its cells. 
#' There are many ways of computing genetic relatedness some of which are 
#' included in this package:
#' * astle: using the covariance between the scaled SNP-scores (see e.g. equation (2.2) in Astle and Balding (2009))
#' * IBS: Identity by State (see e.g. equation (2.3) in Astle and Balding (2009))
#' * vanRaden: using the formula by VanRaden (2008)
#' * kin: user-defined.
#'
#' Minor Allele Frequency
#' By default all SNPs with a MAF lower than 0.01 are excluded from the analysis. 
#' This can be controlled by the parameter MAF.

GWAS <- runSingleTraitGwas(gData, 
                           traits  = 'ASC_Score', 
                           thrType = 'fixed', 
                           LODThr  = 3.5, 
                           kinshipMethod = 'vanRaden')


#' GWAResult is a data.table that has the following columns:
#
#' * trait:	   trait name
#' * snp:	     SNP name
#' * chr:	     chromosome on which the SNP is located
#' * pos:  	   position of the SNP on the chromosome
#' * allFreq:	 allele frequency of the SNP
#' * pValue: 	 P-value for the SNP
#' * effect: 	 effect of the SNP on the trait value
#' * effectSe: standard error of the effect of the SNP on the trait value
#' * RLR2:   	 likelihood-ratio-based R2 as defined in Sun et al. (2010)
#' * LOD:    	 LOD score (logarithm of the odds) for the SNP, defined as -log10(pValue)

View(GWAS$GWAResult$phenotypic_data)


#' Note that the estimated effect is computed for a single allele. Its direction 
#' depends on the coding of the markers in the gData object. In this example the 
#' minor allele was used as reference allele, so the effects are the estimated 
#' effects for the minor allele.
#
#' signSnp is a data.tables containing the significant SNPs. Optionally also the 
#' SNPs close to the significant SNPs are included in the data.table. 
#' The data.table in signSnp consist of the same columns as those in GWAResult 
#' described above. Two extra columns are added:
#
#' * snpStatus:	either significant SNP or within ... of a significant SNP
#' * propSnpVar:	proportion of the variance explained by the SNP

View(GWAS$signSnp$phenotypic_data)


################################################################################
#' get best marker name (highest LOD value)
snp <- GWAS$signSnp$phenotypic_data[LOD == max(LOD), snp]
#snp <- 'SCa5_29561931'


#' Manual check using t.test (old school way)
df <- as.data.frame(marker_matrix[,snp])

df$genotype  <- rownames(marker_matrix)
colnames(df) <- c('snp', 'genotype')

#' merge phenotyping data
df <- merge(df, phenotypic_data, id = 'genotype')

#' # remove all records with missing values
df <- na.omit(df) 

#' frequency table by snp or by score
table(df$snp)

#' remove all records with heterozygous SNPs
df <- df[df$snp %in% c('C', 'G', 'A', 'T'),]

t.test(ASC_Score ~ snp, data = df)
boxplot(ASC_Score ~ snp, data = df)


################################################################################
#' The kinship matrix used in the GWAS analysis.
#' The runSingleTraitGwas function has an argument kinshipMethod, which defines 
#' the kinship matrix used for association mapping. There are four options:
#
#' * astle (the default): see e.g. equation (2.2) in Astle and Balding (2009)
#' * IBS (Identity by State): see e.g. equation (2.3) in Astle and Balding (2009)
#' * vanRaden: using the formula by VanRaden (2008)
#' * User-defined, in which case the parameter kin needs to be specified.
#
#' By default, the same kinship matrix is used for testing all SNPs 
#' (GLSMethod = "single"). When GLSMethod = "multi", the kinship matrix is 
#' chromosome-specific. As shown by Rincent et al. (2014), this often gives a 
#' considerable improvement in power.

View(GWAS$kinship)


#' kinship matrix has nothing to do with the phenotyping data
# write.csv(GWAS$kinship, file = 'kinship_matrix.csv')

heatmap(GWAS$kinship)


#' performs a principal components analysis

pca <- prcomp(GWAS$kinship)


#' Cumulative Proportion

summary(pca)$importance[,1:2]


#' biplot shows correlations between traits

plot(pca$x[,1:2], pch = 20,
     main = paste0('Kinship PCA (', round(100*summary(pca)$importance[3,2],1), '%)'),
     xlab = paste0('PC1 (', round(100*summary(pca)$importance[2,1],1), '%)'),
     ylab = paste0('PC2 (', round(100*summary(pca)$importance[2,2],1), '%)'))


#' Additional information on the analysis, e.g. the call and the type of 
#' threshold used. GWASInfo$inflationFactor returns the inflation factor 
#' (Devlin and Roeder 1999). Ideally this factor should be 1, meaning there 
#' is no inflation at all. If the values are further away from 1, the 
#' inflation can be corrected for by setting genomicControl = TRUE in 
#' runSingleTraitGwas.

GWAS$GWASInfo


#' GWAS Summary
#' For a quick overview of the results, e.g. the number of significant SNPs, 
#' use the summary function.

summary(GWAS)


#' GWAS Plots, for more options: help(plot.GWAS)

#' A QQ-plot of the observed against the expected log10(p) values. Most of the 
#' SNPs are expected to have no effect, resulting in P-values uniformly 
#' distributed on (0,1), and leading to the identity function (y=x) 
#' on the log10(p) scale. Deviations from this line should only occur on the 
#' right side of the plot, for a small number of SNPs with an effect on the 
#' phenotype (and possibly SNPs in LD). 

plot(GWAS, plotType = 'qq', 'ASC_Score')


#' A Manhattan plot of GWAS. Significant SNPs are marked in red.
#' Plot only 5% of SNPs with a LOD below 2 (lod = 2)

plot(GWAS, plotType = 'manhattan', 'ASC_Score')

#' A qtl plot of GWAS. The significant SNPs are marked by circles at their 
#' genomic positions, with diameter proportional to the estimated effect size. 
#' Colors indicate the direction of the effect: green when the allele increases 
#' the trait value, and blue when it decreases the value.
plot(GWAS, plotType = 'qtl')


################################################################################
#' References:
#' 
#' Astle, William, and David J. Balding. 2009. Population Structure and Cryptic 
#' Relatedness in Genetic Association Studies. Statistical Science 24 (4): 
#' 451-71. https://doi.org/10.1214/09-sts307
#' 
#' Browning, Sharon R., and Brian L. Browning. 2007. Rapid and Accurate Haplotype 
#' Phasing and Missing-Data Inference for Whole-Genome Association Studies by Use 
#' of Localized Haplotype Clustering. The American Journal of Human Genetics 81 
#' (5): 1084-97. https://doi.org/10.1086/521987
#' 
#' Devlin, B., and Kathryn Roeder. 1999. Genomic Control for Association Studies. 
#' Biometrics 55 (4): 997-1004. https://doi.org/10.1111/j.0006-341x.1999.00997.x
#' 
#' Kang, Hyun Min, Noah A. Zaitlen, Claire M. Wade, Andrew Kirby, David Heckerman, 
#' Mark J. Daly, and Eleazar Eskin. 2008. Efficient Control of Population Structure 
#' in Model Organism Association Mapping. Genetics 178 (3): 1709-23. 
#' https://doi.org/10.1534/genetics.107.080101
#' 
#' Rincent, Renaud, Laurence Moreau, Herve Monod, Estelle Kuhn, Albrecht E. 
#' Melchinger, Rosa A. Malvar, Jesus Moreno-Gonzalez, et al. 2014. Recovering 
#' Power in Association Mapping Panels with Variable Levels of Linkage 
#' Disequilibrium. Genetics 197 (1): 375-87. 
#' https://doi.org/10.1534/genetics.113.159731
#' 
#' Sun, G, C Zhu, M H Kramer, S-S Yang, W Song, H-P Piepho, and J Yu. 2010. 
#' Variation Explained in Mixed-Model Association Mapping. Heredity 105 (4): 
#' 333-40. https://doi.org/10.1038/hdy.2010.11
#' 
#' Tunnicliffe, G Wilson. 1989. On the Use of Marginal Likelihood in Time Series 
#' Model Estimation. JRSS 51 (1): 15-27.
#' 
#' VanRaden, P. M. 2008. Efficient Methods to Compute Genomic Predictions. 
#' Journal of Dairy Science 91 (11): 4414-23. https://doi.org/10.3168/jds.2007-0980
