# libraries ----
# Packages in Cran
cran_packages <- c(
  "tidyverse", "readxl", "writexl", "readr", "qqman", "vcfR", "QBMS", "adegenet",
  "ade4", "ggiraph", "ggpubr", "plotly", "poppr", "reactable",
  "rnaturalearth", "scatterpie", "snpReady", "viridis", "tibble",
  "ggplot2", "reshape2", "forcats", "dplyr", "sp", "scales", "htmltools", 
  "ASRgenomics", "statgenGWAS", "gplots", "spdep", "adespatial", "DT", "rrBLUP"
)

# Bioconductor Packages
bioc_packages <- c("rrBLUP", "LEA")

# Installing Cran Packages
installed <- rownames(installed.packages())
missing <- setdiff(cran_packages, installed)
if (length(missing)) {
  message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

# Installing Bioconductor Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
installed <- rownames(installed.packages())
missing <- setdiff(bioc_packages, installed)
if (length(missing)) {
  message("Installing missing Bioconductor packages: ", paste(missing, collapse = ", "))
  BiocManager::install(missing)
}

lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)

# aesthetics ----
clist <- list(
  "shiny"=c("#1D72F5","#DF0101","#77CE61", "#FF9326","#A945FF","#0089B2","#FDF060","#FFA6B2","#BFF217","#60D5FD","#CC1577","#F2B950","#7FB21D","#EC496F","#326397","#B26314","#027368","#A4A4A4","#610B5E"),
  "strong"=c("#11A4C8","#63C2C5","#1D4F9F","#0C516D","#2A2771","#396D35","#80C342","#725DA8","#B62025","#ED2224","#ED1943","#ED3995","#7E277C","#F7EC16","#F8941E","#8C2A1C","#808080"), 
  "funky" = c("#A6CEE3", "#3B8ABE", "#72B29C", "#84C868", "#4F9F3B","#EC9A91", "#E93E3F", "#F06C45", "#FDAC4F",  "#FB820F", "#D1AAB7", "#8C66AF", "#A99099", "#EEDB80", "#B15928"))

options(reactable.theme = reactableTheme(
  style = list(fontFamily = "Arial", fontSize = 12 ),
  color = "black"
))

# functions Java----
csvDownloadButton <- function(tableId, label = "Download as CSV", filename = "data.csv") {
  htmltools::tags$button(
    label,
    onclick = sprintf("Reactable.downloadDataCSV('%s', '%s')", tableId, filename)
  )
}

create_dt <- function(data, caption = NULL) {
  DT::datatable(data,
                caption = caption,
                options = list(pageLength = 10,
                               autoWidth = TRUE))
}

# functions R----
# from vcf
vcfToNumericMatrix <- function(vcfMatrix){
  X <- vcfMatrix
  
  matnumSNP <- ifelse(X == "0/0", as.numeric(0), X) # homozygote AA
  matnumSNP <- ifelse(matnumSNP == "0/1", as.numeric(1) ,matnumSNP) # heterozygote
  matnumSNP <- ifelse(matnumSNP == "1/0", as.numeric(1) ,matnumSNP) # heterozygote
  matnumSNP <- ifelse(matnumSNP == "1/1", as.numeric(2) ,matnumSNP) # homozygote BB
  
  matnumSNP <- matrix(as.numeric(matnumSNP), ncol = ncol(matnumSNP), dimnames = list(rownames(matnumSNP),colnames(matnumSNP)))
  return(matnumSNP)
}

## Quality control
# matrix is nxm matrix with n markers and m individuals in numeric format
# returns nxm matrix with n individuals and m markers
filterData <- function(matrix, call_rate = NULL, maf = NULL, na_ind = NULL, stats = FALSE){
  X <- matrix
  
  if(is.null(maf))
    maf = 0
  if(is.null(call_rate))
    call_rate = 0
  if(is.null(na_ind))
    na_ind = 0
  
  # starting dimensions
  n_start <- nrow(X)
  m_start <- ncol(X)
  
  X <- X[rowMeans(!is.na(X)) > call_rate,] #call rate
  stats_call_rate <- n_start - nrow(X)
  n_after_callrate <- nrow(X)
  
  X <- X[, colMeans(!is.na(X)) > na_ind] #ind missing data
  stats_na_ind <- m_start - ncol(X)
  
  mafFreq <- apply(X, 1, function(row) {
    row <- row[!is.na(row)]
    if(length(row) == 0) return(NA)
    maf <- sum(row) / (2 * length(row))
    min(maf, 1 - maf)
  })
  stats_maf <- nrow(X) - length(which(!is.na(mafFreq) & mafFreq > 0.01))
  X <- X[which(!is.na(mafFreq) & mafFreq > 0.01), ]
  
  if(stats){
    stats_df <- data.frame(
      parameter = c("call rate", "na ind", "maf"),
      value = c(stats_call_rate, stats_na_ind, stats_maf)
    )
    return(list(matrix = X, stats = stats_df))
  } else {
    return(X)
  }
}

## Marker info
# marker dataframe
# matrix is nxm matrix with n individuals and m markers
# returns a SNP density plot and a MAF value histogram
markerPlots <- function(markers, matrix, chrom = NULL){
  SNPinfo <- markers
  colnames(SNPinfo) <- c("id", "alleles", "chrom", "position")
  SNPinfo$position <- as.numeric(SNPinfo$position)
  SNPinfo <- SNPinfo[SNPinfo$id %in% colnames(matrix),]
  
  if(!is.null(chrom)){
    SNPinfo <- SNPinfo[SNPinfo$chrom == chrom,]
  }
  
  # SNP plot
  SNPPlot <- ggplot(data = SNPinfo) + 
    geom_density(aes(x = position), bw = 150000, linewidth = 0.4,  color = "dodgerblue4") +
    labs(title = "Location of SNPs", x = "SNP Position", y = "Density", color = "Chromosome") +
    scale_x_continuous(expand = c(0, 0), breaks = waiver(), n.breaks = 10) + 
    scale_y_continuous(expand = c(0, 0)) +
    geom_point(aes(x = position, y = rep(0, nrow(SNPinfo)), color = factor(chrom)), size = 1, shape = 3)  +
    scale_color_viridis(discrete = TRUE, option = "D", direction = 1) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.key.size = unit(0.2, 'cm'), legend.title = element_text(size=8), legend.text = element_text(size=6)) +
    guides(color = guide_legend(ncol = 8, bycol = TRUE))
  
  # MAF
  alleleFreq <- apply(matrix, 2, table)
  mafFreq <- lapply(alleleFreq, function(x) x/nrow(matrix)) # absolute frequencies
  mafValue <- lapply(mafFreq, function(x) min(x))
  mafDF <- data.frame(matrix(unlist(mafValue), nrow = length(mafValue), byrow = TRUE), stringsAsFactors = FALSE)
  names(mafDF) <- "maf"
  
  # MAF plot
  mafPlot <- ggplot(mafDF, aes(x = maf)) +
    geom_histogram(aes(y = after_stat(count)), bins = 10, boundary = 0, color = "black", fill = "#A6CEE3") +
    scale_x_continuous(breaks = waiver(), n.breaks = 10) +
    scale_y_continuous(breaks = waiver(), n.breaks = 7) +
    theme_classic() +
    labs(x = "Minor Allelic Frequency (MAF)", y = "Number of SNPs")
  
  # missingData
  missingData <- data.frame("missing" = as.vector(colMeans(is.na(matrix))))
  
  # missing data plot
  missingDataPlot <- ggplot(missingData, aes(x = missing)) +
    geom_histogram(aes(y = after_stat(count)), bins = 10, boundary = 0, color = "black", fill = "#A6CEE3") +
    scale_x_continuous(breaks = waiver(), n.breaks = 10) +
    scale_y_continuous(breaks = waiver(), n.breaks = 7) +
    theme_classic() +
    labs(x = "Missing data", y = "Number of SNPs") +
    theme(axis.text.x = element_text(size = 7))
  
  # Plots
  return(list(markerPlot = ggplotly(SNPPlot, tooltip = c("x","color")), mafPlot = ggplotly(mafPlot), missingData = ggplotly(missingDataPlot)))
}

## SNP ready
# geno is nxm numeric matrix with n individuals and m markers
# returns diversity parameters calculated with SNP ready package
genDivSNPReady <- function(geno, plots = FALSE){
  SNPDivStats <- popgen(geno, plot=FALSE)
  
  # diversity parameters by markers
  DivMarkers <- SNPDivStats[["whole"]][["Markers"]]
  DivMarkers <- rownames_to_column(DivMarkers,"Marker")
  
  DivMarkersTable <- htmltools::browsable(
    tagList(
      csvDownloadButton("DivMarkers", "Download as CSV", filename = "DivMarkers.csv"),
      
      reactable(DivMarkers, filterable = TRUE, searchable = TRUE, elementId = "DivMarkers")
    )
  )
  
  # diversity parameters by accessions
  AccessionsHo <- as.data.frame(SNPDivStats[["whole"]][["Genotypes"]][,1])
  AccessionsHo <- rownames_to_column(AccessionsHo, "Accession")
  colnames(AccessionsHo)[2] <- "Observed heterozygosity"
  
  AccessionsHoTable <- htmltools::browsable(
    tagList(
      csvDownloadButton("DivAccessions", "Download as CSV", filename = "DivAccessions.csv"),
      
      reactable(AccessionsHo, filterable = TRUE, searchable = TRUE, elementId = "DivAccessions")
    )
  )
  
  if(plots == FALSE){
    return(list(markers = DivMarkersTable, accessions = AccessionsHoTable))
  } else if(plots == TRUE)
  {
    HePlot <- ggplot(DivMarkers, aes(He)) +
      geom_histogram(aes(y = after_stat(count)), bins = 10, boundary = 0, color = "black", fill = "#A6CEE3") +
      scale_x_continuous(breaks = waiver(), n.breaks = 10) +
      scale_y_continuous(breaks = waiver(), n.breaks = 7) +
      theme_classic() +
      labs(x = "Expected Heterozigosity (He) ", y = "Number of SNPs") +
      theme(axis.text.x = element_text(size = 7))
    
    HoPlot <-ggplot(DivMarkers, aes(Ho)) +
      geom_histogram(aes(y = after_stat(count)), bins = 10, boundary = 0, color = "black", fill = "#A6CEE3") +
      scale_x_continuous(breaks = waiver(), n.breaks = 10) +
      scale_y_continuous(breaks = waiver(), n.breaks = 7) +
      theme_classic() +
      labs(x = "Observed Heterozigosity (Ho) ", y = "Number of SNPs") +
      theme(axis.text.x = element_text(size = 7))
    
    PICPlot <- ggplot(DivMarkers, aes(x = PIC)) +
      geom_histogram(aes(y = after_stat(count)), bins = 10, boundary = 0, color = "black", fill = "#A6CEE3") +
      scale_x_continuous(breaks = waiver(), n.breaks = 10) +
      scale_y_continuous(breaks = waiver(), n.breaks = 7) +
      theme_classic() +
      labs(x = "Polymorphic Information Content (PIC)", y = "Number of SNPs") +
      theme(axis.text.x = element_text(size = 7))
    
    MAFPlot <- ggplot(DivMarkers, aes(x = MAF)) +
      geom_histogram(aes(y = after_stat(count)), bins = 10, boundary = 0, color = "black", fill = "#A6CEE3") +
      scale_x_continuous(breaks = waiver(), n.breaks = 10) +
      scale_y_continuous(breaks = waiver(), n.breaks = 7) +
      theme_classic() +
      labs(x = "Minor Allele Frequency (MAF)", y = "Number of SNPs") +
      theme(axis.text.x = element_text(size = 7))
    plots <- subplot(HePlot, HoPlot, PICPlot, MAFPlot, nrows = 2, titleX = TRUE, titleY = TRUE, margin = 0.1)
    return(list(markers = DivMarkersTable, accessions = AccessionsHoTable, plots = ggplotly(plots)))
  }
}

## He by subgroups
# geno is nxm matrix with n individuals and m markers in "1/1" format
# subgroups is factor vector of n length
# returns He by groups, including optional plot
HeBySubgroups <- function(geno, subgroups, plot = FALSE){
  genInd <- df2genind(geno, sep = "/", ploidy = 2, type = "codom", pop = subgroups)
  
  # expected heterozygosity within populations
  HePop <- Hs(genInd)
  HePop <- data.frame(HePop)
  names(HePop) <- "He"
  
  HePopTable <- htmltools::browsable(
      tagList(
        csvDownloadButton("HePop", "Download as CSV", filename = "HePop.csv"),
        
        reactable(HePop, elementId = "HePop")
      )
    )
  
  if(plot == FALSE){
    return(HePopTable)
  } else if(plot == TRUE){
    HePopPlot <- ggplot(HePop, aes(x = rownames(HePop), y = He)) +
      geom_col(fill = "dodgerblue4") + 
      labs(title="Expected heterozygosity by subgroup", x = "Subgroups") +
      scale_x_discrete(labels = wrap_format(3)) +
      theme_classic() +
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 6))
    return(list(df = HePopTable, plot = ggplotly(HePopPlot, tooltip = "y")))
  }
}

## PCoA of genetic distances between populations
# geno is nxm numeric matrix with n individuals and m markers
# subgroups is factor vector of n length
# returns a PCoA object from matrix and genetic distances
genDistPop <- function(geno, subgroups, method = 1, PCoA = FALSE){
  genInd <- df2genind(geno, sep = "/", ploidy = 2, type = "codom", pop = subgroups)
  genPop <- genind2genpop(genInd)
  
  genDist <- dist.genpop(genPop, method = method, diag = TRUE, upper = FALSE)
  
  if(PCoA == TRUE){
    pco <- dudi.pco(genDist, nf = 3, scannf = FALSE)
    
    Subgroups <- str_wrap(levels(popSet), 10)
    PCoAplot <- ggplot(data = pco$li, aes(x = A1, y = A2, colour = Subgroups)) + 
      geom_point() + 
      labs(title = "PCoA", x = paste("PCo1"), y = paste("PCo2"), colour = "Subgroups") +
      scale_color_viridis(discrete = TRUE, direction = -1, option = "D", aesthetics = "colour") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(hjust = 0.5), legend.key.size = unit(0.1, 'cm'), legend.position = "bottom", legend.text = element_text(size = 6))
    Plot <- ggplotly(PCoAplot, tooltip = c("colour"))
    
    return(list(genDist = genDist, PCoA = pco, PCoAPlot = Plot))
  }else{
    return(genDist)
  }
}

## AMOVA
# geno is nxm numeric matrix with n individuals and m markers
# subgroups is factor vector of n length
AMOVA <- function(geno, subgroups){
  genInd <- df2genind(geno, sep = "/", ploidy = 2, type = "codom", pop = subgroups)
  genInd.mean <- missingno(genInd, type = "mean")
  strata(genInd.mean) <- data.frame(pop = pop(genInd.mean)) #single stratum
  amovaResult <- poppr.amova(genInd.mean, within = FALSE, ~pop)
  return(amovaResult)
}

## PCA
# geno is nxm numeric matrix with n individuals and m markers
# subgroups is factor vector of n length
# returns a PCA object from matrix
PCAFromMatrix <- function(geno, subgroups = NULL){
  if(is.null(subgroups)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno))
  } else if(!is.null(subgroups)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno), pop = subgroups)
  }
  
  # replace NAs
  x.genInd <- scaleGen(genInd, NA.method = "mean", scale = FALSE)
  
  # pca
  pca <- dudi.pca(x.genInd, center = FALSE, scale = FALSE, scannf = FALSE, nf = round(ncol(geno)/3, 0))
  
  # variance
  PCAvariance <- data.frame(PC = 1:length(pca[["eig"]]), Variance = (pca[["eig"]]/sum(pca[["eig"]])) * 100, CumulativeVar = cumsum((pca[["eig"]]/sum(pca[["eig"]])) * 100)) 
  
  # plot
  if(is.null(subgroups)){
    PCAplot <- ggplot(data = pca$li, aes(x = Axis1, y = Axis2, text = rownames(pca$tab))) + 
      geom_point() + 
      labs(title = "PCA", x = paste("PC1 (", round(PCAvariance$Variance[1], 2),")"), y = paste("PC2 (", round(PCAvariance$Variance[2], 2),")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))
    plot <- ggplotly(PCAplot, tooltip = "text")
  } else if(!is.null(subgroups)){
    Subgroups <- str_wrap(subgroups, 10)
    PCAplot <- ggplot(data = pca$li, aes(x = Axis1, y = Axis2, text = rownames(pca$tab), colour = Subgroups)) + 
      geom_point() + 
      labs(title = "PCA", x = paste("PC1 (", round(PCAvariance$Variance[1], 2),")"), y = paste("PC2 (", round(PCAvariance$Variance[2], 2),")"), colour = "Subgroups") +
      scale_color_viridis(discrete = TRUE, direction = -1, option = "D", aesthetics = "colour") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.title = element_text(size = 8, hjust = 0.5), legend.key.size = unit(0.1, 'cm'), legend.position = "bottom", legend.text = element_text(size = 6))
    plot <- ggplotly(PCAplot, tooltip = c("text", "colour"))
  }
  return(list(pca = pca, var = PCAvariance, plot = plot))
}

## DAPC
# geno is nxm numeric matrix with n individuals and m markers
# pca is an S3: pca, dudi type object
# maxK is a numerical value (maximum number of K groups)
# stat is type of test statistic to calculate ("BIC", "AIC", "WCC")
# returns test statistic results for different K values
kMeansStats <- function(geno, pca = NULL, maxK, stat){
  if(is.null(pca)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno)) 
    # replace NAs
    x.genInd <- scaleGen(genInd, NA.method = "mean", scale = FALSE)
    # pca
    pca <- dudi.pca(x.genInd, center = FALSE, scale = FALSE, scannf = FALSE, nf = round(ncol(geno)/3, 0))
  } else if(!is.null(pca)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno))
    pca <- pca
  }
  
  kmeans <- find.clusters(genInd, n.pca = round(ncol(geno)/3, 0), stat = stat, choose.n.clust = FALSE, max.n.clust = maxK, dudi = pca)
  
  plot <- ggplot(data = as.data.frame(kmeans$Kstat), aes(x = 1:maxK, y = kmeans$Kstat)) + 
    geom_point(shape = 21, color = "blue") +
    theme_classic() +
    labs(title = paste("Detection based on ", stat), x = "Number of clusters (K)", y = stat) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list(kmeans = kmeans, statPlot = plot))
}

# geno is nxm numeric matrix with n individuals and m markers
# pca is an S3: pca, dudi type object
# maxK is a numerical value (maximum number of K groups)
# subgroups is factor vector of n length
# returns dapc and optional dapc and composition plots
DAPC <- function(geno, krange, subgroups = NULL, pca = NULL, dapcPlot = FALSE){
  if(is.null(subgroups)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno))
  } else if(!is.null(subgroups)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno), pop = subgroups)
  }
  
  if(is.null(pca)){
    # replace NAs
    x.genInd <- scaleGen(genInd, NA.method = "mean", scale = FALSE)
    # pca
    pca <- dudi.pca(x.genInd, center = FALSE, scale = FALSE, scannf = FALSE, nf = round(ncol(geno)/3, 0))
  } else if(!is.null(pca)){
    pca <- pca
  }
  
  #dapc
  DAPCgrouplist <- vector(mode = "list", length = length(krange))
  DAPC <- vector(mode = "list", length = length(krange))
  
  for(i in 1:length(DAPC)){
    set.seed(10)
    DAPCgrouplist[[i]] <- find.clusters(genInd, n.pca = round(ncol(geno)/3, 0), n.clust = krange[i], dudi = pca, n.iter = 100)
    DAPC[[i]] <- dapc(genInd, pop = DAPCgrouplist[[i]]$grp, n.pca = round(ncol(geno)/3, 0), n.da = 5, dudi = pca, parallel = TRUE, var.contrib = TRUE, var.loadings = TRUE)
  }
  
  #dapc variance
  DAPCvar <- vector(mode = "list", length = length(DAPC))
  for(i in 1:length(DAPCvar)){
    DAPCvar[[i]] <- data.frame(DA = 1:length(DAPC[[i]][["eig"]]), Variance = (DAPC[[i]][["eig"]] / sum(DAPC[[i]][["eig"]])) * 100, CumulativeVar = cumsum((DAPC[[i]][["eig"]] / sum(DAPC[[i]][["eig"]])) * 100)) 
  }
  
  #dapc plot
  if(dapcPlot == TRUE){
    dapcPlot <- ggplot(data.frame(DAPC[[length(DAPC)]]$ind.coord, Group = DAPC[[length(DAPC)]]$grp), aes(x = LD1, y = LD2, color = Group, fill = Group, text = rownames(DAPC[[1]][["ind.coord"]]))) +
      geom_point(size = 2, shape = 21) +
      theme_classic() +
      scale_color_manual(values=clist$shiny) +
      scale_fill_manual(values=c(paste(clist$shiny, "66", sep = "")))
    plot <- ggplotly(dapcPlot, tooltip = "text")
    return(list(dapc = DAPC, var = DAPCvar, dapcPlot = plot))
  } else{
    return(list(dapc = DAPC, var = DAPCvar))
  }
}

# DAPC is DAPC object
# geno is nxm numeric matrix with n individuals and m markers
# subgroups is factor vector of n length
DAPCCompoPlot <- function(DAPC, geno, krange, subgroups = NULL){
  if(!is.null(subgroups)){
    # DAPC data frame for plot
    DAPCtemp <- as.data.frame(DAPC[[1]]$posterior)
    DAPCtemp$K <- krange[1]
    DAPCtemp$Individual <- rownames(DAPCtemp)
    # merge with passport data
    DAPCtemp <- merge(DAPCtemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE) 
    DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
    DAPCDF <- DAPCtemp
    
    for(i in 2:length(DAPC)){
      DAPCtemp <- as.data.frame(DAPC[[i]]$posterior)
      DAPCtemp$K <- krange[i]
      DAPCtemp$Individual <- rownames(DAPCtemp)
      DAPCtemp <- merge(DAPCtemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE)
      DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
      DAPCDF <- rbind(DAPCDF, DAPCtemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(DAPCDF, aes(x = Individual, y = Probability, fill = Group)) +
      geom_bar(stat = "identity") + 
      facet_grid(K ~ `Subgroup`, scales = "free_x", space = "free", labeller = labeller(K = grp.labs, `Subgroup` = label_wrap_gen(5))) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "DAPC Composition Plot", y = "Membership probability") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 4, angle = 90), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 4), legend.title = element_text(size = 5))
  } else if(is.null(subgroups)){
    # DAPC data frame for plot
    DAPCtemp <- as.data.frame(DAPC[[1]]$posterior)
    DAPCtemp$K <- krange[1]
    DAPCtemp$Individual <- rownames(DAPCtemp)
    DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
    DAPCDF <- DAPCtemp
    
    for(i in 2:length(DAPC)){
      DAPCtemp <- as.data.frame(DAPC[[i]]$posterior)
      DAPCtemp$K <- krange[i]
      DAPCtemp$Individual <- rownames(DAPCtemp)
      DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
      DAPCDF <- rbind(DAPCDF, DAPCtemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(DAPCDF, aes(x = Individual, y = Probability, fill = Group)) +
      geom_bar(stat = "identity") + 
      facet_grid(rows = vars(K), scales = "free_x", space = "free", labeller = labeller(K = grp.labs)) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "DAPC Composition Plot", y = "Membership probability") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 4), legend.title = element_text(size = 5))
  }
  return(CompoPlot)
}

# DAPC is DAPC object
# geno is nxm numeric matrix with n individuals and m markers
# subgroups is factor vector of n length
DAPCCompoPlotInt <- function(DAPC, geno, krange, subgroups = NULL){
  if(!is.null(subgroups)){
    # DAPC data frame for plot
    DAPCtemp <- as.data.frame(DAPC[[1]]$posterior)
    DAPCtemp$K <- krange[1]
    DAPCtemp$Individual <- rownames(DAPCtemp)
    # merge with passport data
    DAPCtemp <- merge(DAPCtemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE) 
    DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
    DAPCDF <- DAPCtemp
    
    for(i in 2:length(DAPC)){
      DAPCtemp <- as.data.frame(DAPC[[i]]$posterior)
      DAPCtemp$K <- krange[i]
      DAPCtemp$Individual <- rownames(DAPCtemp)
      DAPCtemp <- merge(DAPCtemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE)
      DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
      DAPCDF <- rbind(DAPCDF, DAPCtemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(DAPCDF, aes(x = Individual, y = Probability, fill = Group, tooltip = Individual)) +
      geom_bar_interactive (stat = "identity", color="white",size=0.02) + 
      facet_grid_interactive(K ~ `Subgroup`, scales = "free_x", space = "free", labeller = labeller(K = grp.labs, `Subgroup` = label_wrap_gen(5))) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "DAPC Composition Plot", y = "Membership probability") +
      theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 4, angle = 90), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 6))
    
  } else if(is.null(subgroups)){
    # DAPC data frame for plot
    DAPCtemp <- as.data.frame(DAPC[[1]]$posterior)
    DAPCtemp$K <- krange[1]
    DAPCtemp$Individual <- rownames(DAPCtemp)
    DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
    DAPCDF <- DAPCtemp
    
    for(i in 2:length(DAPC)){
      DAPCtemp <- as.data.frame(DAPC[[i]]$posterior)
      DAPCtemp$K <- krange[i]
      DAPCtemp$Individual <- rownames(DAPCtemp)
      DAPCtemp <- melt(DAPCtemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
      DAPCDF <- rbind(DAPCDF, DAPCtemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(DAPCDF, aes(x = Individual, y = Probability, fill = Group, tooltip = Individual)) +
      geom_bar_interactive(stat = "identity", color="white",size=0.02) + 
      facet_grid(rows = vars(K), scales = "free_x", space = "free", labeller = labeller(K = grp.labs, `Subgroup` = label_wrap_gen(5))) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "DAPC Composition Plot", y = "Membership probability") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 6))
  }
  plot <- girafe(ggobj = CompoPlot)
  return(plot)
}

# DAPC is DAPC object
# subgroups is factor vector of n length
# k is number of clusters to plot
DAPCFreqPlot <- function(DAPC, subgroups){
  dapc <- DAPC
  k <- dapc[["n.da"]]
  par(mar=c(1,1,1,1))
  return(table.value(table(subgroups, dapc[["assign"]]), row.labels = rownames(table(subgroups, dapc[["assign"]])), col.labels = 1:k, clabel.row = 0.5, clabel.col = 0.7, clegend = 0.7))
}

## sPCA
# coordinate conversion functions from https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
# xy 2 column longitude latitude vector or df
# geno is nxm numeric matrix with n individuals and m markers
# pca is an S3: pca, dudi type object
# maxK is a numerical value (maximum number of K groups)
# subgroups is factor vector of n length
# returns sPCA and optional eigenvalue plot, global and local tests, and map plot of results

sPCA <- function(geno, subgroups = NULL, xy, eigenPlot = TRUE, tests = TRUE){
  #coordinate conversion functions
  find_UTM_zone <- function(longitude, latitude) {
    #Special zones for Svalbard and Norway
    if (latitude >= 72.0 && latitude < 84.0 ) 
      if (longitude >= 0.0  && longitude <  9.0) 
        return(31);
    if (longitude >= 9.0  && longitude < 21.0)
      return(33)
    if (longitude >= 21.0 && longitude < 33.0)
      return(35)
    if (longitude >= 33.0 && longitude < 42.0) 
      return(37)
    (floor((longitude + 180) / 6) %% 60) + 1
  }
  find_UTM_hemisphere <- function(latitude) {
    ifelse(latitude > 0, "north", "south")
  }
  # returns a DF containing the UTM values, the zone and the hemisphere
  longlat_to_UTM <- function(id, long, lat, units = 'm') {
    df <- data.frame(
      id = id, 
      x = long, 
      y = lat
    )
    sp::coordinates(df) <- c("x", "y")
    hemisphere <- find_UTM_hemisphere(lat)
    zone <- find_UTM_zone(long, lat)
    sp::proj4string(df) <- sp::CRS("EPSG:4326")
    CRSstring <- paste0(
      "+proj=utm +zone=", zone,
      " +ellps=WGS84",
      " +", hemisphere,
      " +units=", units)
    if (dplyr::n_distinct(CRSstring) > 1L) 
      stop("multiple zone/hemisphere detected")
    
    res <- sp::spTransform(df, sp::CRS(CRSstring[1L])) %>%
      as.data.frame()
    
    res
  }
  
  if(is.null(subgroups)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno))
  } else if(!is.null(subgroups)){
    genInd <- new("genind", geno, indNames = rownames(geno), locNames = colnames(geno), pop = subgroups)
  }
  
  #coordinates
  xy <- data.frame(id = rownames(geno), longitude = xy[,1], latitude = xy[,2])
  xycoord <- do.call(rbind, apply(xy, 1, function(z) longlat_to_UTM(z[1], as.numeric(z[2]), as.numeric(z[3]))))
  xyjit <- jitter(cbind(xycoord$coords.x1, xycoord$coords.x2))
  
  # replace NAs
  x.genInd <- scaleGen(genInd, NA.method = "mean", scale = FALSE)
  
  #spca
  SPCA <- spca(x.genInd, xyjit, ask = FALSE, scannf = FALSE, type = 1, nfposi = 5, nfnega = 5, plot.nb = FALSE, printres = FALSE)
  
  if(eigenPlot == TRUE){
    #spca eigenvalue plot
    EigenPlot <- ggplot(data = data.frame(x = 1:length(SPCA$eig), eig = SPCA$eig), aes(x = x, y = eig, fill = spectral(length(eig)))) + 
      geom_col() +
      labs(title = "Eigenvalues of sPCA", y = "Eigenvalues") +
      theme_classic() + 
      scale_x_continuous(expand = c(0, 1)) + 
      scale_y_continuous(expand = c(0, 1)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), line = element_blank()) +
      guides(fill= FALSE)
  }
  
  if(tests == TRUE){
    #global
    Gtest <- global.rtest(x.genInd, SPCA$lw, nperm=99)
    GtestPlot <- plot(Gtest)
    
    #local
    Ltest <- local.rtest(x.genInd, SPCA$lw, nperm=99)
    LtestPlot <- plot(Ltest)
  }
  
  if(eigenPlot == TRUE & tests == TRUE){
    return(list(spca = SPCA, eigenPlot = EigenPlot, globalTest = list(test = Gtest, plot = GtestPlot), localTest = list(tes = Ltest, plot = LtestPlot)))
  } else if(eigenPlot == TRUE & tests == FALSE){
    return(list(spca = SPCA, eigenPlot = EigenPlot))
  } else if(eigenPlot == FALSE & tests == TRUE){
    return(list(spca = SPCA, globalTest = list(test = Gtest, plot = GtestPlot), localTest = list(tes = Ltest, plot = LtestPlot)))
  } else if(eigenPlot == FALSE & tests == FALSE & mapPlot == FALSE){
    return(list(spca = SPCA))
  }
}

## spca Plots
# spca is an spca object
# xy 2 column longitude latitude vector or df
# axis is the number of axis the user wants to plot
sPCAMapPlot <- function(spca, geno, xy, axis = 1, pos = TRUE){
  xy <- data.frame(id = rownames(geno), longitude = xy[,1], latitude = xy[,2])
  if(pos == TRUE){
    MapPlot <- ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf() + 
      geom_point(data = xy , aes(x = longitude, y = latitude, text = id, size = spca$li[,axis]^2, fill = factor(sign(spca$li[,axis])), 
                                 color = factor(sign(spca$li[,axis]))), shape = 22) +
      scale_fill_manual(values = c("white", "grey", "black"),
                        breaks = c("-1", "0", "1")) +
      scale_color_manual(values = c("black", "grey", "white"),
                         breaks = c("-1", "0", "1")) +
      coord_sf(xlim = range(xy$longitude), ylim = range(xy$latitude), expand = TRUE) +
      theme_classic() +
      labs(title = "sPCA - Positive Principal Component", x = "Longitude", y = "Latitude") +
      theme(plot.title = element_text(hjust = 0.5, size = 10), axis.line = element_line(color = "black", linewidth = 0.4), 
            panel.grid.major = element_line(color = "grey", linewidth = 0.2), axis.title = element_text(size = 10), 
            legend.key.size = unit(0.4, 'cm'), legend.title = element_text(size = 8), legend.text = element_text(size = 6), 
            axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom") +
      guides(color = "none", fill = "none", size = "none")
  }else{
    MapPlot <- ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf() + 
      geom_point(data = xy , aes(x = longitude, y = latitude, text = id, size = spca$li[,length(spca$li)+1-axis]^2, fill = factor(sign(spca$li[,length(spca$li)+1-axis])), 
                                 color = factor(sign(spca$li[,length(spca$li)+1-axis]))), shape = 22) +
      scale_fill_manual(values = c("white", "grey", "black"),
                        breaks = c("-1", "0", "1")) +
      scale_color_manual(values = c("black", "grey", "white"),
                         breaks = c("-1", "0", "1")) +
      coord_sf(xlim = range(xy$longitude), ylim = range(xy$latitude), expand = TRUE) +
      theme_classic() +
      labs(title = "sPCA - Negative Principal Component", x = "Longitude", y = "Latitude") +
      theme(plot.title = element_text(hjust = 0.5, size = 10), axis.line = element_line(color = "black", linewidth = 0.4), 
            panel.grid.major = element_line(color = "grey", linewidth = 0.2), axis.title = element_text(size = 10), 
            legend.key.size = unit(0.4, 'cm'), legend.title = element_text(size = 8), legend.text = element_text(size = 6), 
            axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "bottom") +
      guides(color = "none", fill = "none", size = "none")
  }
  plot <- ggplotly(MapPlot, tooltip = "text")
  return(plot)
}

## sNMF (LEA)
# xy 2 column longitude latitude vector or df

# subgroups is factor vector of n length

# geno is nxm numeric matrix with n individuals and m markers
# write geno type file modified function
write.geno.mod <- function(geno, output.file) 
{
  test_character <- function(name, param, default)
  {
    if(missing(param)) {
      if(is.null(default)) {
        p = paste("'",name,"' argument is missing.", sep="");
        stop(p)
      } else 
        return(default);
    } else {
      if(!is.character(param)) {
        p = paste("'",name,"' argument has to be of type character.", 
                  sep="");
        stop(p);
      }
    }
    return(param);
  }
  
  if(missing(geno))
    stop("'Geno' argument is missing.")
  else if (!(is.matrix(geno) || is.data.frame(geno) || is.vector(geno)))
    stop("'Geno' argument has to be of type matrix, data.frame or vector.")
  else if (is.vector(geno))
    geno = matrix(geno,ncol=1,nrow=length(geno))
  else if (is.data.frame(geno))
    geno = as.matrix(geno)
  
  output.file = test_character("output.file", output.file, NULL)
  
  geno[which(is.na(geno))] = 9
  geno[which(is.nan(geno))] = 9
  
  write.table(t(geno), output.file, col.names=FALSE,row.names=FALSE,sep="");
  return(output.file);
}

# maxK is a numerical value (maximum number of K groups)
# geno is nxm numeric matrix with n individuals and m markers
# file is the geno type file address
sNMFFunction <- function(geno, file, maxK, subgroups = NULL, cePlot = TRUE){
  
  snmfObjectKvar <- snmf(file, K=1:maxK, ploidy = 2, alpha = 100, entropy = TRUE, project = "new")
  snmfCrossEntr <- data.frame(matrix(nrow = maxK, ncol = 2))
  
  for(i in 1:maxK){
    snmfCrossEntr[i,1] <- i
    snmfCrossEntr[i,2] <- snmfObjectKvar@runs[[i]]@crossEntropy
  }
  
  colnames(snmfCrossEntr) <- c("K", "Cross-entropy")
  
  krange <- 2:maxK
  sNMFmatrix <- vector(mode = "list", length = length(krange))
  for(i in 1:length(sNMFmatrix)){
    sNMFmatrix[[i]] <- Q(snmfObjectKvar, K = krange[i])
    colnames(sNMFmatrix[[i]]) <- c(1:krange[i])
  }
  
  if(cePlot == TRUE){
    CrossEntrPlot <- ggplot(data = snmfCrossEntr, aes(x = K, y = `Cross-entropy`)) + 
      geom_point(shape = 21, color = "blue") +
      scale_x_continuous(breaks=seq(0, maxK, 1)) +
      theme_classic() +
      labs(title = "Detection based on cross-entropy", x = "Number of clusters (K)", y = "Cross-entropy") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  if(cePlot == TRUE){
    return(list(snmf = snmfObjectKvar, qmatrix = sNMFmatrix, crossEntropy = snmfCrossEntr, crossEntropyPlot = CrossEntrPlot))
  }else {
    return(list(snmf = snmfObjectKvar, qmatrix = sNMFmatrix, crossEntropy = snmfCrossEntr))
  }
}

sNMFCompoPlot <- function(sNMFmatrix, geno, krange, subgroups = NULL){
  if(!is.null(subgroups)){
    # NMF data frame
    sNMFTemp <- as.data.frame(sNMFmatrix[[1]])
    sNMFTemp$K <- krange[1]
    sNMFTemp$Individual <- rownames(geno)
    sNMFTemp <- merge(sNMFTemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE)
    sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
    sNMFdf <- sNMFTemp
    
    for(i in 2:length(sNMFmatrix)){
      sNMFTemp <- as.data.frame(sNMFmatrix[[i]])
      sNMFTemp$K <- krange[i]
      sNMFTemp$Individual <- rownames(geno)
      sNMFTemp <- merge(sNMFTemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE)
      sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
      sNMFdf <- rbind(sNMFdf, sNMFTemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(sNMFdf, aes(x = Individual, y = Probability, fill = Group, tooltip = Individual)) +
      geom_bar(stat = "identity") + 
      facet_grid(K ~ `Subgroup`, scales = "free_x", space = "free", labeller = labeller(K = grp.labs, `Subgroup` = label_wrap_gen(5))) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "sNMF Composition Plot", y = "Membership probability") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 4, angle = 90), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 4), legend.title = element_text(size = 5))
  }
  if(is.null(subgroups)){
    # NMF data frame
    sNMFTemp <- as.data.frame(sNMFmatrix[[1]])
    sNMFTemp$K <- krange[1]
    sNMFTemp$Individual <- rownames(geno)
    sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
    sNMFdf <- sNMFTemp
    
    for(i in 2:length(sNMFmatrix)){
      sNMFTemp <- as.data.frame(sNMFmatrix[[i]])
      sNMFTemp$K <- krange[i]
      sNMFTemp$Individual <- rownames(geno)
      sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
      sNMFdf <- rbind(sNMFdf, sNMFTemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(sNMFdf, aes(x = Individual, y = Probability, fill = Group, tooltip = Individual)) +
      geom_bar(stat = "identity") + 
      facet_grid(vars(K), scales = "free_x", space = "free", labeller = labeller(K = grp.labs)) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "sNMF Composition Plot", y = "Membership probability") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 6, angle = 90), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 4), legend.title = element_text(size = 5))
  }
  return(CompoPlot)
}

sNMFCompoPlotInt <- function(sNMFmatrix, geno, krange, subgroups = NULL){
  if(!is.null(subgroups)){
    # NMF data frame
    sNMFTemp <- as.data.frame(sNMFmatrix[[1]])
    sNMFTemp$K <- krange[1]
    sNMFTemp$Individual <- rownames(geno)
    sNMFTemp <- merge(sNMFTemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE)
    sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
    sNMFdf <- sNMFTemp
    
    for(i in 2:length(sNMFmatrix)){
      sNMFTemp <- as.data.frame(sNMFmatrix[[i]])
      sNMFTemp$K <- krange[i]
      sNMFTemp$Individual <- rownames(geno)
      sNMFTemp <- merge(sNMFTemp, data.frame(Individual = rownames(geno), Subgroup = subgroups), by="Individual", all.x=TRUE)
      sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K", "Subgroup"))
      sNMFdf <- rbind(sNMFdf, sNMFTemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(sNMFdf, aes(x = Individual, y = Probability, fill = Group, tooltip = Individual)) +
      geom_bar_interactive(stat = "identity") + 
      facet_grid(K ~ `Subgroup`, scales = "free_x", space = "free", labeller = labeller(K = grp.labs, `Subgroup` = label_wrap_gen(5))) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "sNMF Composition Plot", y = "Membership probability") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 6, angle = 90), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 6))
  }
  if(is.null(subgroups)){
    # NMF data frame
    sNMFTemp <- as.data.frame(sNMFmatrix[[1]])
    sNMFTemp$K <- krange[1]
    sNMFTemp$Individual <- rownames(geno)
    sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
    sNMFdf <- sNMFTemp
    
    for(i in 2:length(sNMFmatrix)){
      sNMFTemp <- as.data.frame(sNMFmatrix[[i]])
      sNMFTemp$K <- krange[i]
      sNMFTemp$Individual <- rownames(geno)
      sNMFTemp <- melt(sNMFTemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
      sNMFdf <- rbind(sNMFdf, sNMFTemp)
    }
    
    grp.labs <- paste("K =", krange)
    names(grp.labs) <- krange
    
    CompoPlot <- ggplot(sNMFdf, aes(x = Individual, y = Probability, fill = Group, tooltip = Individual)) +
      geom_bar_interactive(stat = "identity") + 
      facet_grid(vars(K), scales = "free_x", space = "free", labeller = labeller(K = grp.labs, `Subgroup` = label_wrap_gen(5))) +
      scale_fill_manual(values=clist$shiny) +
      labs(title = "sNMF Composition Plot", y = "Membership probability") +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.x = element_text(size = 6, angle = 90), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 6))
  }
  plot <- girafe(ggobj = CompoPlot)
  return(plot)
}


# xlim optional 2 value vector establishing min and max longitude values for visualization
# ylim optional 2 value vector establishing min and max latitude values for visualization
sNMFMapPlot <- function(geno, sNMFObjectVar, xy, k, Xlim = NULL, Ylim = NULL){
  LEAplotDF <- data.frame(id = rownames(geno), longitude = xy[,1], latitude = xy[,2])
  sNMFmatrix <- Q(sNMFObjectVar, K = k)
  LEAplotDF<- cbind(LEAplotDF, sNMFmatrix)
  groups <- paste("G", 1:k, sep = "")
  names(LEAplotDF)<- c("id", "longitude", "latitude", groups)
  
  if(is.null(Xlim) & is.null(Ylim)){
    mapPlot <- ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf() +
      geom_scatterpie(data = LEAplotDF, aes(x = longitude, y = latitude), cols = groups, pie_scale = 0.5, size = 0.1) +
      coord_sf(xlim = range(LEAplotDF$longitude), ylim = range(LEAplotDF$latitude), expand = TRUE) +
      scale_fill_manual(values = clist$shiny) +
      theme_classic() +
      labs(title = "LEA map plot", x = "Longitude", y = "Latitude") +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.line = element_line(color = "black", linewidth = 0.4), panel.grid.major = element_line(color="grey", linewidth = 0.2), axis.title = element_text(size = 10), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 10), legend.text = element_text(size = 8))
  } else{
    A <- diff(range(LEAplotDF$longitude))*diff(range(LEAplotDF$latitude))
    
    mapPlot <- ggplot(data = ne_countries(scale = "medium", returnclass = "sf")) +
      geom_sf() +
      geom_scatterpie(data = LEAplotDF, aes(x = longitude, y = latitude), cols = groups, pie_scale = 5 * (diff(range(Xlim))*diff(range(Ylim)))/A, size = 0.1) +
      coord_sf(xlim = Xlim, ylim = Ylim, expand = TRUE) +
      scale_fill_manual(values = clist$shiny) +
      theme_classic() +
      labs(title = "LEA map plot", x = "Longitude", y = "Latitude") +
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.line = element_line(color = "black", linewidth = 0.4), panel.grid.major = element_line(color="grey", linewidth = 0.2), axis.title = element_text(size = 10), legend.key.size = unit(0.5, 'cm'), legend.title = element_text(size = 10), legend.text = element_text(size = 8))
  }
  
  return(mapPlot)
}

kOrderedPlot <- function(df){
  # sortInd() function from PopHelper
  
  #' @title Internal: Handles individual sorting
  #' @description Internal: Handles individual sorting
  #' @param dframe A q-matrix dataframe
  #' @param grplab A dataframe with group labels
  #' @param selgrp A single character denoting selected group label set.
  #' See details.
  #' @param ordergrp A logical indicating if individuals must be ordered by all
  #' group labels
  #' @param sortind A character indicating how individuals are sorted. Default is
  #' NA (Same order of individuals as in input file). Other options are 'all'
  #' (sorting by values of all clusters), by any one cluster (eg. 'Cluster1') or
  #' 'labels' (sorting by individual labels).
  #' @param grplabpos A numeric indicating y-axis position of labels
  #' @param linepos A numeric indicating y-axis position of label line
  #' @param corder Current order of individuals
  #' @return Returns a list with ordered q-matrix dataframe and ordered grplab.
  #' @details When multiple group label sets are in use, \code{selgrp} defines
  #' which group label set is used for group ordering (\code{ordergrp}),
  #' subsetting (\code{subsetgrp}) and group mean (\code{grpmean}). \code{selgrp}
  #' is also used for plotting divider lines and and sorting (\code{sortind}).
  #' If \code{selgrp} is not specified, the first group label set is used by
  #' default.
  #' @noRd
  #' @keywords internal
  #'
  sortInd <- function(dframe=NULL,grplab=NA,selgrp=NA,ordergrp=FALSE,sortind=NA,grplabpos=NA,linepos=NA,corder=NA)
  {
    if(is.null(dframe)) stop("sortInd: Argument 'dframe' is empty.")
    if(!all(is.na(grplab)))
    {
      if(is.na(selgrp)) selgrp <- names(grplab)[1]
      if(!is.character(selgrp)) stop("sortInd: Argument 'selgrp' must be a character datatype.")
      if(length(selgrp)>1) stop("sortInd: Argument 'selgrp' must be a character datatype of length 1.")
      if(!any(selgrp %in% names(grplab))) stop(paste0("sortInd: Argument 'selgrp' contains (",selgrp,") which is not in the 'grplab' titles (",paste0(names(grplab),collapse=", "),")."))
    }
    if(all(!is.na(sortind)))
    {
      if(length(sortind) > 1) stop("sortInd: Argument 'sortind' must be of length 1. Use 'all','label' or a cluster name like 'Cluster1'.")
      #if(sortind != "all" && sortind != "label" && !grepl("Cluster[0-9+]",sortind)) stop("sortInd: Argument 'sortind' must be set to 'all', 'label' or a cluster like 'Cluster1'.")
    }
    if(is.na(grplabpos)) grplabpos <- 0.25
    if(is.na(linepos)) linepos <- 0.75
    
    # get original order
    if(any(is.na(corder))) corder <- 1:nrow(dframe)
    if(length(corder)!=nrow(dframe)) stop("grpLabels: Length of 'corder' not equal to number of individuals.")
    # sorting without grplab
    if(any(is.na(grplab)))
    {
      if(!is.na(sortind))
      {
        if(sortind=="all")
        {
          dftemp <- dframe
          dftemp$maxval <- as.numeric(apply(dframe,1,max))
          dftemp$matchval <- as.numeric(apply(dframe,1,FUN=function(x) match(max(x),x)))
          dframe$corder <- corder
          dframe <- dframe[with(dftemp,order(matchval,-maxval)),,drop=FALSE]
          corder <- dframe$corder
          dframe$corder <- NULL
          rm(dftemp)
        }
        
        if(sortind=="label")
        {
          dframe$corder <- corder
          dframe <- dframe[order(rownames(dframe)),,drop=FALSE]
          corder <- dframe$corder
          dframe$corder <- NULL
        }
        
        if(sortind!="all" && sortind!="label")
        {
          if(!(sortind %in% colnames(dframe))) stop(paste0("sortInd: 'sortind' value (",sortind,") not found in file header (",paste0(colnames(dframe),collapse=", "),")."))
          dframe$corder <- corder
          dframe <- dframe[order(dframe[[sortind]]),,drop=FALSE]
          corder <- dframe$corder
          dframe$corder <- NULL
        }
        
      }
      label_position <- NA
      marker_position <- NA
    }
    
    # sorting with grplab
    if(!all(is.na(grplab)))
    {
      gnames <- names(grplab)
      onames <- setdiff(gnames,selgrp)
      
      if(!is.na(sortind))
      {
        # sort by all
        if(sortind=="all")
        {
          dftemp <- dframe
          # find the max value for each individual
          dftemp$maxval <- as.numeric(apply(dframe,1,max))
          # pick cluster with max value
          dftemp$matchval <- as.numeric(apply(dframe,1,FUN=function(x) match(max(x),x)))
          dftemp$corder <- corder
          
          if(length(intersect(colnames(dftemp),colnames(grplab)))!=0) stop(paste0("sortInd: One or more header labels in the run file are duplicated in grplab header. Change labels to be unique. Following are the duplicate label(s): (",paste0(intersect(colnames(dftemp),colnames(grplab)),collapse=", "),")."))
          dftemp <- cbind(dftemp,grplab)
          
          if(ordergrp)
          {
            sort_asc <- c(selgrp,onames,"matchval")
            sort_desc <- "maxval"
            dframe <- dftemp[do.call(order,c(as.list(dftemp[sort_asc]),lapply(dftemp[sort_desc],function(x) -xtfrm(x)))),]
          }else{
            rle1 <- rle(as.character(unlist(grplab[selgrp])))
            grplabnames <- rle1$values
            tovec <- cumsum(rle1$lengths)
            fromvec <- (tovec - rle1$lengths)+1
            dftemplist <- vector("list",length=length(grplabnames))
            for(k in 1:length(tovec))
            {
              dftemp1 <- dftemp[fromvec[k]:tovec[k],,drop=FALSE]
              dftemp1$grp <- NULL
              dftemplist[[k]] <- dftemp1[with(dftemp1,order(matchval,-maxval)),,drop=FALSE]
            }
            dframe <- do.call("rbind",dftemplist)
          }
          
          corder <- dframe$corder
          dframe$corder <- NULL
          grplab <- dframe[,gnames,drop=FALSE]
          dframe[,gnames] <- NULL
          dframe$maxval <- NULL
          dframe$matchval <- NULL
        }
        
        # sort by label
        if(sortind=="label")
        {
          if(length(intersect(colnames(dframe),colnames(grplab)))!=0) stop(paste0("sortInd: One or more header labels in the run file are duplicated in grplab header. Change labels to be unique. Following are the duplicate label(s): (",paste0(intersect(colnames(dframe),colnames(grplab)),collapse=", "),")."))
          dftemp <- cbind(dframe,grplab)
          dftemp$corder <- corder
          
          if(ordergrp)
          {
            dftemp$label <- rownames(dftemp)
            sort_asc <- c(selgrp,onames,"label")
            dframe <- dftemp[do.call(order,dftemp[,sort_asc]),]
            dframe$label <- NULL
          }else{
            rle1 <- rle(as.character(unlist(grplab[selgrp])))
            grplabnames <- rle1$values
            tovec <- cumsum(rle1$lengths)
            fromvec <- (tovec - rle1$lengths)+1
            dftemplist <- vector("list",length=length(grplabnames))
            for(k in 1:length(tovec))
            {
              dftemp1 <- dftemp[fromvec[k]:tovec[k],,drop=FALSE]
              dftemp1$grp <- NULL
              dftemplist[[k]] <- dftemp1[order(rownames(dftemp1)),,drop=FALSE]
            }
            dframe <- do.call("rbind",dftemplist)
          }
          
          corder <- dframe$corder
          dframe$corder <- NULL
          grplab <- dframe[,gnames,drop=FALSE]
          dframe[,gnames] <- NULL
        }
        
        # sort by cluster
        if(sortind!="all" && sortind!="label")
        {
          if(!(sortind %in% colnames(dframe))) stop(paste0("sortInd: 'sortind' value (",sortind,") not found in file header (",paste0(colnames(dframe),collapse=", "),")."))
          # checks if sortind variable is a column in dframe
          if(length(intersect(colnames(dframe),colnames(grplab)))!=0) stop(paste0("sortInd: One or more header labels in the run file are duplicated in grplab header. Change labels to be unique. Following are the duplicate label(s): (",paste0(intersect(colnames(dframe),colnames(grplab)),collapse=", "),")."))
          dftemp <- cbind(dframe,grplab)
          dftemp$corder <- corder
          
          if(ordergrp)
          {
            sort_asc <- c(selgrp,onames,sortind)
            dframe <- dftemp[do.call(order,dftemp[,sort_asc]),]
          }else{
            rle1 <- rle(as.character(unlist(grplab[selgrp])))
            grplabnames <- rle1$values
            tovec <- cumsum(rle1$lengths)
            fromvec <- (tovec - rle1$lengths)+1
            dftemplist <- vector("list",length=length(grplabnames))
            for(k in 1:length(tovec))
            {
              dftemp1 <- dftemp[fromvec[k]:tovec[k],,drop=FALSE]
              dftemp1$grp <- NULL
              dftemplist[[k]] <- dftemp1[order(dftemp1[[sortind]]),,drop=FALSE]
            }
            dframe <- do.call("rbind",dftemplist)
          }
          
          corder <- dframe$corder
          dframe$corder <- NULL
          grplab <- dframe[,gnames,drop=FALSE]
          dframe[,gnames] <- NULL
        }
      }
      
      # create label_position and marker_position
      gnames <- names(grplab)
      marker_position_list <- vector("list",length=length(gnames))
      label_position_list <- vector("list",length=length(gnames))
      intlablist <- vector("list",length=length(gnames))
      for(k in seq_along(gnames))
      {
        rlegrp <- rle(unlist(grplab[gnames[k]]))
        label_position_df <- data.frame(label=rlegrp$values,freq=rlegrp$lengths,stringsAsFactors=FALSE)
        marker_position_df <- data.frame(markerxpos=c(0,cumsum(label_position_df$freq)),stringsAsFactors=FALSE)
        label_position_df$labxpos <- round((diff(marker_position_df$markerxpos)/2)+marker_position_df$markerxpos[1:length(marker_position_df$markerxpos)-1],1)
        label_position_df$labypos <- rep(grplabpos,nrow(label_position_df))
        #marker_position_df$temp <- factor(rep(1,nrow(marker_position_df)))
        marker_position_df$markerypos <- rep(linepos,nrow(marker_position_df))
        
        marker_position_df$count <- gnames[k]
        marker_position_list[[k]] <- marker_position_df
        label_position_df$count <- gnames[k]
        label_position_list[[k]] <- label_position_df
      }
      
      marker_position <- do.call("rbind",marker_position_list)
      marker_position$count <- factor(marker_position$count,levels=gnames)
      label_position <- do.call("rbind",label_position_list)
      label_position$count <- factor(label_position$count,levels=gnames)
      
      rownames(marker_position) <- 1:nrow(marker_position)
      rownames(label_position) <- 1:nrow(label_position)
      
      #adjust divider position
      marker_position$divxpos <- marker_position$markerxpos+0.5
    }
    return(as.data.frame(dframe))
  }
  
  df <- df
  dfSorted <- sortInd(df, sortind="all")
  
  dftemp <- dfSorted
  dftemp$K <- ncol(dftemp)
  dftemp$Individual <- rownames(dfSorted)
  dftemp <- melt(dftemp, variable.name = "Group", value.name = "Probability", id = c("Individual", "K"))
  df1 <- dftemp
  
  OrderedPlot <- ggplot(df1, aes(x = fct_inorder(Individual), y = Probability, fill = Group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values=clist$funky) +
    labs(title = "Ordered Composition Plot", y = "Membership probability") +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), strip.text.y = element_text(size = 6), axis.text.y = element_text(size = 6), legend.text = element_text(size = 6))
  return(OrderedPlot)
}

kinshipMatrix <- function(geno, method = "vanRaden", save = FALSE){
  # mean value imputation
  SNPmeans <- apply(geno, 2, mean, na.rm=TRUE) 
  
  # imputing median value with NA  
  matrixMean <- geno
  for(i in 1:ncol(matrixMean)){
    matrixMean[,i][is.na(matrixMean[,i])] <- SNPmeans[i] 
  }
  
  kinshipMat <- kinship(matrixMean, method)
  
  if(save == TRUE){
    write.csv(kinshipMat, file = "kinship.csv")
  }
  
  return(kinshipMat)
}

# geno is a m x n allele matrix
# threshold is a similarity percentage
# method can be "astle", "IBS", "vanRaden", "identity"
kinshipDuplicates <- function(geno, threshold, method = "vanRaden", save = FALSE, kinship = NULL){
  if(!is.null(kinship)){
    kinshipMat <- kinship
  }else{
    # mean value imputation
    SNPmeans <- apply(geno, 2, mean, na.rm=TRUE) 
    
    # imputing median value with NA  
    matrixMean <- geno
    for(i in 1:ncol(matrixMean)){
      matrixMean[,i][is.na(matrixMean[,i])] <- SNPmeans[i] 
    }
    
    kinshipMat <- kinship(matrixMean, method)
  }
  
  kinshipInfo <- kinship.diagnostics(K = kinshipMat, duplicate.thr = threshold, plots = TRUE, message = FALSE)
  duplicateDF <- kinshipInfo$list.duplicate
  
  if(save == TRUE){
    write.csv(kinshipMat, file = "kinship.csv")
  }
  
  return(list(kinshipMatrix = kinshipMat, potentialDuplicates = duplicateDF, plots = list(kinshipInfo$plot.diag, kinshipInfo$plot.offdiag)))
}


kinshipHeatmap <- function(kinship, file){
  png(file, width = 3000, height = 3000)
  par(mar=c(1,1,1,1))
  heatmap.2(kinship, trace = "none", keysize = 0.5)
  dev.off()
}

kinshipFilter <- function(matrix, duplicatesdf, kinship, save = FALSE){
  kMatrix <- matrix[!rownames(matrix) %in% unique(duplicatesdf$Indiv.A),]
  cKinship <- kinship[!rownames(kinship) %in% unique(duplicatesdf$Indiv.A),!colnames(kinship) %in% unique(duplicatesdf$Indiv.A)]
  
  if(save == TRUE){
    write.csv(kMatrix, file = "cleanMatrix.csv")
    write.csv(cKinship, file = "cleanKinship.csv")
  }
  
  return(list(markerMatrix = kMatrix, kinship = cKinship))
}

phyloTree <- function(geno, treeType, distanceType, samples, path){
  glSNP <- new("genlight", geno, indNames = rownames(geno), locNames = colnames(geno), parallel=FALSE)
  alleleFreq <- tab(glSNP, freq = TRUE)
  tree <- aboot(alleleFreq, tree = treeType, distance = distanceType, sample = samples, showtree = FALSE)
  ape::write.tree(tree, file = path)
  return(tree)
}