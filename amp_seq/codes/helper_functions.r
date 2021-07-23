


#######################################################################################################################
message("loading following functions:
        format_ps(), 
        copy_number_correction(),
        makeAlphaDiversity(),
        plot_beta_time(), 
        rarefy_even_sampling_depth()")
## Function for formatting the pseq object from dada2
format_ps <- function(ps){
  message("using the 'format_ps' function from helper_functions.r")
  taxic <- as.data.frame(ps@tax_table) 
  taxic$OTU <- rownames(taxic)  # Add the OTU ids from OTU table into the taxa table at the end.
  colnames(taxic) <- c("Domain",
                       "Phylum",
                       "Class", 
                       "Order", 
                       "Family",
                       "Genus", 
                       "Species",
                       "OTU")
  colnames(taxic)
  taxmat <- as.matrix(taxic)  # convert it into a matrix.
  new.tax <- tax_table(taxmat)  # convert into phyloseq compatible file.
  tax_table(ps) <- new.tax  # incroporate into phyloseq Object
  colnames(tax_table(ps))
  ps.sp <- aggregate_taxa(ps, "Species")
  head(tax_table(ps.sp)[1:5])
  
  message("manual check for sanity the file is stored in output/tables")
  write_phyloseq(ps.sp, 
                 "OTU", 
                 "output")
  
  table(tax_table(ps.sp)[, "Species"], 
        exclude = NULL)
  
  message("Remove unassigned ASVs") 
  ps <- subset_taxa(ps.sp, 
                    Species != "Unknown")
  return(ps)
}

# copy_number_correction

# @parm ps phyloseq object
# @parm column_with_ids Genus_species i.e column with names to match
# @parm copy_num_tab table with copynumbers 
copy_number_correction <- function(ps, column_with_ids, copy_num_tab) {
  counts <- copy_num_file <- column_with_ids <- corrected_tab <- tax.names <- NULL
  # Extract strain counts from phyloseq object
  counts <- as.data.frame(abundances(ps))
  copy_num_file <- copy_num
  column_with_ids = "Genus_species"
  intersect(rownames(counts), copy_num_file[, column_with_ids])
  tax.names <- rownames(counts)
  copy_num_file <- copy_num_file[order(match(copy_num_file$Genus_species, rownames(counts))), ]
  # The copy_number file needs to have a column matching the rownames of the otu_table.
  if (rownames(counts) != copy_num_file[, column_with_ids]) {
    print(
      "The copy_number file needs to have a column with matching values to the rownames of the otu_table."
    )
  }
  rownames(copy_num_file) <- copy_num_file[, column_with_ids]
  
  copy_num_file <- copy_num_file[tax.names, ]
  
  corrected_tab <- ceiling((counts) / copy_num_file[, "copy_number"])
  message("Check uncorrected table in output/corrected_tab.csv")
  write.csv(otu_table(ps), "output/not_corrected_tab.csv")  
  
  otu_table(ps) <- otu_table(as.matrix(corrected_tab), taxa_are_rows = T)
  message("Check corrected table in output/corrected_tab.csv")
  write.csv(corrected_tab, "output/corrected_tab.csv")
  #DT::datatable(otu_table(ps1))
  return(ps)
}



#######################################################################################################################

# rarefaction to even sampling depth #
# author: Doris Vandeputte           #
######################################
# this script doesn't include copy number correction, a function for copy number correction is included in RDP classifier 2.12 
# this script uses function rarefy_even_depth from phyloseq 1.20.0, it needs package phyloseq to be installed and loaded in order to work.
# with cnv_corrected_abundance_table: a copy number variation corrected abundance table with sample-identifiers as rows, copy number corrected taxa-abundances as columns
# with cell_counts_table: a table with sample-identifiers as rows, cell counts as columns 
library(phyloseq)
rarefy_even_sampling_depth <- function(cnv_corrected_abundance_table, cell_counts_table) 
{
  try(if(all(row.names(cnv_corrected_abundance_table) == row.names(cell_counts_table))==FALSE) stop("Cnv_corrected_abundance_table and cell_counts_table do not have the same sample-names, Please check!"))
  cnv_corrected_abundance_table = ceiling(cnv_corrected_abundance_table) # data values are rounded up in order to make use of integer values during the calculations
  cell_counts_table = t(cell_counts_table[order(row.names(cnv_corrected_abundance_table)),]) # make sure the order of the samples is the same in both files  
  sample_sizes = rowSums(cnv_corrected_abundance_table) # sample size of each sample (total nr of reads)
  sampling_depths = sample_sizes / cell_counts_table # sampling depth of each sample (total nr of reads divided by the cell count)
  minimum_sampling_depth = min(sampling_depths) # minimum of all sampling depths
  rarefy_to = cell_counts_table * minimum_sampling_depth # nr of reads to rarefy in each sample in order to get to an even sampling depth over all samples
  cnv_corrected_abundance_table_phyloseq = otu_table(cnv_corrected_abundance_table, taxa_are_rows = FALSE) # convert to phyloseq otutable
  rarefied_matrix=matrix(nrow = nrow(cnv_corrected_abundance_table_phyloseq), ncol = ncol(cnv_corrected_abundance_table_phyloseq), dimnames = list(rownames(cnv_corrected_abundance_table_phyloseq), colnames(cnv_corrected_abundance_table_phyloseq)))
  for (i in 1:nrow(cnv_corrected_abundance_table_phyloseq))
  {
    x <- rarefy_even_depth(cnv_corrected_abundance_table_phyloseq[i,], sample.size = rarefy_to[i], rngseed = 711, replace = FALSE, trimOTUs = F, verbose = FALSE)
    rarefied_matrix[i,] = x
  }
  normalised_rarefied_matrix = rarefied_matrix/rowSums(rarefied_matrix)
  QMP = normalised_rarefied_matrix*cell_counts_table[1,]
  return(QMP)
}

# Example
a = matrix( c(4,4,2,1,8,5,2,0,3,5,3,1,10,8,3,0,0,6,4,3), nrow=5, ncol=4, byrow = TRUE, dimnames = list(c("Sample A", "Sample B", "Sample C", "Sample D", "Sample E"),c("taxa1", "taxa2", "taxa3", "taxa4"))) # my cnv_corrected_abundance_table
b = matrix(c(10,20,34,21,12), nrow=5, ncol=1, byrow = TRUE, dimnames = list(c("Sample A", "Sample B", "Sample C", "Sample D", "Sample E"),c("#")))*100000 # my cell_counts_table
rarefy_even_sampling_depth(a,b)












##########################################################################################################
makeAlphaDiversity = function(abund, coverage) {
  relAb = abund / sum(abund)
  relAb = sort(relAb, decreasing = TRUE)
  cumprop = cumsum(relAb)
  min(which(cumprop >= coverage))
}

plot_diversity_coverage = function(abund, samples, coverage, prepost_cols) {
  diversity = apply(abund, 1, makeAlphaDiversity, coverage)
  ggplot(data.frame(diversity, sample_data(ps1.sub)))+ 
    geom_point(aes(x = Time_hr_num, color = Fasting, y =diversity, shape= Acetate_Feed), size=4, alpha=0.7) +
    facet_wrap(~ StudyIdentifier, nrow = 3) + #geom_smooth(aes(x = Time_hr_num, color = Fasting, y =diversity))  
    scale_color_manual(values = fasting_cols) +
    theme(panel.border = element_rect("transparent", size = 0.15)) +
    labs("x" = "Time (hr)", "col" = "Fasting") + ylab("No. of species") + scale_y_continuous(breaks=seq(0,16,1), limits = c(11,16))
  #+      
}

#######################################################################################################################
plot_beta_time <- function(ps=NA, title = NULL){
  #timepoint = NA, sample = "sample", t.start= NA, 
  require(dplyr)
  ps.sub <- df.x <- df <- s0 <- a <- st <- beta <- p1 <- NULL
  #sample_data(ps)$time <- sample_data(ps)[,timepoint]
  #colnames(sample_data(ps)) <- gsub("Time_hr_num", "time",colnames(sample_data(ps)) )
  #ps <- ps1.sub.b5
  #sample_data(ps)$sample <- sample_data(ps)[,sample]
  
  ps.sub <- subset_samples(ps, time >= t.start)
  
  df.x <-  meta(ps.sub)
  
  df <- df.x %>% arrange(time)
  
  s0 <- subset(df, time == t.start)$sample
  #s0 <- s0[,sample ]
  
  #sample_data(s0)$sample
  # Let us transform to relative abundance for Bray-Curtis calculations
  a <- abundances(microbiome::transform(ps.sub, "compositional")) 
  # Calculate the beta diversity between each time point and
  # the baseline (first) time point
  beta <- c() # Baseline similarity
  for (tp in df$time[-1]) {
    # Pick the samples for this subject
    # If the same time point has more than one sample,
    # pick one at random
    st <- sample(subset(df, time == tp)$sample, 1)
    #st <- st[,sample ]
    # Beta diversity between the current time point and baseline
    b <- vegdist(rbind(a[, s0], a[, st]), method = "bray")
    # Add to the list
    beta <- rbind(beta, c(tp, b))
  }
  colnames(beta) <- c("time", "beta")
  beta <- as.data.frame(beta)
  
  #p1 <- ggplot(beta, aes(x = time, y = beta)) + 
  #  geom_point() + geom_line()  + geom_smooth() #+ theme_bw()
  #p1 <- p1 + labs(x = 'Time (hr)', title = title,y = paste0("Community dissimilarity \ncompared to ",  t.start, "hr"), tag = 'a  ') 
  
  #+ ggtitle(title) + labs(y= ) #+ geom_label_repel()
  
  #return(p1)   
  return(beta)   
  
  
}
#######################################################################################################################

#https://github.com/ShadeLab/microbiome_trait_succession/blob/master/succession_custom_functions.R
#https://github.com/ShadeLab/microbiome_trait_succession/blob/master/succession_analysis.R


meltdist <- function(x, y, tree = NULL, method = c('bray','euclidean','unifrac', 'canberra')) {
  #meltdist takes a list of samples (x) and an otu table (y) and calculates intersample dissimilarity using (meth)
  
  #prune the original otu table
  otumat <- y[y$Sample %in% x, !names(y) %in% c('Sample','subject','Time_hr_num')] %>% as.data.frame()
  row.names(otumat) <- y$Sample[y$Sample %in% x]
  
  if (method == 'unifrac') {
    
    #create phyloseq object
    j <- phyloseq(otu_table(otumat, taxa_are_rows = FALSE),
                  drop.tip(tree, tree$tip.label[!tree$tip.label %in% names(otumat)]))
    
    #calculate UniFrac
    j <- UniFrac(j, weighted = TRUE)
    
  } else {
    
    #calculate bray-curtis or euclidean dissimilarity
    j <- vegdist(otumat, method = method)
    
  }
  
  #melt into a dataframe, append times
  j <- data.frame(melt(as.matrix(j)), stringsAsFactors = FALSE) %>%
    transmute(
      sample1 = Var1,
      sample2 = Var2,
      t1 = y$Time_hr_num[match(sample1, y$Sample)],
      t2 = y$Time_hr_num[match(sample2, y$Sample)],
      dist = value) %>%
    mutate_if(is.factor, as.character)
  
  return(j)
  
}

#######################################################################################################################

## For Taylor's law
prep_for_analysis <- function(pobj, normalize = c(TRUE,FALSE)){
  
  otu.mat <- otu.mat2 <- otu.mat3 <- NULL
  
  otu.mat <- as.data.frame(abundances(pobj)) 
  
  require(dplyr)
  colnames(otu.mat) <- sample_data(pobj)$Time_hr
  #head(otu.mat)
  otu.mat2 <- as.data.frame(t(otu.mat))
  #head(otu.mat2)
  otu.mat2$Time <- rownames(otu.mat2)
  otu.mat2$Time <- gsub("X", "", otu.mat2$Time)
  otu.mat2$Time <- as.numeric(otu.mat2$Time)
  otu.mat2 <- otu.mat2 %>% arrange(Time)
  rownames(otu.mat2) <- otu.mat2$Time
  otu.mat3 <- otu.mat2[,-17]
  
  if (normalize == TRUE){
    otu.mat3 <- normalize(t(otu.mat3), removeZero = TRUE)# compositional for Whisker neutrailty test
  } else {
    otu.mat3 <- otu.mat3
  }
  
  return(otu.mat3)
}

################################################################################
library("plotrix")
library("vegan")

titleSize <- 0.95

getMembers <- function(mergeTable, idx) {
  if (idx < 0) {
    return(-idx)
  } else {
    
    aList <- getMembers(mergeTable, mergeTable[idx, 1])
    bList <- getMembers(mergeTable, mergeTable[idx, 2])
    
    return(c(aList, bList))
  }
}

getReferencePhaseEnd <- function(data, begin, end) {
  myData <- data[data[1] >= begin & data[1] <= end,]
  rownames(myData) <- myData[[1]]
  myData <- myData[-1]
  myD <- dist(myData, method="canberra")
  hc <- hclust(myD)
  
  plot(hc, cex=0.6, main="Reference state detection", xlab=NA, sub=NA, font.main=1, cex.main=titleSize)
  
  startTime <- min(as.numeric(hc$labels))
  
  index <- match(startTime, hc$labels)
  
  nMerger <- nrow(hc$merge)
  startIdx <- 0
  
  for (i in 1:nMerger) {
    if (hc$merge[i,1] == -index) {
      startIdx <- hc$merge[i,2]
      break
    } else {
      if (hc$merge[i,2] == -index) {
        startIdx <- hc$merge[i,1]
        break
      }
    }
  }
  
  if (startIdx == 0) {
    print("Error!")
  }  
  
  myList <- getMembers(hc$merge, startIdx)
  
  value <- max(as.numeric(hc$labels[myList]))
  
  return(value)
}

#37482 F5T0
#299851
#counts <- as.data.frame(abundances(ps))
#copy_num_file <- copy_num
#column_with_ids = "Genus_species"
#intersect(rownames(counts), copy_num_file[, column_with_ids])
#tax.names <- rownames(counts)
#copy_num_file <- copy_num_file[order(match(copy_num_file$Genus_species, rownames(counts))), ]
# The copy_number file needs to have a column matching the rownames of the otu_table.
#if (rownames(counts) != copy_num_file[, column_with_ids]) {
#  print(
#    "The copy_number file needs to have a column with matching values to the rownames of the otu_table."
#  )
#}
#rownames(copy_num_file) <- copy_num_file[, column_with_ids]

#copy_num_file <- copy_num_file[tax.names, ]

#corrected_tab <- ceiling((counts) / copy_num_file[, "copy_number"])
#message("Check uncorrected table in output/corrected_tab.csv")
#write.csv(otu_table(ps), "output/not_corrected_tab.csv")  

#otu_table(ps) <- otu_table(as.matrix(corrected_tab), taxa_are_rows = T)
#message("Check corrected table in output/corrected_tab.csv")
#write.csv(corrected_tab, "output/corrected_tab.csv")
#DT::datatable(otu_table(ps1))

#dir.create("output/rds")
# add copy number corrected