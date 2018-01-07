###################################################################
### This script is designed to analyze Ribomap stats output files
###################################################################

# Step 1 - Set working directory and increase memory size for xlsx file writing
setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Pulsed_silac_Ozlem/Ribosome_profiling")



# Step 2 - Read in stats file 
file_names = grep("^GTML5", dir(), value = TRUE)
raw.data = lapply(file_names, #sapply preserves file names
                  function(x) read.csv(x, header = FALSE, colClasses = "character"))



# Step 3 - Data clean-up
## Restructure data into data frame format
require(tidyr) #data wrangling
require(dplyr) #data wrangling
raw.data = lapply(raw.data,
                  function(x) {
                    # number each block of transcript data
                    x$element = rep(1:(nrow(x)/5), each = 5) 
                    
                    # Rearrange data into columns
                    x = x %>%
                      group_by(element) %>%
                      summarise(data = paste(V1, collapse = ",")) %>%   #one transcript per row
                      separate(data, into = c("refID", "tid", "rabd", "tabd", "te"), 
                               sep = ",")
                    
                    # remove key terms in front of characters
                    x = as.data.frame(sapply(x, function(y) sub("^.{2,5}: ", "", y))) 
                    
                    # separate transcript headers into columns
                    x = separate(x, tid, into = c("transcript_id", "gene_id", "HAVANA_gene", 
                                                  "HAVANA_transcript", "symbol1", "symbol2",
                                                  "transcript_length", "architecture"),
                                 sep = "\\|", extra = "merge") 
                    
                    # remove terminal pipes in architecture column
                    x$architecture = sub("\\|$", "", x$architecture) 
                    
                    # convert factors to numeric
                    id = c(2, 11:13)
                    x[,id] = as.numeric(as.character(unlist(x[,id])))
                    return(x[-1]) #removes element column
                  })


## Assign sample-specific suffix to column names (DMSO, Rapa, M1071)
require(stringr) #extract string from pattern
file_labels = str_extract(file_names, "(?<=GTML5_).*-[123]")
raw.data = lapply(1:length(file_labels),
                  function(i) {
                    df = raw.data[[i]]
                    suffix = paste0("_", file_labels[i])
                    names(df) = paste0(names(df), c(rep("", 9), rep(suffix, 3)))
                    return(df)
                  })
names(raw.data) = file_labels


## Remove data with ribosome footprint < 128 
data = lapply(raw.data, function(x) x[x$rabd >= 128, ])



# Step 4 - Data-specific manipulations
## Organize data by treatment group
merge_treatment = function(pattern, data) {
  # pattern = name of treatment matching to names(data)
  # data = list of data frames containing pre-processed data
  name = grep(pattern, names(data), value = TRUE)
  data_list = data[name]
  df = Reduce(function(df1, df2) inner_join(df1, df2, 
                                            by = names(df1[1:9])),
              data_list)
  return(df)
}
DMSO = merge_treatment("DMSO", data)
M1071 = merge_treatment("M1071", data)
Rapa = merge_treatment("Rapamycin", data)
all = merge_treatment(".*", data)


## Calculate log2 and replicate averages of ribomap statistics
log2_compute = function(df, pattern = ".*(?=-[[:digit:]])") {
  # df = data frame containing rabd, tabd, te data
  # pattern = regex expression describing how to extract the replicates
  label = grep("(rabd_|tabd_|te_)", names(df), value = TRUE)
  log2.label = paste0("LOG2.", label)
  df[, log2.label] = sapply(df[, label], 
                            function(x) log2(x) - median(log2(x))) #compute normalized log2
  
  # Average log2 by replicates
  unique.label = unique(str_extract(log2.label, pattern))
  avg.label = paste0("MEAN.", unique.label)
  df[, avg.label] = sapply(unique.label, function(x) {
    replicate.label = grep(x, log2.label, value = TRUE) #vector of replicate labels
    return(apply(df[, replicate.label], 1, mean)) #compute average across replicates
  })
  return(df)
}
DMSO = log2_compute(DMSO)
M1071 = log2_compute(M1071)
Rapa = log2_compute(Rapa)
all = log2_compute(all)


## Split all data by information type (tabd, rabd, te)
rabd = list(raw = cbind(all[1:9], all[grep("rabd_", names(all), value = TRUE)]))
tabd = list(raw = cbind(all[1:9], all[grep("tabd_", names(all), value = TRUE)]))
te = list(raw = cbind(all[1:9], all[grep("te_", names(all), value = TRUE)]))


## LOG2(fold change) for each BIOLOGICAL replicate
log2FC_ratio = function(df, diff = c("M1071", "DMSO"), bioRep = c("-1", "-2", "-3")) {
  # Function only works for pairwise comparisons!
  # df = data frame containing log2 intensity data
  # diff = find the Ratio by computing element1 - element2
  # bioRep = biological replicate identifiers
  log2.names = grep("^LOG2", names(df), value = TRUE)   #Extract all LOG2. names
  ratio.names = paste0("RATIO.", diff[1], "/", diff[2], bioRep)
  df[ratio.names] = lapply(bioRep, function(x) {
    bR.names = grep(x, log2.names, value = TRUE)
    first = grep(diff[1], bR.names, value = TRUE)
    second = grep(diff[2], bR.names, value = TRUE)
    return(df[[first]] - df[[second]])
  })
  return(df)
}
rabd$Rapa.DMSO = log2FC_ratio(rabd$raw, c("Rapa", "DMSO"))
rabd$M1071.DMSO = log2FC_ratio(rabd$raw, c("M1071", "DMSO"))
rabd$Rapa.M1071 = log2FC_ratio(rabd$raw, c("Rapa", "M1071"))
tabd$Rapa.DMSO = log2FC_ratio(tabd$raw, c("Rapa", "DMSO"))
tabd$M1071.DMSO = log2FC_ratio(tabd$raw, c("M1071", "DMSO"))
tabd$Rapa.M1071 = log2FC_ratio(tabd$raw, c("Rapa", "M1071"))
te$Rapa.DMSO = log2FC_ratio(te$raw, c("Rapa", "DMSO"))
te$M1071.DMSO = log2FC_ratio(te$raw, c("M1071", "DMSO"))
te$Rapa.M1071 = log2FC_ratio(te$raw, c("Rapa", "M1071"))


## Welch's two-tailed paired t-test
welch_test = function(df, diff = c("M1071", "DMSO"), bioRep = c("-1", "-2", "-3")) {
  # Function only works for pairwise comparisons
  # df = data frame containing log2 intensity data
  # diff = find the Ratio by computing element1 - element2
  # bioRep = biological replicate identifiers
  log2.names = grep("^LOG2", names(df), value = TRUE)   #Extract all LOG2. names
  group1 = grep(diff[1], log2.names, value = TRUE)
  group2 = grep(diff[2], log2.names, value = TRUE)
  group1 = sapply(bioRep, function(x) grep(x, group1, value = TRUE))  #Orders elements in group1
  group2 = sapply(bioRep, function(x) grep(x, group2, value = TRUE))  #Orders elements in group2
  
  # Perform Welch's two-sided two-sample t-test
  result = data.frame(Pvalue = numeric(0), MEAN.RATIO = numeric(0))
  for (i in 1:nrow(df)) {
    stats = t.test(unlist(df[i, group1]), unlist(df[i, group2]), alternative = "two.sided",
                   paired = TRUE, var.equal = FALSE)
    result[i, ] = c(stats$p.value, stats$estimate)
  }
  result$LOG.Pvalue = -log10(result$Pvalue)
  df = cbind(df, result)
  
  return(df)
}
rabd$Rapa.DMSO = welch_test(rabd$Rapa.DMSO, c("Rapa", "DMSO"))
rabd$M1071.DMSO = welch_test(rabd$M1071.DMSO, c("M1071", "DMSO"))
rabd$Rapa.M1071 = welch_test(rabd$Rapa.M1071, c("Rapa", "M1071"))
tabd$Rapa.DMSO = welch_test(tabd$Rapa.DMSO, c("Rapa", "DMSO"))
tabd$M1071.DMSO = welch_test(tabd$M1071.DMSO, c("M1071", "DMSO"))
tabd$Rapa.M1071 = welch_test(tabd$Rapa.M1071, c("Rapa", "M1071"))
te$Rapa.DMSO = welch_test(te$Rapa.DMSO, c("Rapa", "DMSO"))
te$M1071.DMSO = welch_test(te$M1071.DMSO, c("M1071", "DMSO"))
te$Rapa.M1071 = welch_test(te$Rapa.M1071, c("Rapa", "M1071"))


## Write output files
if(FALSE) {
  print_csv = function(lst, prefix) {
    # lst = list of data frames containing raw data and log2FC
    # prefix = character denoting the output file name PREFIX
    file_names = names(lst)
    invisible(lapply(1:length(file_names), function(i) {
      write.csv(lst[[i]], paste0(prefix, "_", file_names[i], ".csv"), 
                row.names = FALSE)
    }))
  }
  print_csv(rabd, "rabd")
  print_csv(tabd, "tabd")
  print_csv(te, "te")
  write.csv(all, "all_treatment_merged_ribomap_output.csv", row.names = FALSE)
  write.csv(DMSO, "DMSO_replicates_merged_ribomap_output.csv", row.names = FALSE)
  write.csv(M1071, "M1071_replicates_merged_ribomap_output.csv", row.names = FALSE)
  write.csv(Rapa, "Rapa_replicates_merged_ribomap_output.csv", row.names = FALSE)
}


## Save workspace
save.image("ribosome_profiling_workspace.RData")

