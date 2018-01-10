###################################################################
### This script is designed to analyze Ribomap stats output files
###################################################################

# Step 1 - Set working directory and increase memory size for xlsx file writing
#setwd("C:/Set/Working/Directory")



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


## Rename columns with addition of DMSO, Rapa, M1071 tags
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



# Step 4 - Data pre-processing
## Remove data with ribosome footprint < 128 
pro.data = lapply(raw.data, function(x) x[x$rabd >= 128, ])


## Merge all samples into one data frame
merge_treatment = function(pattern, data) {
  # pattern = name of treatment matching to names(data)
  # data = list of data frames containing pre-processed data
  name = grep(pattern, names(data), value = TRUE)
  data_list = data[name]
  df = Reduce(function(df1, df2) full_join(df1, df2, #Full join finds the intersection of the data frames
                                           by = names(df1[1:9])),
              data_list)
  return(df)
}
all = merge_treatment(".*", pro.data)


## Filter for at least 2 valid values in each condition (tailored for ribomap dataset)
filter_NA = function(df, min_count = 2) {
  # df = data frame containing LOG2 data for filtering
  # min_count = minimum number of quantified replicates for each condition for retention
  conditions = c("DMSO", "M1071", "Rapamycin")
  log.names = grep("^rabd", names(df), value = TRUE)   #extract columns starting with rabd (could be tabd or te because either all three info is available or none at all)
  replicates = lapply(conditions, function(x) grep(x, log.names, value = TRUE)) #group technical replicates
  rep.filter = sapply(replicates, function(y) rowSums(!is.na(df[y]))) #count number quantified for each condition
  keep = apply(rep.filter, 1, function(x) all(x >= min_count))
  
  ####Uncomment to require all valid values in DMSO replicates
  #df = df[keep, ]
  #x = grep("DMSO", log.names, value = TRUE)
  #keep = rowSums(!is.na(df[x])) == 3
  
  return(df[keep, ])   #filter rows
}
all.filtered = filter_NA(all)


## Calculate normalized log2 and averages of ribomap statistics
log2_compute = function(df, pattern = ".*(?=-[[:digit:]])") {
  # df = data frame containing rabd, tabd, te data
  # pattern = regex expression describing how to extract the replicates
  label = grep("(rabd_|tabd_|te_)", names(df), value = TRUE)
  log2.label = paste0("LOG2.", label)
  df[, log2.label] = sapply(df[, label], #compute normalized log2
                            function(x) log2(x) - median(log2(x), na.rm = TRUE)) 
  
  # Average log2 by replicates
  unique.label = unique(str_extract(log2.label, pattern))
  avg.label = paste0("MEAN.", unique.label)
  df[, avg.label] = sapply(unique.label, function(x) {
    replicate.label = grep(x, log2.label, value = TRUE) #vector of replicate labels
    return(apply(df[, replicate.label], 1, function(x) mean(x, na.rm = TRUE))) #compute average across replicates
  })
  return(df)
}
all.filtered = log2_compute(all.filtered)


## Split all data by information type (tabd, rabd, te)
rabd = list(raw = cbind(all.filtered[1:9], 
                        all.filtered[grep("rabd_", names(all.filtered), value = TRUE)]))
tabd = list(raw = cbind(all.filtered[1:9], 
                        all.filtered[grep("tabd_", names(all.filtered), value = TRUE)]))
te = list(raw = cbind(all.filtered[1:9], 
                      all.filtered[grep("te_", names(all.filtered), value = TRUE)]))


## Welch's two-tailed unpaired t-test
welch_test = function(df, diff = c("M1071", "DMSO")) {
  # Function only works for pairwise comparisons
  # df = data frame containing log2 intensity data
  # diff = find the Ratio by computing element1 - element2
  log2.names = grep("^LOG2", names(df), value = TRUE)   #Extract all LOG2. names
  group1 = grep(diff[1], log2.names, value = TRUE)
  group2 = grep(diff[2], log2.names, value = TRUE)
  
  # Perform Welch's two-sided two-sample t-test
  result = data.frame(Pvalue = numeric(0), LOG2.RATIO = numeric(0))
  for (i in 1:nrow(df)) {
    stats = t.test(unlist(df[i, group1]), unlist(df[i, group2]), alternative = "two.sided",
                   paired = FALSE, var.equal = FALSE)   #var.equal = FALSE means Welch's test
    #result[i, ] = c(stats$p.value, stats$estimate) #for paired
    result[i, ] = c(stats$p.value, stats$estimate[1] - stats$estimate[2]) #for unpaired
  }
  result$FDR = p.adjust(result$Pvalue, method = "fdr")   #ADJUST P-VALUE USING BENJAMINI-HOCHBERG FDR
  result$LOG.Pvalue = -log10(result$Pvalue)
  df = cbind(df, result)
  
  return(df)
}
rabd$Rapa.DMSO = welch_test(rabd$raw, c("Rapa", "DMSO"))
rabd$M1071.DMSO = welch_test(rabd$raw, c("M1071", "DMSO"))
rabd$Rapa.M1071 = welch_test(rabd$raw, c("Rapa", "M1071"))
tabd$Rapa.DMSO = welch_test(tabd$raw, c("Rapa", "DMSO"))
tabd$M1071.DMSO = welch_test(tabd$raw, c("M1071", "DMSO"))
tabd$Rapa.M1071 = welch_test(tabd$raw, c("Rapa", "M1071"))
te$Rapa.DMSO = welch_test(te$raw, c("Rapa", "DMSO"))
te$M1071.DMSO = welch_test(te$raw, c("M1071", "DMSO"))
te$Rapa.M1071 = welch_test(te$raw, c("Rapa", "M1071"))


## Write output files
if(TRUE) {
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
  write.csv(all.filtered, "all_treatment_merged_ribomap_output.csv", row.names = FALSE)
}


## Save workspace
save.image("ribosome_profiling_workspace.RData")
