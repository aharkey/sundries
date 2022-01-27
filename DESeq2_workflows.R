# Andria Harkey

# Functions and workflow for simplifying results from multiple DESeq2 comparisons




### Functions #####

# Takes output dataframe from combined DESeq2 results (a) and:
  # removes lfcSE columns
  # removes all but first baseMean column
  # and changes log2FoldChange column name to logFC

  # txt (text) option will remove any additional text from all column labels
  # containing it (such as a prefix or suffix)
simpleout <- function(a, txt)
{
  b <- a[,-c(grep("lfcSE", colnames(a)))]
  colnames(b)[grep("baseMean", colnames(b))][1] <- "baseMean"
  b <- b[,-grep("baseMean", colnames(b))[-1]]
  colnames(b) <- gsub("log2FoldChange", "logFC", colnames(b))
  if(missing(txt))
  {}
  else
  {colnames(b) <- gsub(txt, "", colnames(b))}
  return(b)
}

# Useful functions for summarizing what genes are doing across multiple samples:
# Given logFC & pvalue and cutoffs for both values, output is U for genes with significant positive logFC, D for genes with significant negative logFC, and n for genes with non-significant results
upordown <- function(pval, lfc, pcut, lfccut)
{if(pval < pcut & abs(lfc) > lfccut)
{if(lfc > 0)
{out <- "U"}
  else
  {out <- "D"}}
  else
  {out <- "n"}
  return(out)}

# Given a range of values, which value has the largest absolute value
maxabs <- function(a)
{ifelse(max(a)>abs(min(a)), max(a), min(a))}








### Work flow #####

# Given a set of data frames (df1, df2, df3, etc.) with DESeq2 results

# Combine all pairwise comparison results into one dataframe
## Give columns unique names
colnames(df.1) <- paste0("df.1.", colnames(df.1))
colnames(df.1) <- paste0("df.2.", colnames(df.2))
colnames(df.3) <- paste0("df.2.", colnames(df.3))

## Combine by gene IDs
df.1$genes <- rownames(df.1)
df.2$genes <- rownames(df.2)
df.3$genes <- rownames(df.3)

df.combine <- plyr::join_all(lapply(list(df.1,
                                         df.2,
                                         df.3), as.data.frame), by = "genes")

rownames(df.combine) <- df.combine$genes
df.combine <- df.combine[,-6]

df.combine <- simpleout(df.combine, "df.")








### IQR workflow #####

# For removing low IQR genes from counts (cts) table
# BEFORE running DESeq2

# Set IQR cutoff
iqrcut <- 5

# Calculate IQRs
iqrs <- NA
for(i in 1:nrow(cts))
{
  iqrs[i] <- IQR(cts[i,])
}

# Remove genes with IQR less than 5
cts <- cts[which(iqrs > iqrcut),]