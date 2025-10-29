# Install CellChatDB if you haven't already
# devtools::install_github("sqjin/CellChatDB")

# Load libraries
library(CellChat)
library(dplyr)
library(tidyr)

# 1. Load the CellChatDB tables
interactions_df <- CellChatDB.human$interaction
complex_df <- CellChatDB.human$complex

# 2. Helper function to get subunits for any name
# It checks the complex_df. If the name is not there, it's a single gene.
get_subunits <- function(name, complex_db) {
  if (name %in% rownames(complex_db)) {
    # It's a complex. Get its subunits from the columns.
    subunits <- unlist(complex_db[name, ])
    # Filter out any NAs or empty strings
    subunits <- subunits[!is.na(subunits) & subunits != ""]
    return(as.character(subunits))
  } else {
    # It's a single gene. Return the name itself as a vector.
    return(as.character(name))
  }
}

# 3. Create an empty list to store all 1-to-1 pairs
all_subunit_pairs <- list()

# 4. Iterate through each original interaction
for (i in 1:nrow(interactions_df)) {
  interaction_row <- interactions_df[i, ]
  
  ligand_name <- interaction_row$ligand
  receptor_name <- interaction_row$receptor
  
  # Get all single-gene subunits for the ligand
  ligand_subunits <- get_subunits(ligand_name, complex_df)
  
  # Get all single-gene subunits for the receptor
  receptor_subunits <- get_subunits(receptor_name, complex_df)
  
  # 5. Create all possible 1-to-1 combinations (all-vs-all)
  # expand.grid is perfect for this
  subunit_grid <- expand.grid(ligand_subunit = ligand_subunits, 
                              receptor_subunit = receptor_subunits, 
                              stringsAsFactors = FALSE)
  
  all_subunit_pairs[[i]] <- subunit_grid
}

# 6. Combine all pairs into one big data frame
final_pairs_df <- do.call(rbind, all_subunit_pairs)

# 7. Clean up, remove duplicates, and rename columns
final_pairs_df <- final_pairs_df %>% 
  dplyr::distinct() %>%
  dplyr::rename(ligand = ligand_subunit, receptor = receptor_subunit)

# 8. Save the final file
write.table(
  final_pairs_df,
  file = "results/cellchat_rl_pairs.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

print(paste("Successfully created cellchat_rl_pairs.tsv with", nrow(final_pairs_df), "1-to-1 gene pairs."))
