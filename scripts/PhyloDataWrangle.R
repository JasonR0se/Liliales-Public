### Importing Libraries ###
library("phytools")
library('tidyverse')
library('stringr')

outputdirectory = PhyloDataWrangle
# Directs to Reserach Folder #
setwd("C:/Users/ijaso/Desktop/Research/Jason/Liliales")

## Reading in Data ##
raw_Data <- read.csv("data/TraitDataCollection/NewDataCollection.csv",na.strings=c(""," ","NA"))
tree <- read.tree("data/big_constraint_fixed_dated_Liliales.tre")

### How Many Unique References? - What does this mean again?
cat(length( unique(raw_Data$Reference.1)),": unique references!" )
print("\n")

### Species Names and Character State 
clean_Data <- data.frame(raw_Data$Species..in.Phylo., #What's the formatting for this?
                         raw_Data$InfloresenceCoded) %>%
  drop_na() 

# Print Characer States Facts
unique_names <- as.character(unique(clean_Data$raw_Data.InfloresenceCoded))
cat("Number of Character States:", length(unique_names),
    "\nCharacter States:", unique_names )
print("\n")

#Checking for duplicates
allDup <- function (value) {
  duplicated(value) | duplicated(value, fromLast = TRUE)
}

clean_Data[allDup(clean_Data$raw_Data.Species..in.Phylo.),]

#Duplicate Removal - I had no duplicates so I commented this out.
clean_Data_duplicatesRemoved <- clean_Data[!duplicated(clean_Data), ]
conflict <- clean_Data_duplicatesRemoved[allDup(clean_Data_duplicatesRemoved$raw_Data.Species..in.Phylo.),]
print("These Species have conflicting Scoring")
(conflict) #Clean up by using the cat funtion. See above for examples. Maybe wirte a for loop?

clean_Data_Final <- clean_Data_duplicatesRemoved[!(clean_Data_duplicatesRemoved$raw_Data.Species..in.Phylo. %in%
                                                     conflict$raw_Data.Species..in.Phylo. ),]


########## NOTE: THIS REMOVE ROW THAT HAVE TAXA SCORED AS SINGLE - How exactly does this work?
#clean_Data_Final <- clean_Data_Final[ ! clean_Data_Final$raw_Data.InfloresenceCoded %in% "single", ]

# Remove any species in the phylogeny data not in database
tips2keep <- as.character(clean_Data_Final$raw_Data.Species..in.Phylo.)

#Generates an error. (Returing Null) but does not stop me. - It stops me now.
trimtree = drop.tip(tree, setdiff(tree$tip.label, tips2keep))

clean_Data_TreeMatch <- clean_Data_Final[!(clean_Data_Final$raw_Data.Species..in.Phylo. %in%
                                             trimtree$tip.label ),]
clean_Data_ForNexus <- clean_Data_Final[!(clean_Data_Final$raw_Data.Species..in.Phylo. %in%
                                            clean_Data_TreeMatch$raw_Data.Species..in.Phylo. ),]

#Nexus Formatting - Narrowed down my problem but I still don't know how to fix it really.
for_Nexus <- transform(clean_Data_ForNexus, id = as.numeric(factor(clean_Data_ForNexus$raw_Data.InfloresenceCoded))) %>%
  t() %>%
  as.data.frame()

# Put colnames back
names(for_Nexus) <- lapply(for_Nexus[1, ], as.character)
for_Nexus   <- for_Nexus[-1, ] # Removes the numbers at the top
nexus_Ready <- for_Nexus[-1, ] # Removes the infloresnces


# Write csv data
write.csv(t(nexus_Ready), file = "data/TraitDataCollection/output/NewCSVFile.csv",quote = FALSE)

# write nexus data
write.nexus.data(nexus_Ready, file ="data/TraitDataCollection/output/cleanTraitDataContiouse.nex", format = "standard")

# The follwing removes the values 0,6,7,8,9  from that section
#tx  <- readLines("Jason/data/TraitDataCollection/output/cleanTraitDataContiouse.nex")
#tx2  <- gsub(pattern = "\".*?1", replace = "\"1", x = tx)
#tx3  <- gsub(pattern = "4.*?\"", replace = "4\"", x = tx2)
#writeLines(tx3, con="Jason/dataTraitDataCollection/output/cleanTraitDataContiouse.nex")

# Write trimed Phylogeny and Data .tree and .nexus 
write.tree(trimtree, file = "TraitDataCollection/output/trimed_big_constraint_fixed_dated.tre")

