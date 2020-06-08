### Altered function .get_sum_comp to enable kendrick mass defect calculations. 
### Code now requires the input of m/z values for all input names.
### Input data is no longer filtered for unique names, but each row is treated as a separate feature.
### Adducts and lipid classes are extracted. 
.get_sum_comp_kmd <- function(names, mzinput) {
  .en <- function(...) paste0("(", ..., ")")
  csep = "[/-_]"
  dbond_config = "\\((?:\\d{1,2}[ZE][,]*)+\\)"
  chain_notes = "\\((?:\\d{1,2}[^)]*)+\\)"
  chain_p = paste0("([dthOP]?(?:methyl)?-?)(\\d{1,2}):(\\d{1,2})(",.en("?:", chain_notes), "*)")
  chain_multi_p = paste0(.en("?:", chain_p, csep, "?"), "+")
  
  names <- names
  mzinput <- mzinput
  fas <- tibble(names = names, mz = mzinput, match=str_extract(names, chain_multi_p))
  sum_comp <- fas %>% filter(!is.na(match)) %>%
    mutate(chains = str_match_all(match, chain_p)) %>%
    group_by(match) %>%
    mutate(
      links = .collapse_char(chains[[1]][,2]),
      total_cl = sum(as.numeric(chains[[1]][,3])),
      total_cs = sum(as.numeric(chains[[1]][,4])),
      modifs = .collapse_char(chains[[1]][,5]),
      odd_chain = any(as.numeric(chains[[1]][,3]) %% 2 > 0)
    ) %>% 
    ungroup() %>%
    #### Added to extract adducts from input names > Seems to work on all rows in results table from Shiny App. 
    select(-chains)  %>% 
    mutate(adduct = sub(".*)", "", names), class = sub("\\(.*", "", names)) %>%
    mutate(
      links = sub(",+$", "", links), modifs = sub(",+$", "", modifs),
      sum_composition = paste0(links, total_cl, ":", total_cs, modifs)
    ) 
  
}


### New function to calculate kendrick mass defects for input data. 
### This function directly accepts the output from .get_sum_comp_kmd
### This function is required to create the ref_kmd_table below
.get_kmd <- function(mol_info) {
  kendrick_mass <- mol_info$mz * (14 / 14.01565)
  kendrick_mass_defect <- kendrick_mass - as.integer(kendrick_mass)
  
  data.frame(mol_info$names, mol_info$mz, kendrick_mass, kendrick_mass_defect, mol_info$class, mol_info$total_cl, mol_info$total_cs, mol_info$adduct)
}


### This function accepts the output of .get_kmd to calculate the difference between the chain saturation and expected chain saturation
### Load required table containing reference kendrick mass defect per lipid class, per adduct. (in LipID/data)
### Add reference kendrick mass defects, calculate difference from measured KMD, calculate multiples of -0.013399 (error of unsaturation), calculate difference expected cs and cs.
.get_cs_residuals <- function(kmd_df) {
  load("../data/ref_kmd_table.rda")
  
  left_join(kmd_df, ref_kmd_table) %>% mutate(cl_expected = (kendrick_mass_defect - meankmd) / -0.013399) %>% mutate(cs_residual = abs(mol_info.total_cs - cl_expected))
  
}


### I created a reference kendrick mass defect table for all lipid classes / adduct combinations in the LipidMatch_Libraries folder 
### A list of files is created and used to loop over all individual files.
### Some lipid classes will produce an error, as no chain length / unsaturation information is present. tryCatch is used to ingore those.
### The meankmd refers to the kendrick mass defect of a lipid class would have if it is completely saturated.
### This code only needs to be run once to produce the library, but is saved in case of future releases of new databases.
libfilelist <- list.files("../data-raw/LipidMatch_Libraries/", full.names = TRUE)
kmdlib <- list()

for(i in 1:length(libfilelist)){
  libfile <- read.csv(libfilelist[i])
  print(libfilelist[i])
  tryCatch({
  temp <-  .get_kmd(mol_info = .get_sum_comp_kmd(libfile[,1], libfile[,2])) %>% mutate(classkmd = kendrick_mass_defect + (mol_info.total_cs * 0.013399)) %>% group_by(mol_info.class, mol_info.adduct) %>% summarize(meankmd = mean(classkmd)) %>% as.data.frame()
  kmdlib[[i]] <- temp
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ref_kmd_table <- plyr::compact(kmdlib) %>% 
                do.call(rbind, .) %>%
                filter(is.na(meankmd) == FALSE) %>% 
                group_by(mol_info.class, mol_info.adduct) %>% 
                summarize(meankmd = mean(meankmd))

save(ref_kmd_table, file = "ref_kmd_table.rda")




### Test the functions using data
testdata <- read.csv("../data/annotated_features.rda", header = TRUE, stringsAsFactors = FALSE)

testdata_kmd <- .get_sum_comp_kmd(testdata$name, testdata$mz) %>% 
                .get_kmd() %>% 
                .get_cs_residuals()


### If needed, this testdata_kmd table can be joined with the original input data
testdata <- left_join(testdata, testdata_kmd, by = c("mz" = "mol_info.mz", "name" = "mol_info.names")) %>% distinct()

### Visualize the cs_residuals for the lipids marked as best_match of not (We expect lower residuals for the best_match lipids)
ggplot(testdata_plot, aes(x = mol_info.class, y = cs_residual)) + geom_boxplot() + facet_grid(rows = vars(best_match))


