# extract credentials of the app owner
extract_credentials <- function(path=getwd()){
  lastElement <- function(Avector){
    Avector[length(Avector)]
  }
  AppID <- unlist(lapply(strsplit(path,"/"),FUN = lastElement))
  if (Universal.code) {
    AppOwerCredentials <- readRDS("/srv/shiny-server/DataHub/initiated.account.and.password.rds")[["Universal"]][,c(1:6)]
  }else{
    AppOwerCredentials <- readRDS("/srv/shiny-server/DataHub/initiated.account.and.password.rds")[[AppID]][,c(1:6)]
  }
  AdminCredentials <- readRDS("/srv/shiny-server/DataHub/initiated.account.and.password.rds")[["Admin"]][,c(1:6)]
  return(rbind(AppOwerCredentials,AdminCredentials))
}

# check data.frame all columns levels, return a logic vector
check_df_level <- function(df, cutoff = ifelse(exists("SplitBy.levels.max"), SplitBy.levels.max, 15)){
  # as.factor columns in meta.data with levels less than defined by 'SplitBy.levels.max'
  unname(sapply(df,function(x){length(unique(as.character(x)))}))  <= cutoff
}
# check if data.frame all columns are factor, return a logic vector
check_df_factor <- function(df){
  unname(sapply(df,is.factor))
}

color.vector <- colorRampPalette(c("#638ED0", "#ffff33","#ff3300"))(30) # color used in featureplot


# define function - Chcek the input gene, return the revised gene
CheckGene <- function(InputGene, GeneLibrary){
  InputGene <- trimws(InputGene)
  if (InputGene %in% GeneLibrary) { # when input gene is absolutely right
    return(InputGene)
  }else if(tolower(InputGene) %in% tolower(GeneLibrary)){ # when not case sensitive
    return(GeneLibrary[tolower(GeneLibrary) %in% tolower(InputGene)]) # gene list length can > 1
  }else{ # when not match
    return(NA)
  }
}