# Concatenating External Data Matrices
# Allometry & Sexual Selection Project MSU
# Emerald Bender Fall 2024/Spring 2025
# this file defines functions used in merging data:
# 'add_taxonomy', 'df_merge', 'column_merge'

# define function for adding taxonomy
# function 'add.taxonomy' requires an input dataframe with species name in column 'binomial'
# output is a dataframe with standardized species names found by rotl in column 'binomial_rotl'
# and higher taxonomy in columns 'genus', 'family', 'order', and 'class'
# outputs warnings about non-unique column names when making temp3
# ^ these are not a problem! disregard

add_taxonomy<- function(df){  # this function does not work if there are duplicated binomials that rotl can't id
  
  temp.names<- as.vector(df$binomial) # requires column named binomial in df input
  taxlevels<- c('genus', 'family', 'order', 'class') # which taxonomy columns do I retrieve? 
  
  niter<- ceiling(length(temp.names)/2500) 
  loops <-list() 
  
  for (i in 1:niter) { # open niter.loop; since rotl is limited on number of matches, loop over max name vec length
    print(paste("loop", i, "of", niter)) # print loop number
    
    if (length(temp.names)> 2500) { # open names.ifelse
      names.num<- c(1, 2500) + (i-1)*2500
      if (names.num[2]> length(temp.names)) { # open names.ifelse.2
        temp<- tnrs_match_names(names=temp.names[c(names.num[1]:length(temp.names))])
      } else {
        temp<- tnrs_match_names(names=temp.names[c(names.num[1]:names.num[2])])
      } # close names.ifelse.2
    } else {temp<- tnrs_match_names(names=temp.names)} # close names.ifelse
    
    print("ok") # checkpoint 1
    
    NAs<-which(is.na(temp$ott_id))
    NAs.names<- temp$search_string[NAs]
    print(paste(length(NAs), "NAs:", paste(NAs.names, collapse = "; ")))
    if (length(NAs)>0) {
      temp2<- temp[-NAs,] # removing rows that didn't match any names in rotl
    } else temp2<- temp
    
    temp3<- tax_lineage(taxonomy_taxon_info(temp2$ott_id, include_lineage=TRUE)) # finding higher taxonomy info
    print("ok2") # checkpoint 2
    
    temp4<- lapply(temp3, function(x){ # reformatting taxonomy data to df with columns in 'taxlist'
      x<- t(x)
      x<- data.frame(row_to_names(x[c(1:2),], row_number=1))
      x2<- x[intersect(taxlevels, names(x))] # calling taxonomy columns by name
      x2[setdiff(taxlevels, names(x))]<- NA # fill missing columns with NA
      x2<- x2[,c(taxlevels)] # reorder columns to ascending tax level
      colnames(x2)<- paste(taxlevels, "rotl", sep= "_")
      x2
    }
    )
    print("ok3") # checkpoint 3
    
    temp5<- data.frame(do.call(rbind, temp4))%>% mutate( # making list temp3 into dataframe
      binomial_rotl = as.vector(temp2$unique_name)
    )
    print("ok4") # checkpoint 4
    
    # duplicated rows containing taxonomy data are identical and can be safely removed at this point
    temp6<- right_join(distinct(temp5), temp, join_by( binomial_rotl ==  unique_name)) # adding searched name back into df
    loops[[i]]<- temp6
    
  } # close niter.loop
  
  final.tax<- do.call(rbind, loops) %>% mutate( # putting dfs from loop together
    search_str_sentence= str_to_sentence(search_string)
  )
  
  final.tax2<- left_join(data.frame(temp.names), final.tax, join_by(temp.names == search_str_sentence)) %>% rename(
    data_binom_orig = temp.names # merging original name string list with the search strings in rotl output
  ) # so that the taxonomy data can be merged with the trait data easily 
  
  print("ok4") #checkpoint 4
  
  final<- left_join(df, final.tax2, join_by(binomial == data_binom_orig), multiple= "first") 
  
  final<- subset(final, select= -c(search_string)) # removing redundant column
  # rotl doesn't use 'reptilia' as a class
  # but I want to...
  reptile.orders<- c("Squamata", "Testudines", "Crocodylia", "Rhynchocephalia")
  final$class_rotl[final$order_rotl %in% reptile.orders]<- "Reptilia"
  
  return(final)
} # close function add.taxonomy


#######
#######
#######
# function to merge dataframes after they have been read in and cleaned
#   and run through add_taxonomy
# feed in list of databases to be merged into single df
# this also tags column names with the database that they come from
# MUST have column 'binomial_rotl' in every input df, used as key to merge
df_merge<- function(db.list){
  
  remove.cols<- c("approximate_match", "binomial", 
                  "Class", 'Common_name','common_name',
                   "family", "Family","flags","Genus", 
                  "is_synonym", "intercept",
                   "number_matches", "Order", "ott_id", 
                  "References", "score", "slope", "Species")
  names.db.list<- names(db.list)
  names.db.list<- gsub(".edit.tax", "", names.db.list)
  
  # homogenizing naming conventions for my sanity... 
  db.list<- lapply(db.list, function(x){
    colnames(x) <- gsub( ".", "_", colnames(x), fixed=TRUE)
    x$binomial <- gsub( " ", "_", x$binomial)
    x<- x[!is.na(x$binomial_rotl),]
    x
  }
  )
  
  db.list<- lapply(db.list, function(x)
    x<- x[,-which(names(x) %in% remove.cols)])
  all_colnames<- unlist(sapply(db.list[c(2:length(db.list))], function(x) colnames(x)))
  
  # finding columns unique to first db
  # to tag them with db name below
  db1_unique<- setdiff(colnames(db.list[[1]]), all_colnames)
  #print(db1_unique)
  
  for (i in 1:length(db.list)) {
    
    if (i==1) {
      temp.df <- db.list[[1]]
      temp.df<- setNames(temp.df, paste0(
        names(temp.df),
        ifelse(names(temp.df) %in% db1_unique, paste0(".",names.db.list[i]), "")
      ))
      
    } else{ 
      temp.df.prev<- temp.df
      #shared_col<- intersect(names(temp.df.prev), names(db.list[[i]]))
      #print(shared_col)
      
      temp.df<- merge(temp.df.prev, db.list[[i]], by = c("binomial_rotl", "genus_rotl", "family_rotl", "class_rotl", "order_rotl"),
                      all.x=T, all.y=T, suffixes = c("", names.db.list[i]))
      setdiff (names(db.list[[i]]), names(temp.df.prev))
      temp.df<- setNames(temp.df, paste0(
        names(temp.df), 
        ifelse(names(temp.df) %in% setdiff(names(db.list[[i]]), names(temp.df.prev)), paste0(".", names.db.list[i]), "")
      ))
    }
  }
  return(temp.df)
} # close function 'df.merge'

#######
#######

# If there are multiple values for one species, unique values are placed in new column
# result is multiple columns for continuous traits if there are multiple measurements for a single species
# The if genus_rotl statement in the middle allows removal of duplicated continuous trait values
#     Do NOT apply genus_rotl section to categorical traits

column_merge<- function(dat, datname, cols=NULL){
  
  temp.colnames<- colnames(dat)
  temp.colnames.cont<- temp.colnames[!temp.colnames %in% c('binomial_rotl', 'genus_rotl')]
  
  #identifying duplicated continuous traits
  if("genus_rotl" %in% colnames(dat)){
    
    dat.format<- data.frame(sapply(dat, as.character))
    temp.pivot<- pivot_longer(dat.format, cols=all_of(temp.colnames.cont))
    temp.pivot <- na.omit(temp.pivot)
    temp.pivot2<- which(duplicated(temp.pivot[c('genus_rotl', 'value')])|duplicated(temp.pivot[c('genus_rotl', 'value')], fromLast = TRUE))
    temp.pivot3<- temp.pivot[temp.pivot2,]
    temp.pivot4 <- which(duplicated(temp.pivot3[c('binomial_rotl', 'value')])) #|duplicated(temp.pivot3$binomial_rotl, fromLast = TRUE))
    temp.pivot5<- temp.pivot3[-temp.pivot4,]
    temp.pivot6<- which(duplicated(temp.pivot5[c('genus_rotl', 'value')])|duplicated(temp.pivot5[c('genus_rotl', 'value')], fromLast = TRUE))
    temp.pivot.final<- temp.pivot5[temp.pivot6,]
    print("ok")
    
    dat$genus_rotl<- NULL
  }
  
  
  temp<- lapply(split(dat[-1], dat$binomial_rotl), function(x) x[!is.na(x)])
  
  
  zeros<- which(lengths(temp)==0)
  temp[zeros]<- "PLACEHOLDER"
  
  df <- as.data.frame(do.call(rbind, temp))
  
  df$V2[which(df$V1 == df$V2)]<- NA
  
  if (ncol(df)==3) {
    
    df$V3[which(df$V1 == df$V3)]<- NA
    NAs<- which(is.na(df$V2) & !is.na(df$V3))
    df$V2[NAs]<- df$V3[NAs]
    df$V3[which(df$V2 == df$V3)]<- NA
  }
  print(paste(datname, seq(from=1, to=(ncol(df))), sep="_"))
  colnames(df)<- paste(datname, seq(from=1, to=(ncol(df))), sep="_")
  print("ok")
  df[df == "PLACEHOLDER"]<- NA
  
  if('genus_rotl' %in% temp.colnames){
    df[rownames(df) %in% temp.pivot.final$binomial_rotl,]<-NA
  }
  
  df
} # close function 'column_merge'
