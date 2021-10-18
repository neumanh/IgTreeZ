update_region_object <- function(cdr3_end, region = IMGT_V) {
  # Updates the IMGT_V region and adds the CDR3 region to the region definition
  
  udpt_reg <- region
  
  if (!(is.na(cdr3_end)) && (udpt_reg@seqLength < cdr3_end)){
    
    # Get the CDR3 length
    cdr3_len <- cdr3_end - udpt_reg@seqLength
    
    # Update the seqLength slot
    udpt_reg@seqLength <- cdr3_end 
    
    # Defing new factor for the CDR3 region
    short_bndrs <- factor(rep('cdr', cdr3_len), levels = c('cdr', 'fwr'))
    
    # Concatinating the two factors
    long_bndrs <- unlist(list(udpt_reg@boundaries, short_bndrs)) 
    
    
    # Update the boundaries slot
    udpt_reg@boundaries <- long_bndrs
    
    # Update the description slot
    udpt_reg@description <- paste0(udpt_reg@description, ' including CDR3 region.')
    
  } else {
    if (is.na(cdr3_end))
      cat('CDR3 end column contains NA.\n')
    else 
      cat('CDR3 end position is shorter than FWR3 end position. No change was done.\n')
  }
  
  return(udpt_reg)
}


calc_expected_on_one_row <- function(dfr, 
                                     sequenceColumn="clonal_sequence",
                                     germlineColumn="clonal_germline",
                                     targetingModel=HH_S5F) {
  # Updates the region definition by the CDR3 length, and runs expectedMutations
  cdr3_end = as.integer(dfr[['cdr3_end']])
  #cat('****** cdr3_end   ', cdr3_end, ' ********\n')
  #cat('****** GL length  ', nchar(as.character(dfr[[germlineColumn]][1])), ' ********\n')
  
  updt_reg = update_region_object(cdr3_end = cdr3_end)
  
  dfr <- dfr
  sink("/dev/null") 
  #sink("NUL")
  
  exp = expectedMutations(dfr,
                    sequenceColumn=sequenceColumn,
                    germlineColumn=germlineColumn,
                    targetingModel=HH_S5F,
                    regionDefinition=updt_reg)
  sink()  
  
  return(exp)
}


update_col_names <- function(filename) {
  # For the new shazam version
  
  # TODO: add check if columns exist
  df = read.csv(filename)
  
  df_names = names(df)
  
  if ('clone_id' %in% df_names)
    df_names = replace(df_names, df_names=='clone_id', 'original_clone_id')
  
  if ('GERMLINE_IMGT_D_MASK' %in% df_names)
    df_names = replace(df_names, df_names=='GERMLINE_IMGT_D_MASK', 'clonal_germline')
  
  if ('SAMPLE' %in% df_names)
    df_names = replace(df_names, df_names=='SAMPLE', 'sample')
  
  
  names(df) =  df_names
  
  return(df)
}

