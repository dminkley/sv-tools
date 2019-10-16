
# Preliminary setup ------------------------------------------------------------------

# Load libraries
library(rbamtools)  # Note that rbamtools is 0-based
library(GenomicAlignments)
library(GenomicRanges)
library(IRanges)
library(reshape)
library(scales)
library(svMisc)
library(futile.logger)

# Get arguments from the command line
flog.info("Processing arguments from the command line.")
args <- commandArgs(trailingOnly=TRUE)
flog.info("%s arguments detected:", length(args))
args_cwd <- args[1]
args_read_length_distro_fn <- args[2]
args_sv_fn <- args[3]
args_sample_to_bamfile_map_fn <- args[4]
args_readpair_info_fn <- args[5]
args_flank_size <- as.numeric(args[6])


#flog.warn("Overriding any command-line arguments with hard-coded arguments!")
#args_cwd <- "/scratch/dminkley/omyk_paper_results/step11_structural_variant_filtering/step09RERUN_visualization/step02_final_sv_set_visualization/all_svs"
#args_read_length_distro_fn <- "omyk_full_set_no_gmo_1_read_fragment_length_info.txt"
#args_sv_fn <- "omyk_full_set_no_gmo_1.oct_2019_filtered.formatted_for_R_vis.txt"
#args_sample_to_bamfile_map_fn <- "sample_name_to_bam_file_mapping.no_Ref_Ind.txt"
#args_readpair_info_fn <- "test_readpair_info_file.txt"
#args_flank_size <- 1000

# Report the used argument values
flog.info("The following arguments are being used:")
flog.info("\targs_cwd: %s", args_cwd)
flog.info("\targs_read_length_distro_fn: %s", args_read_length_distro_fn)
flog.info("\targs_sv_fn: %s", args_sv_fn)
flog.info("\targs_sample_to_bamfile_map_fn: %s", args_sample_to_bamfile_map_fn)
flog.info("\targs_readpair_info_fn: %s", args_readpair_info_fn)
flog.info("\targs_flank_size: %s", args_flank_size)

# Set base directory
setwd(args_cwd)

# Load read length distribution info
flog.info("Reading read length distribution information file")
read_len_distro_col_names <- c("sample_id", "bam_file", "median", "perc_90", "perc_95", "perc_99", "perc_99.5", "perc_99.9", "perc_99.99")
read_len_distro_info.df <- read.table(args_read_length_distro_fn, col.names=read_len_distro_col_names)


# Load in variants to be visualized.
flog.info("Reading structural variant information file")
flog.debug("IMPORTANT: check.names is off when reading the sv info file; subtracting could be dangerous")
vars.df <- read.table(args_sv_fn, header = TRUE, check.names=FALSE)

# Get a mapping of sample_ids to bamfiles
flog.info("Reading in sample_id to bam file map")
samples_and_fns.df <- read.table(args_sample_to_bamfile_map_fn, header=FALSE, stringsAsFactors=FALSE,
                                 col.names=c("sample_id", "filename"))
samples_and_fns.df$basename <- basename(samples_and_fns.df$filename)

# Sanity check that vars.df ind_names are in the same order as the samples in samples_and_fns.df
if (!identical(as.character(samples_and_fns.df$sample_id), as.character(colnames(vars.df[10:ncol(vars.df)])))) {
  flog.fatal("The sample ids in the sample_id to bam file map are either in a different order or different altogether from those in SV info file")
  stop("A fatal error occurred, see log message for details")
}

# Function definitions ----------------------------------------------------
flog.trace("Loading function definitions")

# Given a rbamtools cigarData dataframe, return the original CIGAR string
getCigarString <- function(cigar_data) {
  cigar_string <- ""
  for (i in 1:nrow(cigar_data)) {
    cigar_string <- paste(cigar_string, cigar_data[i,]$Length, cigar_data[i,]$Type, sep="")
  }
  return(cigar_string)
}

# Function to determine the 'type' for a given read pair
# Note that internally, this converts to 1-based coordinates (from 0-based used by rbamtools)
get_readpair_type <- function(good_read_pair, max_normal_insert_size) {
  
  # Read one will always be the lowest of the two reads
  if (rbamtools::position(good_read_pair[[1]]) <= rbamtools::position(good_read_pair[[2]])) {
    r1 <- good_read_pair[[1]]
    r2 <- good_read_pair[[2]]    
  } else {
    r1 <- good_read_pair[[2]]
    r2 <- good_read_pair[[1]]   
  }
  
  # Get read 1 values
  r1_start <- position(r1) + 1
  r1_end <- position(r1) + cigarWidthAlongReferenceSpace(getCigarString(cigarData(r1)))
  if (reverseStrand(r1)) 
    r1_strand <- "-"
  else
    r1_strand <- "+"
  
  # Get read 2 values
  r2_start <- position(r2) + 1
  r2_end <- position(r2) + cigarWidthAlongReferenceSpace(getCigarString(cigarData(r2)))
  if (reverseStrand(r2)) 
    r2_strand <- "-"
  else
    r2_strand <- "+"  
  
  # Assign type to the read pair
  if (r1_strand == "+" && r2_strand == "-") {
    if ((r2_end - r1_start) <= max_normal_insert_size) {
      type <- "normal"
    } else {
      type <- "long"
    }
  } else if (r1_strand == "+" && r2_strand == "+") {
    type <- "both_positive"
  } else if (r1_strand == "-" && r2_strand == "-") {
    type <- "both_negative"
  } else {
    type <- "outward_facing"
  }
  return(type)
}

is_read_of_interest <- function(align) {
  # Check if read fails any basic QC checks, and see if its a secondary or supplementary alignment:
  if (failedQC(align) || pcrORopt_duplicate(align) || secondaryAlign(align) || suppAlign(align)) {
    return(FALSE)
  }
  
  # Check that quality is good
  if (mapQuality(align) < 20) {
    return(FALSE)
  } 
  
  # Check if read is paired and both it and its mate are mapped
  if (!paired(align) || mateUnmapped(align) || unmapped(align)) {
    return(FALSE)
  }
  
  # If it's gotten this far it's a read of interest!
  return(TRUE)
}


# Main Plotting Loop ------------------------------------------------------

# Initialize bamReaders for every sample's bam file and store them in an environment/hash table
flog.info("Initializing bam file readers for all samples...")
bam_reader_env <- new.env()
for (i in 1:nrow(samples_and_fns.df)) {
  #progress(i, nrow(samples_and_fns.df))
  curr_sample <- as.character(samples_and_fns.df$sample_id[i])
  curr_bam.fn <- as.character(samples_and_fns.df$filename[i])
  bam_reader_env[[curr_sample]] <- bamReader(curr_bam.fn, indexname = paste(curr_bam.fn, ".bai", sep=""))
}

flog.info("Starting to loop through SVs...")
for (i in 1:nrow(vars.df)) {
  #i <- 1
  
  # Get initial stats about the SV
  # 1-based genome coordinates of the SV
  sv_type <- as.character(vars.df[i,]$sv_type)
  sv_ncbi_acc <- as.character(vars.df[i,]$ncbi_acc)
  sv_start <- vars.df[i,]$start_pos
  sv_end <- vars.df[i,]$end_pos
  sv_len <- vars.df[i,]$sv_len
  
  flog.info("Processing Structural Variant # %s: (%s) %s:%s-%s", i, sv_type, sv_ncbi_acc, sv_start, sv_end)

  # There are potentially two seperate intervals, or one if they overlap.
  # IRanges can merge two intervals easily
  regions_of_interest.ir <- IRanges(c(sv_start-1-args_flank_size, sv_end-1-args_flank_size),
                                    c(sv_start-1+args_flank_size, sv_end-1+args_flank_size))
  merged_regions_of_interest.ir <- reduce(regions_of_interest.ir)
  
  # Create an R dataframe for easy lookup of genotypes from a sample_id
  flog.trace("\tSetting up sample_id -> genotype lookup dataframe")
  genotypes.v <- unlist(vars.df[i, 10:ncol(vars.df)])
  genotype_lookup.df <- data.frame(sample_id=samples_and_fns.df$sample_id, genotype=genotypes.v)

  # The counts of each readpair type (occurring coincident with a 0/1 or 1/1 SV genotype) will be stored in readtype_counts_for_inds.df
  # This is using the setNames function on a matrix so that the names are hard-coded in when rows are added
  readpair_types <- c("both_negative", "both_positive", "outward_facing", "long", "normal")
  readtype_counts_for_inds.df <- setNames(data.frame(matrix(ncol=5, nrow=0)), readpair_types)
  
  # Initial setup for this SV is finished; now begin looping through all individuals.
  flog.info("\tProcessing individuals.")
  for (j in 1:nrow(samples_and_fns.df)) {
    #j <- 1
  
    if (j == nrow(samples_and_fns.df)) {
      flog.info("\t\t...%s...done", j)
    } else if (j %% 20 == 0) {
      flog.info("\t\t...%s", j)
    }
    
    #if (flog.threshold() %in% c("INFO")) {
    #  progress(j, nrow(samples_and_fns.df))
    #}

    # Basic info about the current (jth) sample
    curr_sample <- as.character(samples_and_fns.df$sample_id[j])
    curr_bam.fn <- as.character(samples_and_fns.df$filename[j])
    curr_bam_fn_base <- as.character(samples_and_fns.df$basename[j])
    
    flog.trace("\t\tIteration: %s: %s", j, curr_sample)
    flog.trace("\t\t\tBAM File: %s", curr_bam.fn)

    # Initialize a bamReader to process the bam file for this individual
    bam_reader <- bam_reader_env[[curr_sample]]
    bam_ref_data <- getRefData(bam_reader)
    
    max_normal_insert_size <- read_len_distro_info.df[read_len_distro_info.df$bam_file==curr_bam_fn_base,]$perc_99.5
    flog.trace("\t\t\tThe maximum insert size to be considered \"normal\" is: %s", max_normal_insert_size)
    
    # Reads of interest will be added to readpairs.l, forming lists indexed by their read ID.  Ultimately,
    #  only IDs associated with a list of two reads (the two ends) will be further processed.
    readpairs.l <- list()
    curr_contig_refid <- bam_ref_data[bam_ref_data$SN==sv_ncbi_acc,]$ID
    
    # First, get good reads from the one range that we know exists (the first, or only range)
    flog.trace("\t\t\tBeginning to process reads in the first bam range, from %s to %s", start(first_range.ir), end(first_range.ir))
    first_range.ir <- merged_regions_of_interest.ir[1]
    first_bam_range <- bamRange(bam_reader, c(curr_contig_refid, start(first_range.ir), end(first_range.ir)))
    rewind(first_bam_range) # Neurotic
    align <- getNextAlign(first_bam_range)
    while(!is.null(align)) {
      # Check if this is a read of interest
      # Reads of interest pass basic QC checks, are not secondary or supplemental, have a MAPQ >= 20, and are mapped
      #  along with their mates
      if (is_read_of_interest(align)) {
        if (is.null(readpairs.l[[name(align)]])) {
          readpairs.l[[name(align)]] <- list(align)
        } else {
          readpairs.l[[name(align)]] <- append(readpairs.l[[name(align)]][[1]], align)
        }
      }
      # Get the next read/alignment
      align <- getNextAlign(first_bam_range)
    }
    
    # Next, check if there is a second interval/range that was not merged with the first, and process if if necessary
    if (length(merged_regions_of_interest.ir)==2) {
      flog.trace("\t\t\tThere is a second unmerged bam range from %s to %s! Starting to process.", start(first_range.ir), end(first_range.ir))
      second_range.ir <- merged_regions_of_interest.ir[2]
      second_bam_range <- bamRange(bam_reader, c(curr_contig_refid, start(second_range.ir), end(second_range.ir)))
      rewind(second_bam_range) # Neurotic
      align <- getNextAlign(second_bam_range)
      while(!is.null(align)) {
        # Check if this is a read of interest
        # Reads of interest pass basic QC checks, are not secondary or supplemental, have a MAPQ >= 20, and are mapped
        #  along with their mates
        if (is_read_of_interest(align)) {
          if (is.null(readpairs.l[[name(align)]])) {
            readpairs.l[[name(align)]] <- list(align)
          } else {
            readpairs.l[[name(align)]] <- append(readpairs.l[[name(align)]][[1]], align)
          }
        }
        # Get the next read/alignment
        align <- getNextAlign(second_bam_range)
      }
    }
    flog.debug("\t\t\t%s distinct read IDs were found which had at least one read of interest", length(readpairs.l))
    
    # Having processed all relevent reads, we only want to retain those where both reads in a pair were discovered
    good_pairs.l <- readpairs.l[lapply(readpairs.l, length) == 2]
    flog.debug("\t\t\t%s read pairs were found where both reads were of interest", length(good_pairs.l))
    
    # Get the readpair types of each good read pair
    good_readpair_types <- unlist(lapply(good_pairs.l, get_readpair_type, max_normal_insert_size=max_normal_insert_size))
    
    # Count the number of each readpair type, and only append it to the others if the genotype for the current individual
    #  indicates the presence of at least one SV (non-reference) allele
    if (genotype_lookup.df$genotype[j] %in% c("0/1", "1/0", "1/1")) {
      readtype_counts <- c(both_negative=sum(good_readpair_types=="both_negative"),
                           both_negative=sum(good_readpair_types=="both_positive"),
                           both_negative=sum(good_readpair_types=="outward_facing"),
                           both_negative=sum(good_readpair_types=="long"),
                           both_negative=sum(good_readpair_types=="normal"))
      readtype_counts_for_inds.df[nrow(readtype_counts_for_inds.df)+1,] <- readtype_counts
    }
  }

  flog.info("\tDone processing processing individuals")
  
  # Get numbers of homR, homA and het genotypes
  num_homR_gt <- sum(genotype_lookup.df$genotype == "0/0")
  num_homA_gt <- sum(genotype_lookup.df$genotype == "1/1")
  num_het_gt <- sum(genotype_lookup.df$genotype %in% c("0/1", "1/0"))
  
  # Before we finish this loop, write the readpair information for this SV to file
  readtype_counts_summary <- apply(readtype_counts_for_inds.df, 2, mean)
  readtype_counts_summary_line <- c(c(i, sv_type, sv_len, sv_ncbi_acc, sv_start, sv_end, num_homR_gt, num_homA_gt, num_het_gt, 
                                      nrow(readtype_counts_for_inds.df)), readtype_counts_summary)

  # Finally, reset the output file if necessary, and write/append this SV's results to it
  if (i==1) {
    flog.trace("\tThis is the first SV in this analysis; starting a new output file (this will overwrite an existing file of the same name)")
    header_line <- c("sv_num", "sv_type", "sv_len", "sv_ncbi_acc", "sv_start", "sv_end", "num_homR_gt", "num_homA_gt", "num_het_gt",
                     "n_inds_with_at_least_one_sv_allele", "both_negative", "both_positive", "outward_facing", "long", "normal")
    write(header_line, args_readpair_info_fn, append=FALSE, sep='\t', ncolumns=length(readtype_counts_summary_line))
    append_val <- FALSE
  }
  flog.info("\tAppending results to file")
  write(readtype_counts_summary_line, args_readpair_info_fn, append=TRUE, sep='\t', ncolumns=length(readtype_counts_summary_line))
}

flog.info("Finished processing %s SVs", nrow(vars.df))
flog.info("Done! :-)")
