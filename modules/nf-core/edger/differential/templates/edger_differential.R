#!/usr/bin/env Rscript

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and contains more than just whitespace.
#' It returns TRUE if the input is a non-empty, non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid, non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE

is_valid_string <- function(input) {
    !is.null(input) && nzchar(trimws(input))
}

#' Parse out options from a string without recourse to optparse
#'
#' @param x Long-form argument list like --opt1 val1 --opt2 val2
#'
#' @return named list of options and values similar to optparse

parse_args <- function(x){
    args_list <- unlist(strsplit(x, ' ?--')[[1]])[-1]
    args_vals <- lapply(args_list, function(x) scan(text=x, what='character', quiet = TRUE))

    # Ensure the option vectors are length 2 (key/ value) to catch empty ones
    args_vals <- lapply(args_vals, function(z){ length(z) <- 2; z})

    parsed_args <- structure(lapply(args_vals, function(x) x[2]), names = lapply(args_vals, function(x) x[1]))
    parsed_args[! is.na(parsed_args)]
}

#' Flexibly read CSV or TSV files
#'
#' @param file Input file
#' @param header Boolean. TRUE if first row is header. False without header.
#' @param row.names The first column is used as row names by default.
#' Otherwise, give another number. Or use NULL when no row.names are present.
#'
#' @return output Data frame

read_delim_flexible <- function(file, header = TRUE, row.names = 1, check.names = TRUE){

    ext <- tolower(tail(strsplit(basename(file), split = "\\\\.")[[1]], 1))

    if (ext == "tsv" || ext == "txt") {
        separator <- "\\t"
    } else if (ext == "csv") {
        separator <- ","
    } else {
        stop(paste("Unknown separator for", ext))
    }

    mat <- read.delim(
        file,
        sep = separator,
        header = header,
        row.names = row.names,
        check.names = check.names
    )
}

################################################
################################################
## Parse arguments                            ##
################################################
################################################

# Set defaults and classes

opt <- list(
    prefix             = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),

    # input count matrix
    counts             = '$counts',
    features_id_col    = 'gene_id',            # column name of feature ids

    # comparison groups
    samplesheet        = '$samplesheet',
    obs_id_col         = 'sample',             # column name of observation ids
    contrast_variable  = "$contrast_variable", # column name of contrast variable
    reference_group    = "$reference",         # reference group for contrast variable
    target_group       = "$target",            # target group for contrast variable
    blocking_variables = NULL                  # column name of blocking variables

    # other parameters
    round_digits       = NA,                   # number of digits to round results
    ncores             = as.integer('$task.cpus')
)

opt_types <- list(
    prefix             = 'character',
    counts             = 'character',
    features_id_col    = 'character',
    samplesheet        = 'character',
    obs_id_col         = 'character',
    contrast_variable  = 'character',
    reference_group    = 'character',
    target_group       = 'character',
    blocking_variables = 'character',
    round_digits       = 'numeric',
    ncores             = 'numeric'
)

# Apply parameter overrides

args_ext <- ifelse('$task.ext.args' == 'null', '', '$task.ext.args')
args_opt <- parse_args(args_ext)
for ( ao in names(args_opt)){
    if (! ao %in% names(opt)){
        stop(paste("Invalid option:", ao))
    } else {

        # Preserve classes from defaults where possible
        args_opt[[ao]] <- as(args_opt[[ao]], opt_types[[ao]])

        # handle NA, and avoid errors when NA is provided by user as character
        if (args_opt[[ao]] %in% c('NA', NA)) args_opt[[ao]] <- NA

        # replace values
        opt[[ao]] <- args_opt[[ao]]
    }
}

# Check if required parameters have been provided

required_opts <- c('counts','samplesheet','contrast_variable','reference_group','target_group', 'prefix')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('counts','samplesheet')){
    if (! is.is_valid_string(opt[[file_input]])) {
        stop(paste("Please provide", file_input), call. = FALSE)
    }
    if (! file.exists(opt[[file_input]])){
        stop(paste0('Value of ', file_input, ': ', opt[[file_input]], ' is not a valid file'))
    }
}

################################################
################################################
## Finish loading libraries                   ##
################################################
################################################

library(edgeR)

################################################
################################################
## Load data                                  ##
################################################
################################################

################################################
################################################
## Perform differential expression analysis   ##
################################################
################################################

################################################
################################################
## Generate outputs                           ##
################################################
################################################

################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink(paste0(opt\$prefix, ".R_sessionInfo.log"))
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

edger.version <- as.character(packageVersion('edgeR'))

writeLines(
    c(
        '"${task.process}":',
        paste('    bioconductor-edger:', edger.version)
    ),
'versions.yml')

################################################
################################################
################################################
################################################
