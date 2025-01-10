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
    output_prefix              = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'),

    # input files
    count_file                 = '$counts',
    sample_file                = '$samplesheet',

    # features and sample columns
    features_id_col            = 'gene_id',                  # column name of feature ids
    obs_id_col                 = 'sample',                   # column name of observation ids

    # comparison groups
    contrast_variable          = "$contrast_variable",       # column name of contrast variable
    reference_group            = "$reference",               # reference group for contrast variable
    target_group               = "$target",                  # target group for contrast variable
    blocking_variables         = NULL                        # column name of blocking variables

    # exclude samples
    subset_to_contrast_samples = FALSE,                      # keep only the samples involved in the contrast
    exclude_samples_col        = NULL,                       # exclude samples
    exclude_samples_values     = NULL,                       # exclude samples with given values

    # other parameters
    round_digits               = NA,                         # number of digits to round results
    ncores                     = as.integer('$task.cpus')    # number of cpus
)

opt_types <- list(
    output_prefix              = 'character',
    count_file                 = 'character',
    sample_file                = 'character',
    features_id_col            = 'character',
    obs_id_col                 = 'character',
    contrast_variable          = 'character',
    reference_group            = 'character',
    target_group               = 'character',
    blocking_variables         = 'character',
    subset_to_contrast_samples = 'character',
    exclude_samples_col        = 'character',
    exclude_samples_values     = 'character',
    round_digits               = 'numeric',
    ncores                     = 'numeric'
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

required_opts <- c('count_file','sample_file','contrast_variable','reference_group','target_group', 'prefix')
missing <- required_opts[unlist(lapply(opt[required_opts], is.null)) | ! required_opts %in% names(opt)]
if (length(missing) > 0){
    stop(paste("Missing required options:", paste(missing, collapse=', ')))
}

# Check file inputs are valid

for (file_input in c('count_file','sample_file')){
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
## READ IN COUNTS FILE AND SAMPLE METADATA    ##
################################################
################################################

count.table <-
    read_delim_flexible(
        file = opt\$count_file,
        header = TRUE,
        row.names = opt\$gene_id_col,
        check.names = FALSE
    )
sample.sheet <- read_delim_flexible(file = opt\$sample_file)

# Deal with spaces that may be in sample column
opt\$sample_id_col <- make.names(opt\$sample_id_col)

if (! opt\$sample_id_col %in% colnames(sample.sheet)){
    stop(paste0("Specified sample ID column '", opt\$sample_id_col, "' is not in the sample sheet"))
}

# Sample sheet can have duplicate rows for multiple sequencing runs, so uniqify
# before assigning row names

sample.sheet <- sample.sheet[! duplicated(sample.sheet[[opt\$sample_id_col]]), ]
rownames(sample.sheet) <- sample.sheet[[opt\$sample_id_col]]

# Check that all samples specified in the input sheet are present in the counts
# table. Assuming they are, subset and sort the count table to match the sample
# sheet

missing_samples <-
    sample.sheet[!rownames(sample.sheet) %in% colnames(count.table), opt\$sample_id_col]

if (length(missing_samples) > 0) {
    stop(paste(
        length(missing_samples),
        'specified samples missing from count table:',
        paste(missing_samples, collapse = ',')
    ))
} else{
    # Save any non-count data, will gene metadata etc we might need later
    noncount.table <-
        count.table[, !colnames(count.table) %in% rownames(sample.sheet), drop = FALSE]
    count.table <- count.table[, rownames(sample.sheet)]
}

################################################
################################################
## CHECK CONTRAST SPECIFICATION               ##
################################################
################################################

contrast_variable <- make.names(opt\$contrast_variable)
blocking.vars <- c()

if (!contrast_variable %in% colnames(sample.sheet)) {
    stop(
        paste0(
        'Chosen contrast variable \"',
        contrast_variable,
        '\" not in sample sheet'
        )
    )
} else if (any(!c(opt\$reference_level, opt\$target_level) %in% sample.sheet[[contrast_variable]])) {
    stop(
        paste(
        'Please choose reference and target levels that are present in the',
        contrast_variable,
        'column of the sample sheet'
        )
    )
} else if (!is.null(opt\$blocking_variables)) {
    blocking.vars = make.names(unlist(strsplit(opt\$blocking_variables, split = ';')))
    if (!all(blocking.vars %in% colnames(sample.sheet))) {
        missing_block <- paste(blocking.vars[! blocking.vars %in% colnames(sample.sheet)], collapse = ',')
        stop(
            paste(
                'Blocking variables', missing_block,
                'do not correspond to sample sheet columns.'
            )
        )
    }
}

# Handle conflicts between blocking variables and block

if (!is.null(opt\$block) && !is.null(opt\$blocking_variables)) {
    if (opt\$block %in% blocking.vars) {
        warning(paste("Variable", opt\$block, "is specified both as a fixed effect and a random effect. It will be treated as a random effect only."))
        blocking.vars <- setdiff(blocking.vars, opt\$block)
        if (length(blocking.vars) == 0) {
            opt\$blocking_variables <- NULL
        } else {
            opt\$blocking_variables <- paste(blocking.vars, collapse = ';')
        }
    }
}

# Optionally, subset to only the samples involved in the contrast

if (opt\$subset_to_contrast_samples) {
    sample_selector <- sample.sheet[[contrast_variable]] %in% c(opt\$target_level, opt\$reference_level)
    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

# Optionally, remove samples with specified values in a given field (probably
# don't use this as well as the above)

if ((! is.null(opt\$exclude_samples_col)) && (! is.null(opt\$exclude_samples_values))) {
    exclude_values = unlist(strsplit(opt\$exclude_samples_values, split = ';'))

    if (! opt\$exclude_samples_col %in% colnames(sample.sheet)){
        stop(paste(opt\$exclude_samples_col, ' specified to subset samples is not a valid sample sheet column'))
    }

    print(paste0('Excluding samples with values of ', opt\$exclude_samples_values, ' in ', opt\$exclude_samples_col))
    sample_selector <- ! sample.sheet[[opt\$exclude_samples_col]] %in% exclude_values

    selected_samples <- sample.sheet[sample_selector, opt\$sample_id_col]
    count.table <- count.table[, selected_samples]
    sample.sheet <- sample.sheet[selected_samples, ]
}

################################################
################################################
## Build the Model Formula                    ##
################################################
################################################

# Build the model formula with blocking variables first
model_vars <- c()

if (!is.null(opt\$blocking_variables)) {
    # Include blocking variables (including pairing variables if any)
    model_vars <- c(model_vars, blocking.vars)
}

# Add the contrast variable at the end
model_vars <- c(model_vars, contrast_variable)

# Construct the model formula
model <- paste('~ 0 +', paste(model_vars, collapse = '+'))

# Make sure all the appropriate variables are factors
vars_to_factor <- model_vars  # All variables in the model need to be factors
for (v in vars_to_factor) {
    sample.sheet[[v]] <- as.factor(sample.sheet[[v]])
}

################################################
################################################
## Run edgeR processes                        ##
################################################
################################################

# Create a DGEList object from counts

dge <- DGEList(counts = count.table)

# Normalize counts using TMM

dge <- calcNormFactors(dge, method = "TMM")


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
