# server_tab1.R
# QC & Pre-processing Tool
# Tab 1: Data Input & Validation
# Author: Eren Ada, PhD
# Date: Created using the current date from system

library(shiny)
library(DT)
library(ggplot2)
library(plotly)

#' Server Module for Data Input & Validation
#' @param id Module ID
#' @param parent_session Parent session object
#' @importFrom shiny moduleServer reactive reactiveVal observeEvent req NS icon showNotification
#' @importFrom DT renderDT datatable
#' @importFrom utils read.csv read.delim write.csv
#' @importFrom dplyr %>% mutate select group_by summarise
mod_input_validation_server <- function(id, parent_session) {
  shiny::moduleServer(id, function(input, output, session) {
    
    # Initialize reactive values
    rv <- reactiveValues(
      raw_counts_data = NULL,            # Original count matrix
      raw_metadata = NULL,               # Original metadata
      detected_counts_separator = NULL,   # Detected separator for counts
      detected_metadata_separator = NULL, # Detected separator for metadata
      has_duplicates = FALSE,            # Whether duplicates were found
      duplicate_genes = NULL,            # List of duplicate genes
      files_analyzed = FALSE,            # Whether initial analysis is done
      data_processed = FALSE,            # Whether final processing is done
      validation_messages = list(),      # List of validation messages
      non_integer_warning = FALSE,       # Whether non-integer values were found
      processed_counts = NULL,           # Processed count matrix
      processed_metadata = NULL,         # Processed metadata
      original_gene_ids = NULL,          # Original gene IDs before any processing
      temp_renamed_ids = NULL,           # Temporarily renamed gene IDs
      original_counts_data = NULL        # Original data before filtering
    )
    
    # Helper function: Detect file separator
    detect_separator <- function(file_path) {
      # Read first few lines
      lines <- readLines(file_path, n = 5)
      
      # Count occurrences of potential separators
      comma_count <- max(sapply(lines, function(x) sum(gregexpr(",", x)[[1]] > 0)))
      tab_count <- max(sapply(lines, function(x) sum(gregexpr("\t", x)[[1]] > 0)))
      
      if (comma_count > 0 && comma_count >= tab_count) {
        return(list(sep = ",", type = "CSV"))
      } else if (tab_count > 0) {
        return(list(sep = "\t", type = "TSV"))
      } else {
        return(list(sep = ",", type = "CSV")) # Default to CSV
      }
    }
    
    # Helper function: Check for non-integer values
    check_non_integers <- function(data) {
      numeric_cols <- sapply(data, is.numeric)
      if (any(numeric_cols)) {
        non_integer_check <- sapply(data[, numeric_cols, drop = FALSE], function(x) {
          any(abs(x - round(x)) > .Machine$double.eps^0.5)
        })
        return(any(non_integer_check))
      }
      return(FALSE)
    }
    
    # Helper function: Round non-integer values
    round_non_integers <- function(data) {
      numeric_cols <- sapply(data, is.numeric)
      if (any(numeric_cols)) {
        data[, numeric_cols] <- round(data[, numeric_cols])
      }
      return(data)
    }
    
    # Helper function: Handle duplicates
    handle_duplicates <- function(data, strategy, gene_ids, remove_zeros = FALSE) {
      dupes <- gene_ids[duplicated(gene_ids)]
      if(length(dupes) == 0) {
        return(list(success = TRUE, data = data, message = "No duplicates found"))
      }
      
      result <- switch(strategy,
        "rename" = {
          new_names <- make.unique(gene_ids, sep = "_")
          rownames(data) <- new_names
          list(success = TRUE, 
               data = data, 
               message = "Duplicates renamed with suffix",
               details = data.frame(
                 Original = gene_ids,
                 New = new_names,
                 stringsAsFactors = FALSE
               ))
        },
        "sum" = {
          # Ensure data and gene_ids have the same length
          if(length(gene_ids) != nrow(data)) {
            return(list(
              success = FALSE,
              data = NULL,
              message = sprintf("Dimension mismatch: %d gene IDs but %d data rows", 
                              length(gene_ids), nrow(data))
            ))
          }
          
          # Convert data to matrix if it isn't already
          data_matrix <- as.matrix(data)
          
          # Split data by gene IDs and apply sum
          unique_genes <- unique(gene_ids)
          result_matrix <- t(sapply(unique_genes, function(g) {
            rows <- which(gene_ids == g)
            colSums(data_matrix[rows, , drop = FALSE], na.rm = TRUE)
          }))
          
          # Convert back to data frame
          agg_data <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
          rownames(agg_data) <- unique_genes
          colnames(agg_data) <- colnames(data)
          
          list(success = TRUE, 
               data = agg_data, 
               message = "Counts summed for duplicates",
               details = data.frame(
                 Gene = unique(dupes),
                 Count = table(gene_ids[gene_ids %in% dupes]),
                 stringsAsFactors = FALSE
               ))
        },
        "keep_highest_avg" = {
          avg_expr <- rowMeans(as.matrix(data))
          keep_idx <- tapply(1:length(gene_ids), gene_ids, function(idx) {
            idx[which.max(avg_expr[idx])]
          })
          data_filtered <- data[unlist(keep_idx), , drop = FALSE]
          list(success = TRUE, 
               data = data_filtered, 
               message = "Kept highest average expression for duplicates",
               details = data.frame(
                 Gene = names(keep_idx),
                 AvgExpr = avg_expr[unlist(keep_idx)],
                 stringsAsFactors = FALSE
               ))
        },
        "remove_all" = {
          keep <- !gene_ids %in% dupes
          data_filtered <- data[keep, , drop = FALSE]
          list(success = TRUE, 
               data = data_filtered, 
               message = "Removed all duplicate entries",
               details = data.frame(
                 RemovedGene = dupes,
                 stringsAsFactors = FALSE
               ))
        }
      )
      
      # Remove zero-count genes if requested
      if(remove_zeros && result$success) {
        nonzero_rows <- rowSums(result$data) > 0
        result$data <- result$data[nonzero_rows, , drop = FALSE]
        result$message <- paste(result$message, 
                              sprintf("and removed %d zero-count genes", 
                                    sum(!nonzero_rows)))
      }
      
      return(result)
    }
    
    # Observe Analyze Files button click
    observeEvent(input$analyze_files_btn, {
      req(input$count_matrix_file, input$metadata_file)
      
      # Reset states
      rv$files_analyzed <- FALSE
      rv$validation_messages <- list()
      rv$data_processed <- FALSE
      rv$processed_counts <- NULL
      rv$processed_metadata <- NULL
      
      # Detect separators
      counts_format <- detect_separator(input$count_matrix_file$datapath)
      metadata_format <- detect_separator(input$metadata_file$datapath)
      
      # Read files
      tryCatch({
        # First read without row.names to check for duplicates
        counts_data <- read.delim(input$count_matrix_file$datapath,
                                sep = counts_format$sep,
                                header = TRUE,
                                stringsAsFactors = FALSE)
        
        # Round numeric columns in counts_data immediately after reading
        numeric_cols <- sapply(counts_data, is.numeric)
        if(any(numeric_cols)) {
          has_decimals <- any(sapply(counts_data[, numeric_cols, drop = FALSE], function(x) {
            any(abs(x - round(x)) > .Machine$double.eps^0.5)
          }))
          
          if(has_decimals) {
            counts_data[, numeric_cols] <- round(counts_data[, numeric_cols], 0)
            rv$non_integer_warning <- TRUE
            showNotification(
              "Decimal values detected in count matrix. All values have been rounded to integers.",
              type = "warning",
              duration = 7
            )
          }
        }
        
        # Store original data before any filtering
        rv$original_gene_ids <- counts_data[,1]
        rv$original_counts_data <- counts_data
        
        # Check for zero-count genes immediately
        count_cols <- which(sapply(counts_data, is.numeric))
        zero_count_genes <- rowSums(counts_data[,count_cols, drop = FALSE]) == 0
        n_zero_genes <- sum(zero_count_genes)
        
        if(n_zero_genes > 0) {
          counts_data <- counts_data[!zero_count_genes,]
          gene_ids <- counts_data[,1]  # Update gene_ids after filtering
          
          showNotification(
            sprintf("Removed %d genes with zero counts across all samples (%.1f%% of total genes). You can download the full dataset using the button below.", 
                    n_zero_genes,
                    n_zero_genes/length(rv$original_gene_ids) * 100),
            type = "warning",
            duration = 7
          )
        } else {
          gene_ids <- counts_data[,1]
        }
        
        # Check for duplicates in filtered data
        duplicates <- gene_ids[duplicated(gene_ids)]
        rv$has_duplicates <- length(duplicates) > 0
        rv$duplicate_genes <- duplicates
        
        # Create temporary unique names for duplicates
        if(rv$has_duplicates) {
          temp_names <- make.unique(gene_ids, sep = "_TEMP_")
          rv$temp_renamed_ids <- temp_names
          
          # Now read the data with the temporary unique names
          counts_data[,1] <- temp_names
          rv$raw_counts_data <- counts_data[,-1]
          rownames(rv$raw_counts_data) <- temp_names
          
          # Calculate total expression for duplicates correctly
          duplicate_preview <- data.frame()
          for(gene in unique(duplicates)) {
            # Find all rows for this gene in original data
            gene_indices <- which(gene_ids == gene)
            
            # For each occurrence of the gene
            for(i in seq_along(gene_indices)) {
              # Calculate total expression correctly
              row_data <- as.numeric(counts_data[gene_indices[i], count_cols])
              total_expr <- sum(row_data)
              
              duplicate_preview <- rbind(duplicate_preview, data.frame(
                Gene_ID = gene,
                Entry = sprintf("Entry %d", i),
                Total_Expression = round(total_expr, 0),
                stringsAsFactors = FALSE
              ))
            }
          }
          
          # Store the preview data for the duplicate genes table
          rv$duplicate_preview <- duplicate_preview
          
          showNotification(
            sprintf("Found %d duplicate gene IDs. Temporarily renamed them for preview. Please select how to handle them in the options below.", 
                   length(unique(duplicates))),
            type = "warning",
            duration = 7
          )
        } else {
          # If no duplicates, read normally
          rv$raw_counts_data <- counts_data[,-1]
          rownames(rv$raw_counts_data) <- gene_ids
        }
        
        # Add download handler for original data
        output$download_original_data <- downloadHandler(
          filename = function() {
            "original_count_matrix.csv"
          },
          content = function(file) {
            write.csv(rv$original_counts_data, file, row.names = FALSE)
          }
        )
        
        # Read metadata
        rv$raw_metadata <- read.delim(input$metadata_file$datapath,
                                    sep = metadata_format$sep,
                                    header = TRUE,
                                    stringsAsFactors = FALSE)
        
        # Store detected formats
        rv$detected_counts_separator <- counts_format
        rv$detected_metadata_separator <- metadata_format
        
        # --- START: New Validation Checks ---
        
        # 1. Sample Name Consistency Check
        if (!is.null(rv$raw_counts_data) && !is.null(rv$raw_metadata)) {
          count_sample_names <- colnames(rv$raw_counts_data)
          metadata_sample_ids_for_check <- character(0)
          metadata_source_msg <- "metadata"

          if (ncol(rv$raw_metadata) > 0) {
            # Use the first column of metadata for sample ID comparison by default
            # This can be refined later when user can specify the sample ID column
            metadata_sample_ids_for_check <- as.character(rv$raw_metadata[[1]])
            metadata_source_msg <- "first column of metadata"
          } else {
            metadata_source_msg <- "metadata (no columns found for comparison)"
          }

          n_count_samples <- length(count_sample_names)
          n_metadata_ids <- length(unique(metadata_sample_ids_for_check))
          n_matching <- sum(count_sample_names %in% metadata_sample_ids_for_check)
          
          consistency_msg <- sprintf(
            "Counts vs %s: %d samples in counts, %d unique IDs in %s. %d matching.",
            metadata_source_msg,
            n_count_samples,
            n_metadata_ids,
            metadata_source_msg,
            n_matching
          )

          if (n_matching < n_count_samples || n_matching < n_metadata_ids) {
            counts_not_in_meta <- count_sample_names[!count_sample_names %in% metadata_sample_ids_for_check]
            meta_ids_not_in_counts <- unique(metadata_sample_ids_for_check[!metadata_sample_ids_for_check %in% count_sample_names])
            
            mismatch_details_list <- c()
            if (length(counts_not_in_meta) > 0) {
              if (length(counts_not_in_meta) <= 5) {
                detail <- sprintf("Samples in counts but not in %s: %s", 
                                metadata_source_msg, 
                                paste(counts_not_in_meta, collapse=", "))
              } else {
                detail <- sprintf("Samples in counts but not in %s: %s (and %d more)", 
                                metadata_source_msg, 
                                paste(head(counts_not_in_meta, 5), collapse=", "),
                                length(counts_not_in_meta) - 5)
              }
              mismatch_details_list <- c(mismatch_details_list, detail)
            }
            if (length(meta_ids_not_in_counts) > 0) {
              if (length(meta_ids_not_in_counts) <= 5) {
                detail <- sprintf("IDs in %s but not in counts: %s", 
                                metadata_source_msg, 
                                paste(meta_ids_not_in_counts, collapse=", "))
              } else {
                detail <- sprintf("IDs in %s but not in counts: %s (and %d more)", 
                                metadata_source_msg, 
                                paste(head(meta_ids_not_in_counts, 5), collapse=", "),
                                length(meta_ids_not_in_counts) - 5)
              }
              mismatch_details_list <- c(mismatch_details_list, detail)
            }
            if (length(mismatch_details_list) > 0) {
              consistency_msg <- paste(consistency_msg, "\nDetails:", paste(mismatch_details_list, collapse="\n"))
            }
            if (n_matching == 0 && n_count_samples > 0 && n_metadata_ids > 0) {
              consistency_msg <- paste(consistency_msg, "\nWARNING: No overlap in sample identifiers found!")
            }
          } else if (n_count_samples > 0) {
            consistency_msg <- paste(consistency_msg, "All sample names appear consistent.")
          } else {
            consistency_msg <- paste(consistency_msg, "Not enough data for full consistency check.")
          }
          rv$validation_messages$sample_consistency <- consistency_msg
        }

        # 2. Missing Value (NA) Summary
        if (!is.null(rv$raw_counts_data) && !is.null(rv$raw_metadata)) {
          na_in_counts <- sum(is.na(rv$raw_counts_data))
          na_in_metadata <- sum(is.na(rv$raw_metadata))
          rv$validation_messages$na_summary <- sprintf(
            "%d NAs found in count data. %d NAs found in metadata.",
            na_in_counts,
            na_in_metadata
          )
        }
        # --- END: New Validation Checks ---
        
        # Set analysis complete flag
        rv$files_analyzed <- TRUE
        
      }, error = function(e) {
        showNotification(paste("Error reading files:", e$message), type = "error", duration = 10)
      })
    })
    
    # Render detected settings
    output$detected_settings <- renderUI({
      req(rv$files_analyzed)
      
      # Calculate number of zero-count genes for display
      zero_count_genes <- sum(rowSums(rv$raw_counts_data) == 0)
      
      tagList(
        tags$div(class = "detected-settings",
          tags$p(icon("check"), 
                sprintf("Count Matrix Format: %s (Separator: %s)",
                        rv$detected_counts_separator$type,
                        if(rv$detected_counts_separator$sep == "\t") "Tab" else "Comma")),
          tags$p(icon("check"),
                sprintf("Metadata Format: %s (Separator: %s)",
                        rv$detected_metadata_separator$type,
                        if(rv$detected_metadata_separator$sep == "\t") "Tab" else "Comma")),
          if(rv$non_integer_warning) {
            tags$p(icon("check"), 
                  "Non-integer values were detected and automatically rounded.",
                  class = "text-success")
          },
          # Add zero-count genes option
          if(zero_count_genes > 0) {
            tags$div(
              class = "alert alert-warning",
              tags$p(
                icon("exclamation-triangle"),
                sprintf("Found %s genes with zero counts across all samples.", 
                       format(zero_count_genes, big.mark = ",")),
                style = "margin-bottom: 10px;"
              ),
              checkboxInput(
                session$ns("remove_zero_counts"),
                label = sprintf("Remove zero-count genes (%s genes)", 
                              format(zero_count_genes, big.mark = ",")),
                value = FALSE
              )
            )
          }
        )
      )
    })
    
    # Render data issues
    output$data_issues <- renderUI({
      req(rv$files_analyzed)
      
      issues <- list()
      
      if(rv$has_duplicates) {
        duplicate_count <- length(unique(rv$duplicate_genes))
        total_duplicates <- length(rv$duplicate_genes)
        
        issues$duplicates <- tags$div(
          class = "alert alert-warning",
          tags$div(
            class = "warning-header",
            icon("exclamation-triangle"), 
            tags$strong(" Duplicate Genes Detected")
          ),
          tags$p(
            sprintf("Found %d duplicate gene IDs that need to be resolved before proceeding.", 
                   duplicate_count)
          ),
          tags$div(
            class = "duplicate-handling",
            tags$h4("Duplicate Gene Handling"),
            tags$p("Select handling method:"),
            selectInput(
              session$ns("duplicate_strategy"),
              label = NULL,
              choices = list(
                "Sum Counts for Duplicates" = "sum",
                "Keep Entry with Highest Average Expression" = "keep_highest_avg",
                "Remove All Duplicate Entries" = "remove_all",
                "Rename with Suffix" = "rename"
              ),
              selected = "sum"
            )
          )
        )
      }
      
      do.call(tagList, issues)
    })
    
    # Preview of duplicate handling results
    output$duplicate_preview <- renderUI({
      req(rv$has_duplicates, input$duplicate_strategy)
      
      preview_text <- switch(input$duplicate_strategy,
        "rename" = "Each duplicate gene will be renamed with a suffix (e.g., GENE_1, GENE_2)",
        "sum" = "Counts will be summed for each duplicate gene",
        "keep_highest_avg" = "Only the entry with highest average expression will be kept",
        "remove_all" = "All entries of duplicated genes will be removed"
      )
      
      tags$div(
        class = "alert alert-info",
        tags$div(
          class = "preview-header",
          icon("info-circle"), 
          tags$strong(" Preview of selected strategy:")
        ),
        tags$p(preview_text),
        tags$p(class = "text-muted", "This change will be applied when you click 'Process Data'")
      )
    })
    
    # Observe Process Data button click
    observeEvent(input$process_data_btn, {
      req(rv$files_analyzed)
      
      # Process count matrix
      counts_data <- rv$raw_counts_data
      metadata <- rv$raw_metadata
      
      # Add debugging messages
      message("=== Processing Data Debug ===")
      message(sprintf("Initial counts_data dimensions: %d rows x %d cols", nrow(counts_data), ncol(counts_data)))
      message(sprintf("Number of gene IDs: %d", length(rv$original_gene_ids)))
      
      # Ensure metadata has Sample column
      if (!"Sample" %in% colnames(metadata)) {
        # If first column contains sample IDs, rename it
        colnames(metadata)[1] <- "Sample"
      }
      
      # Ensure sample IDs match between counts and metadata
      sample_ids_counts <- colnames(counts_data)
      sample_ids_metadata <- metadata$Sample
      
      # Check if all count matrix samples are in metadata
      if (!all(sample_ids_counts %in% sample_ids_metadata)) {
        missing_samples <- sample_ids_counts[!sample_ids_counts %in% sample_ids_metadata]
        showNotification(
          sprintf("Warning: Some samples in count matrix are missing from metadata: %s", 
                  paste(head(missing_samples, 3), collapse=", ")),
          type = "warning",
          duration = NULL
        )
      }
      
      # Filter metadata to match count matrix samples and order
      metadata <- metadata[match(sample_ids_counts, sample_ids_metadata), , drop = FALSE]
      rownames(metadata) <- metadata$Sample
      
      # Handle duplicates if needed
      if(rv$has_duplicates) {
        # Get the current gene IDs after any filtering
        current_gene_ids <- if(!is.null(rv$temp_renamed_ids)) {
          rv$temp_renamed_ids
        } else {
          rv$original_gene_ids[!is.na(match(rv$original_gene_ids, rownames(counts_data)))]
        }
        
        message(sprintf("Current gene IDs length: %d", length(current_gene_ids)))
        message(sprintf("Counts data rows: %d", nrow(counts_data)))
        
        duplicate_result <- handle_duplicates(
          counts_data, 
          input$duplicate_strategy, 
          current_gene_ids,
          remove_zeros = !is.null(input$remove_zero_counts) && input$remove_zero_counts
        )
        
        if(!duplicate_result$success) {
          showNotification(duplicate_result$message, type = "error", duration = 10)
          return()
        }
        
        counts_data <- duplicate_result$data
        
        # Show summary of what was done
        showNotification(
          duplicate_result$message,
          type = "default",
          duration = 5
        )
      } else if(!is.null(input$remove_zero_counts) && input$remove_zero_counts) {
        # Handle zero-count removal even if there are no duplicates
        nonzero_rows <- rowSums(counts_data) > 0
        counts_data <- counts_data[nonzero_rows, , drop = FALSE]
        showNotification(
          sprintf("Removed %d zero-count genes", sum(!nonzero_rows)),
          type = "default",
          duration = NULL
        )
      }
      
      # Store processed data
      rv$processed_counts <- counts_data
      rv$processed_metadata <- metadata
      rv$data_processed <- TRUE
      
      # Show success message
      showNotification("Data processed successfully!", type = "default", duration = 5)
    })
    
    # Render Count Matrix Preview
    output$count_matrix_preview <- renderDT({
      req(rv$files_analyzed)
      # Show processed data if available, otherwise show raw data
      preview_data <- if(rv$data_processed && !is.null(rv$processed_counts)) {
        rv$processed_counts
      } else {
        rv$raw_counts_data
      }
      
      # Get numeric columns
      numeric_cols <- which(sapply(preview_data, is.numeric))
      
      dt <- datatable(preview_data, 
                options = list(
                  scrollX = TRUE,
                  pageLength = 10,
                  lengthMenu = list(c(10, 25, 50), c('10', '25', '50')),
                  dom = 'lrtip',
                  searchHighlight = TRUE,
                  search = list(regex = TRUE, caseInsensitive = TRUE)
                ),
                class = 'cell-border stripe')
      
      # Apply integer formatting to all numeric columns
      if(length(numeric_cols) > 0) {
        for(col in numeric_cols) {
          dt <- dt %>% formatRound(col, digits = 0)
        }
      }
      
      dt
    })

    # Render initial Data Summary
    output$data_summary <- renderUI({
      req(rv$files_analyzed)
      
      # Calculate total counts
      total_counts <- sum(rv$raw_counts_data, na.rm = TRUE)
      total_counts_formatted <- if(total_counts > 1e9) {
        sprintf("%.1fB", total_counts/1e9)
      } else if(total_counts > 1e6) {
        sprintf("%.1fM", total_counts/1e6)
      } else {
        format(total_counts, big.mark = ",")
      }
      
      # Calculate number of zero-count genes
      zero_count_genes <- sum(rowSums(rv$raw_counts_data) == 0)
      zero_count_percentage <- round(zero_count_genes/nrow(rv$raw_counts_data) * 100, 1)
      
      # Get metadata column names
      metadata_cols <- colnames(rv$raw_metadata)
      metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample", "sample", "Sample_ID", "sample_id")]
      
      # Validation checks - Use same logic as initial validation
      count_sample_names <- colnames(rv$raw_counts_data)
      metadata_sample_ids <- if (ncol(rv$raw_metadata) > 0) {
        as.character(rv$raw_metadata[[1]])
      } else {
        character(0)
      }
      all_samples_consistent <- length(count_sample_names) > 0 && 
                               length(metadata_sample_ids) > 0 && 
                               all(count_sample_names %in% metadata_sample_ids)
      no_missing_values <- sum(is.na(rv$raw_counts_data)) == 0 && sum(is.na(rv$raw_metadata)) == 0
      
      tagList(
        # Main Summary Box
        div(class = "data-summary-box",
          # Counts Section
          div(class = "summary-section",
            h5(icon("table"), "Count Matrix"),
            tags$ul(
              tags$li(
                strong(format(nrow(rv$raw_counts_data), big.mark = ",")), " genes",
                if(zero_count_genes > 0) {
                  tags$span(class = "text-warning",
                    sprintf(" (including %s zero-count genes)", 
                           format(zero_count_genes, big.mark = ","))
                  )
                }
              ),
              tags$li(strong(ncol(rv$raw_counts_data)), " samples"),
              tags$li(strong(total_counts_formatted), " total counts")
            )
          ),
          
          # Metadata Section
          div(class = "summary-section",
            h5(icon("info-circle"), "Metadata"),
            tags$ul(
              tags$li(
                strong(length(metadata_cols)), " attributes: ",
                tags$span(class = "text-muted", paste(metadata_cols, collapse = ", "))
              )
            )
          ),
          
          # Validation Section
          div(class = "summary-section",
            h5(icon("check-circle"), "Validation Status"),
            tags$ul(
              tags$li(
                if(all_samples_consistent) icon("check", class = "text-success") else icon("times", class = "text-danger"),
                " Sample matching: ",
                if(all_samples_consistent) "All samples matched" else "Some samples mismatched"
              ),
              tags$li(
                if(no_missing_values) icon("check", class = "text-success") else icon("times", class = "text-danger"),
                " Data completeness: ",
                if(no_missing_values) "No missing values" else "Contains missing values"
              )
            ),
            # Add detailed validation messages here
            uiOutput(session$ns("validation_details"))
          )
        )
      )
    })

    # Add new output for detailed validation messages
    output$validation_details <- renderUI({
      req(rv$files_analyzed)
      
      details_list <- list()
      
      # Sample consistency details - only show if there are actual mismatches
      if (!is.null(rv$validation_messages$sample_consistency)) {
        # Check if there are actual mismatches (not just the summary message)
        if (grepl("Details:", rv$validation_messages$sample_consistency)) {
          # Split the message by lines for better formatting
          msg_parts <- strsplit(rv$validation_messages$sample_consistency, "\n")[[1]]
          
          # Create formatted content
          formatted_content <- lapply(msg_parts, function(line) {
            if (grepl("^Details:", line)) {
              div(class = "details-text", line)
            } else if (grepl("WARNING:", line)) {
              div(class = "warning-text", line)
            } else {
              div(class = "sample-mismatch-line", line)
            }
          })
          
          details_list$sample_details <- div(
            class = "alert alert-warning validation-details",
            h6(icon("exclamation-triangle"), " Sample Matching Details"),
            div(
              style = "font-size: 0.9em;",
              formatted_content
            )
          )
        }
      }
      
      # Missing values details
      if (!is.null(rv$validation_messages$na_summary)) {
        na_counts <- sum(is.na(rv$raw_counts_data)) + sum(is.na(rv$raw_metadata))
        if (na_counts > 0) {
          details_list$na_details <- div(
            class = "alert alert-info validation-details",
            h6(icon("info-circle"), " Missing Values Summary"),
            p(rv$validation_messages$na_summary, style = "margin-bottom: 0; font-size: 0.9em;")
          )
        }
      }
      
      # Show the details only if there are issues to report
      if (length(details_list) > 0) {
        do.call(tagList, details_list)
      } else {
        NULL
      }
    })

    # Render library size plot
    output$library_size_plot <- renderPlotly({
      req(rv$files_analyzed)
      
      # Calculate library sizes (sum of counts for each sample)
      library_sizes <- colSums(rv$raw_counts_data)
      
      # Create data frame for plotting
      plot_data <- data.frame(
        Sample = names(library_sizes),
        Size = as.numeric(library_sizes)
      )
      
      # Create the plot
      p <- ggplot(plot_data, aes(x = reorder(Sample, -Size), y = Size, fill = Sample)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none"
        ) +
        labs(
          title = "Library Size per Sample",
          x = "Sample",
          y = "Total Counts"
        ) +
        scale_y_continuous(labels = scales::comma)
      
      # Convert to plotly with explicit height
      ggplotly(p, height = 400) %>%
        layout(margin = list(b = 100))
    })

    # Render Processed Data Summary (restored)
    output$processed_data_summary <- renderUI({
      req(rv$data_processed, rv$processed_counts)
      
      # Calculate total counts for processed data
      total_counts_processed <- sum(rv$processed_counts, na.rm = TRUE)
      total_counts_original <- sum(rv$raw_counts_data, na.rm = TRUE)
      
      # Calculate the changes
      original_genes <- nrow(rv$raw_counts_data)
      processed_genes <- nrow(rv$processed_counts)
      genes_retained_pct <- round(processed_genes / original_genes * 100, 1)
      
      # Calculate zero-count genes
      zero_count_genes <- sum(rowSums(rv$raw_counts_data) == 0)
      zero_removed <- if(!is.null(input$remove_zero_counts) && input$remove_zero_counts) zero_count_genes else 0
      
      # Calculate duplicate stats
      dup_count <- if(rv$has_duplicates) length(unique(rv$duplicate_genes)) else 0
      
      tagList(
        # 1. Processing Results (first)
        div(class = "summary-section",
          h5(icon("cogs"), "Processing Results"),
          tags$ul(
            tags$li(
              strong(format(processed_genes, big.mark = ",")), " genes retained ",
              tags$small(class = "text-muted",
                sprintf("(%.1f%%)", genes_retained_pct)
              )
            ),
            if(zero_removed > 0) {
              tags$li(
                icon("minus-circle"), 
                sprintf(" Removed %s zero-count genes",
                        format(zero_removed, big.mark = ","))
              )
            },
            if(rv$has_duplicates) {
              tags$li(
                icon("clone"),
                sprintf(" Handled %d duplicate genes",
                        dup_count)
              )
            }
          )
        ),
        
        tags$hr(),
        
        # 2. Processed Data Summary (second, without redundant info)
        h5("Processed Data Summary"),
        tags$ul(
          tags$li(
            sprintf("Number of samples: %s (unchanged)", 
                    format(ncol(rv$processed_counts), big.mark = ","))
          ),
          if(rv$has_duplicates) {
            strategy_text <- switch(input$duplicate_strategy,
                                   "rename" = "renamed with suffix",
                                   "sum" = "summed",
                                   "keep_highest_avg" = "kept only highest expression",
                                   "remove_all" = "removed")
            
            tags$li(
              sprintf("Duplicate handling: %s", strategy_text)
            )
          }
        ),
        
        tags$hr(),
        
        # 3. Detailed Statistics (third)
        h6("Detailed Statistics"),
        tags$ul(
          tags$li(sprintf("Original genes: %s", format(original_genes, big.mark = ","))),
          tags$li(sprintf("Final genes: %s", format(processed_genes, big.mark = ","))),
          tags$li(sprintf("Exact total counts: %s", format(total_counts_processed, big.mark = ","))),
          tags$li("Sample characteristics:", 
                 tags$ul(
                   lapply(colnames(rv$raw_metadata)[!colnames(rv$raw_metadata) %in% c("Sample", "sample", "Sample_ID", "sample_id")], function(col) {
                     tags$li(sprintf("%s: %d unique values", 
                                   col, 
                                   length(unique(rv$raw_metadata[[col]]))))
                   })
                 ))
        )
      )
    })

    # Add download handler for duplicate genes statistics
    output$download_duplicate_stats <- downloadHandler(
      filename = function() {
        paste("duplicate_genes_stats_", format(Sys.time(), "%Y%m%d"), ".csv", sep="")
      },
      content = function(file) {
        req(rv$has_duplicates)
        req(rv$original_counts_data)

        # Work off the original uploaded table to avoid index drift after filtering/renaming
        original_df <- rv$original_counts_data
        # Identify numeric count columns
        count_cols <- which(sapply(original_df, is.numeric))
        if (length(count_cols) == 0) {
          # Fallback: treat all but first column as numeric if types were coerced
          count_cols <- setdiff(seq_len(ncol(original_df)), 1L)
        }

        # List of duplicated gene IDs in original space
        dup_genes <- unique(original_df[[1]][original_df[[1]] %in% rv$duplicate_genes])

        # Prepare data for CSV
        result <- data.frame(stringsAsFactors = FALSE)

        for (gene in dup_genes) {
          # Get all rows for this gene in the original data
          gene_rows <- which(original_df[[1]] == gene)
          if (length(gene_rows) == 0) next
          # Sum counts per occurrence across numeric columns
          expressions <- apply(original_df[gene_rows, count_cols, drop = FALSE], 1, function(x) sum(as.numeric(x), na.rm = TRUE))

          # Append one row per occurrence
          for (i in seq_along(expressions)) {
            result <- rbind(result, data.frame(
              Gene_ID = gene,
              Entry = i,
              Total_Expression = round(expressions[i], 0),
              stringsAsFactors = FALSE
            ))
          }
        }

        write.csv(result, file, row.names = FALSE)
      }
    )

    # Add download handler for processed count matrix
    output$download_processed_matrix <- downloadHandler(
      filename = function() {
        paste("processed_count_matrix_", format(Sys.time(), "%Y%m%d"), ".csv", sep="")
      },
      content = function(file) {
        req(rv$data_processed, rv$processed_counts)
        # Add gene IDs as first column
        output_data <- data.frame(
          Gene_ID = rownames(rv$processed_counts),
          rv$processed_counts,
          check.names = FALSE
        )
        write.csv(output_data, file, row.names = FALSE)
      }
    )
    
    # Render Metadata Preview
    output$metadata_preview <- renderDT({
      req(rv$files_analyzed)
      datatable(rv$raw_metadata, 
                options = list(
                  scrollX = TRUE,
                  pageLength = 10,
                  lengthMenu = list(c(10, 25, 50), c('10', '25', '50')),
                  dom = 'lrtip',
                  searchHighlight = TRUE,
                  search = list(regex = TRUE, caseInsensitive = TRUE)
                ),
                class = 'cell-border stripe')
    })
    
    # Render duplicate genes table
    output$duplicate_genes_table <- renderDT({
      req(rv$has_duplicates, rv$duplicate_preview)
      
      # Use the pre-calculated preview data
      datatable(rv$duplicate_preview,
                options = list(
                  pageLength = 10,
                  lengthMenu = list(c(10, 25, 50, -1), c('10', '25', '50', 'All')),
                  scrollX = TRUE,
                  order = list(list(2, 'desc')),  # Sort by Total Expression by default
                  dom = 'lrtip'  # Remove search box
                ),
                rownames = FALSE,
                class = 'cell-border stripe')
    })
    
    # Output flags for UI conditionals
    output$files_analyzed <- reactive({ rv$files_analyzed })
    outputOptions(output, "files_analyzed", suspendWhenHidden = FALSE)
    
    output$has_duplicates <- reactive({ rv$has_duplicates })
    outputOptions(output, "has_duplicates", suspendWhenHidden = FALSE)
    
    output$data_processed <- reactive({ rv$data_processed })
    outputOptions(output, "data_processed", suspendWhenHidden = FALSE)
    
    # Handle Next Step button click
    observeEvent(input$next_step_btn, {
      req(rv$data_processed, rv$processed_counts)
      
      # Navigate to the QC Plots tab
      updateNavbarPage(parent_session, "mainNav", selected = "qc_plots")
      
      # Show notification that we're moving to QC plots
      showNotification(
        "Moving to QC Plots & Summaries tab with processed data",
        type = "default",
        duration = 3
      )
    })
    
    # Processing completion message for proceed section
    output$processing_completion_message <- renderUI({
      source("modules/module1_qc_preprocessing/R/shared_ui_components.R", local = TRUE)
      
      if (rv$data_processed) {
        create_status_message(
          "Data processing complete. You can continue to QC analysis.",
          status = "success"
        )
      } else if (rv$files_analyzed) {
        create_status_message(
          "Please process your data before proceeding to QC analysis.",
          status = "warning"
        )
      } else {
        create_status_message(
          "Please upload and validate your files first.",
          status = "info"
        )
      }
    })
    
    # Return reactive values for use in other modules
    return(list(
      processed_counts = reactive({ rv$processed_counts }),
      processed_metadata = reactive({ rv$processed_metadata }),
      data_processed = reactive({ rv$data_processed })
    ))
  })
} 