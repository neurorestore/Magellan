#' Perform spatial prioritization of the response to a biological perturbation
#'
#' Prioritize spatial locations involved in a complex biological process by 
#' training a machine-learning model to predict sample labels 
#' (e.g., disease vs. control, treated vs. untreated, or time post-stimulus), 
#' and evaluate the performance of the model in cross-validation.
#'
#' If a \code{Seurat} object is provided as input, Magellan will use the default
#' assay (i.e., whatever \link[Seurat]{GetAssayData} returns) as input. To
#' use a different assay, provide the expression matrix and metadata as input
#' separately, using the \code{input} and \code{meta} arguments. 
#' Additionally, Magellan will assume the coordinates of the spatial barcodes
#' can be found in \code{input@images$slice1@coordinates}. To override this, 
#' specify the count matrix, metadata, and coordinates separately. 
#'
#' @param input a matrix, data frame, or \code{Seurat} object containing 
#'   gene expression values (genes in rows, cells in columns) and, 
#'   optionally, metadata about each spatial barcode
#' @param meta optionally, a data frame containing metadata about the
#'   \code{input} gene-by-barcode matrix, at minimum containing the label 
#'   associated with each barcode
#'   (e.g., group, disease, timepoint); can be left as
#'   \code{NULL} if \code{input} is a \code{Seurat} object
#' @param coords optionally, a data frame containing the spatial coordinates 
#'   for each barcode in the \code{input} gene-by-barcode matrix; can be left as
#'   \code{NULL} if \code{input} is a \code{Seurat} object
#' @param k the number of spatial nearest-neighbors to use in the AUC 
#'   calculation; defaults to \code{50}
#' @param label_col the column of the \code{meta} data frame, or the
#'   metadata container in the \code{Seurat} object, that
#'   contains condition labels (e.g., disease, timepoint) for each barcode in the
#'   gene-by-barcode expression matrix; defaults to \code{"label"}
#' @param coord_cols the names of the columns in the \code{coords} data frame, 
#'   or the metadata container in the \code{Seurat} object, that contain the 
#'   coordinates of each spatial barcode in the gene-by-barcode expression 
#'   matrix; defaults to \code{c("coord_x", "coord_y")}
#' @param n_subsamples the number of times to repeat the cross-validation 
#'   procedure for each barcode; defaults to \code{50}.
#'   Set to \code{0} to omit subsampling altogether,
#'   calculating performance on the entire dataset, but note that this may
#'   introduce bias due to cell type or label class imbalance. 
#'   Note that when setting \code{augur_mode = "permute"}, values less than 
#'   \code{100} will be replaced with a default of \code{500}. 
#' @param subsample_size the number of barcodes to randomly sample from among 
#'   the nearest neighbors in each iteration of the cross-validation procedure;
#'   cannot be greater than \code{k}
#' @param folds the number of folds of cross-validation to run; defaults to
#'   \code{3}
#' @param var_quantile the quantile of highly variable genes to retain using 
#'   the variable gene filter (\link[Augur]{select_variance}); 
#'   defaults to \code{0.5}
#' @param feature_perc the proportion of genes that are randomly selected as
#'   features for input to the classifier in each subsample using the
#'   random gene filter (\link[Augur]{select_random}); defaults to \code{0.5}
#' @param n_threads the number of threads to use for parallelization;
#'   defaults to \code{32}.
#' @param show_progress if \code{TRUE}, display a progress bar for the analysis
#'   with estimated time remaining
#' @param augur_mode one of \code{"default"} or \code{"velocity"}. 
#'   Setting \code{augur_mode = "velocity"} disables feature
#'   selection, assuming feature selection has been performed by the RNA
#'   velocity procedure to produce the input matrix
#' @param classifier the classifier to use in calculating area under the curve,
#'   one of \code{"rf"} (random forest) or \code{"lr"} (logistic regression);
#'   defaults to \code{"rf"}, which is the recommended setting
#' @param rf_params for \code{classifier} == \code{"rf"}, a list of parameters
#'   for the random forest models, containing the following items (see
#'   \link[parsnip]{rand_forest} from the \code{parsnip} package):
#'   \describe{
#'     \item{"mtry"}{the number of features randomly sampled at each split
#'       in the random forest classifier; defaults to \code{2}}
#'     \item{"trees"}{the number of trees in the random forest classifier;
#'       defaults to \code{100}}
#'     \item{"min_n"}{the minimum number of observations to split a node in the
#'       random forest classifier; defaults to \code{NULL}}
#'     \item{"importance"}{the method of calculating feature importances
#'       to use; defaults to \code{"accuracy"}; can also specify \code{"gini"}}
#'   }
#' @param lr_params for \code{classifier} == \code{"lr"}, a list of parameters
#'   for the logistic regression models, containing the following items (see
#'   \link[parsnip]{logistic_reg} from the \code{parsnip} package):
#'   \describe{
#'     \item{"mixture"}{the proportion of L1 regularization in the model;
#'       defaults to \code{1}}
#'     \item{"penalty"}{the total amount of regularization in the model;
#'       defaults to \code{"auto"}, which uses \link[glmnet]{cv.glmnet} to set
#'       the penalty}
#'   }
#'
#' @return a list of class \code{"Magellan"}, containing the following items:
#' \enumerate{
#'   \item \code{parameters}: the parameters provided to this function as input
#'   \item \code{results}: the area under the curve for each barcode, in each
#'     fold, in each subsample, in the comparison of interest, as well as a
#'     series of other classification metrics
#'   \item \code{AUC}: a summary of the mean AUC for each barcode (for
#'     continuous experimental conditions, this is replaced by a \code{CCC}
#'     item that records the mean concordance correlation coefficient for each
#'     barcode)
#' }
#'
#' @importFrom dplyr do sample_n group_by ungroup tibble mutate select bind_rows
#'   pull rename n_distinct arrange desc filter summarise row_number left_join n
#'   all_of
#' @importFrom tibble repair_names rownames_to_column
#' @importFrom purrr map map2 map_lgl pmap map2_df
#' @importFrom magrittr %>% %<>% extract extract2 set_rownames set_colnames
#' @importFrom parsnip set_engine logistic_reg rand_forest fit translate
#' @importFrom rsample assessment analysis
#' @importFrom recipes prepper bake recipe
#' @importFrom yardstick metric_set accuracy precision recall sens spec npv
#'   ppv roc_auc ccc huber_loss_pseudo huber_loss mae mape mase
#'   rpd rpiq rsq_trad rsq smape rmse
#' @importFrom stats setNames predict sd
#' @importFrom methods is
#' @importFrom sparseMatrixStats colVars
#' @importFrom pbmcapply pbmclapply
#' @importFrom parallel mclapply
#' @importFrom tester is_numeric_matrix is_numeric_dataframe
#' @importFrom Augur select_variance calculate_auc
#' @importFrom RANN nn2
#' @importFrom tidyr gather
#' @import Matrix
#'
#' @export
navigate_space = function(input,
                          meta = NULL,
                          coords = NULL,
                          k = 50,
                          label_col = "label",
                          coord_cols = c("coord_x", "coord_y"),
                          n_subsamples = 50,
                          subsample_size = 20,
                          folds = 3,
                          var_quantile = 0.5,
                          feature_perc = 0.5,
                          n_threads = 32,
                          show_progress = T,
                          augur_mode = c('default', 'velocity'),
                          classifier = c("rf", "lr"),
                          # random forest parameters
                          rf_params = list(trees = 100,
                                           mtry = 2,
                                           min_n = NULL,
                                           importance = 'accuracy'),
                          # logistic regression parameters
                          lr_params = list(mixture = 1, penalty = 'auto')
) {
  # check arguments
  classifier = match.arg(classifier)
  augur_mode = match.arg(augur_mode)
  # check number of folds/subsample size are compatible (n > 1 in every fold)
  if (n_subsamples > 1 & subsample_size / folds < 2) {
    stop("subsample_size / n_folds must be greater than or equal to 2")
  }
  # subsample_size cannot be larger than k
  if (subsample_size > k) {
    stop('subsample_size cannot be greater than k')
  }
  
  # check glmnet is installed for logistic regression
  if (classifier == 'lr' && !requireNamespace("glmnet", quietly = TRUE)) {
    stop("install \"glmnet\" R package to run Augur with logistic regression ",
         "classifier", call. = FALSE)
  }
  
  # extract cell types and label from metadata
  if ("Seurat" %in% class(input)) {
    # confirm Seurat is installed
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("install \"Seurat\" R package for Augur compatibility with ",
           "input Seurat object", call. = FALSE)
    }
    
    meta = input@meta.data %>%
      droplevels()
    labels = meta[[label_col]]
    expr = Seurat::GetAssayData(input)
    coords = input@images$slice1@coordinates %>%
      dplyr::rename(coord_x = imagecol, coord_y = imagerow) %>%
      dplyr::select(coord_x, coord_y) %>%
      # we need to group by here
      mutate(label = labels,
             barcode = colnames(input),
             idx = row_number())
    
    # print default assay
    default_assay = Seurat::DefaultAssay(input)
    message("using default assay: ", default_assay, " ...")
  } else {
    # data frame or matrix plus metadata data frame
    if (is.null(meta)) {
      stop("must provide metadata if not supplying a Seurat object")
    }
    
    # check if input is sparse matrix or numberic matrix/df
    valid_input = is(input, 'sparseMatrix') ||
      is_numeric_matrix(input) ||
      is_numeric_dataframe(input)
    if (!valid_input)
      stop("input must be Seurat, sparse matrix, numeric matrix, or ",
           "numeric data frame")
    
    # extract columns
    expr = input
    meta %<>% droplevels()
    labels = meta[[label_col]]
    
    # get coordinates
    if (is.null(coords)) {
      stop("must provide coordinate information if not supplying a Seurat object")
    }
    if (any(!coord_cols %in% colnames(coords))) {
      stop("coords must have x and y columns. Please check.")
    }
    coords %<>%
      # we need to group by here
      mutate(label = labels,
             barcode = colnames(expr),
             idx = row_number())
  }
  
  # check dimensions are non-zero
  if (length(dim(expr)) != 2 || !all(dim(expr) > 0)) {
    stop("expression matrix has at least one dimension of size zero")
  }
  
  # check dimensions match
  n_cells1 = nrow(meta)
  n_cells2 = ncol(expr)
  if (n_cells1 != n_cells2) {
    stop("number of barcodes in metadata (", n_cells1,
         ") does not match number of barcodes in expression (", n_cells2, ")")
  }
  
  # check at least two labels
  if (n_distinct(labels) == 1) {
    stop("only one label provided: ", unique(labels))
  }
  
  # check for missing labels
  if (any(is.na(labels))) {
    stop("labels contain ", sum(is.na(labels)), "missing values")
  }
  
  # make sure `label` is not a rowname in `input` (feature in RF)
  if ("label" %in% rownames(expr)) {
    warning("row `label` exists in input; changing ...")
    to_fix = which(rownames(expr) == "label")
    rownames(expr)[to_fix] = paste0("label", seq_along(rownames(expr)[to_fix]))
  }
  
  # fix column names that ranger can't handle
  ## (here, they are still rows)
  rf_engine = "randomForest" ## set to randomForest - don't use ranger
  if (classifier == "rf" && rf_engine == "ranger") {
    invalid_rows = any(grepl("-|\\.|\\(|\\)", rownames(expr)))
    if (invalid_rows) {
      warning("classifier `rf` with engine `ranger` cannot handle characters ",
              "-.() in column names; replacing ...")
      expr %<>% set_rownames(gsub("-|\\.|\\(|\\)", "", rownames(.)))
    }
  }
  
  # remove missing values
  missing = is.na(expr)
  if (any(missing)) {
    stop("matrix contains ", sum(missing), "missing values")
  }
  
  # detect mode
  if (is.numeric(labels)) {
    mode = "regression"
    multiclass = F
    
    # throw a warning if there are only three or less values
    if (n_distinct(labels) <= 3) {
      warning("doing regression with only ", n_distinct(labels),
              " unique values")
    }
  } else {
    mode = "classification"
    
    # check whether we are working with multiclass data
    multiclass = n_distinct(labels) > 2
    # check classifier is compatible with classification problem
    if (multiclass & classifier == "lr") {
      stop("multi-class classification with classifier = 'lr' is currently not ",
           "supported in tidymodels `logistic_reg`")
    }
    
    # make sure y is a factor if doing classification
    if (!is.factor(labels)) {
      warning("coercing labels to factor ...")
      labels %<>% as.factor()
    }
  }
  
  # check if showing progress or not
  if (show_progress == TRUE) {
    apply_fun = pbmclapply
  } else {
    apply_fun = mclapply
  }
  
  # check augur mode
  if (augur_mode == 'velocity') {
    # reset feature selection
    message("disabling feature selection for augur_mode=\"velocity\" ...")
    feature_perc = 1
    var_quantile = 1
  } 
  if (augur_mode != 'velocity') {
    # perform feature selection on the entire spatial transcriptome matrix
    expr %<>% select_variance()
    feature_perc = 1
    var_quantile = 1
  }
  
  # define the augur parameters
  augur_params = list(
    label_col = label_col,
    n_subsamples = n_subsamples,
    subsample_size = subsample_size,
    folds = folds,
    var_quantile = var_quantile,
    feature_perc = feature_perc,
    n_threads = 1,
    show_progress = show_progress,
    augur_mode = augur_mode,
    classifier = classifier,
    # random forest parameters
    rf_params = list(trees = rf_params$trees,
                     mtry = rf_params$mtry,
                     min_n = rf_params$min_n,
                     importance = rf_params$importance),
    # logistic regression parameters
    lr_params = list(mixture = lr_params$mixture, penalty = lr_params$penalty)
  )
  
  # get the nearest neighbors based on spatial coordinates
  nn = coords %>%
    dplyr::select(all_of(coord_cols)) %>%
    nn2(k = ncol(expr)) %>%
    extract("nn.idx") %>%
    as.data.frame() %>%
    mutate(idx = row_number()) %>%
    dplyr::select(-1) %>%
    gather(neighbor, row_index, -idx) %>%
    left_join(coords, by = 'idx') %>%
    # add in the neighbors and re-code everything
    dplyr::select(-idx, -all_of(coord_cols)) %>%
    dplyr::rename(source = barcode,
                  source_group = label,
                  idx = row_index) %>%
    left_join(coords %>% dplyr::select(idx, barcode, label),
              by = 'idx') %>%
    dplyr::rename(target = barcode, target_group = label) %>%
    dplyr::select(-idx) %>%
    group_by(source, target_group) %>%
    mutate(rank = seq_len(n())) %>%
    ungroup()
  
  # iterate over barcodes
  barcodes = colnames(expr)
  res = apply_fun(barcodes, mc.cores = n_threads, function(barcode) {
    # grab the barcodes to keep
    barcodes_keep = nn %>%
      filter(source == barcode) %>%
      group_by(target_group) %>%
      filter(rank %in% seq_len(k)) %>%
      pull(target)
    
    # subset the object and prepare for Augur
    expr0 = expr %>% extract(, barcodes_keep)
    meta0 = meta %>% extract(barcodes_keep, )
    
    # add a dummy channel to ensure compatibility with Augur
    meta0$cell_type = 'magellan'
    
    # by default, Augur will be run with default values
    augur_params$input = expr0
    augur_params$meta = meta0
    augur_params$n_threads = 1
    
    # prioritize this barcode
    augur = do.call(calculate_auc, augur_params)
    
    # reformat the results
    augur$X = NULL
    augur$y = NULL
    augur$feature_importance = NULL
    augur$cell_types = NULL
    augur$parameters$k = k
    augur$results %<>% 
      dplyr::rename(barcode = cell_type) %>%
      mutate(barcode = !!barcode)
    augur$AUC %<>% 
      dplyr::rename(barcode = cell_type) %>%
      mutate(barcode = !!barcode)
    return(augur)
  })
  
  # format the final results and return the Magellan class list
  parameters = res[[1]]$parameters
  results = map(res, ~ extract(., 'results')) %>% bind_rows()
  aucs = map(res, ~ extract(., 'AUC')) %>% bind_rows()
  obj = list(
    parameters = parameters,
    results = results,
    AUC = aucs
  )
  return(obj)
}
