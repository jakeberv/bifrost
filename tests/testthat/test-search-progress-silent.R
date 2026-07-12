test_that("public search disables heartbeats throughout all three stages", {
  testthat::skip_if_not_installed("ape")
  testthat::skip_if_not_installed("phytools")
  testthat::skip_if_not_installed("mvMORPH")

  fit_count <- 0L
  events <- new.env(parent = emptyenv())
  events$sessions <- 0L
  events$created <- character()
  progress_session <- .bifrost_search_progress_session
  testthat::local_mocked_bindings(
    fitMvglsAndExtractGIC.formula = function(formula, tree, trait_data, ...) {
      fit_count <<- fit_count + 1L
      list(
        model = list(corrSt = list(phy = tree)),
        GIC = list(GIC = 1000 - fit_count)
      )
    },
    removeShiftFromTree = function(tree, shift_node, stem = FALSE) tree,
    .bifrost_search_progress_session = function(enabled) {
      events$sessions <- events$sessions + 1L
      renderer <- list(
        create = function(label, total) {
          events$created <- c(events$created, label)
          stop("disabled progress must not create a row")
        },
        update = function(...) stop("disabled progress must not update a row"),
        output = function(...) stop("disabled progress must not emit row output"),
        done = function(...) stop("disabled progress must not finalize a row")
      )
      progress_session(enabled, renderer)
    }
  )

  set.seed(20260714)
  tree <- ape::rtree(8)
  traits <- matrix(stats::rnorm(16), ncol = 2)
  rownames(traits) <- tree$tip.label

  messages <- testthat::capture_messages(
    result <- suppressWarnings(searchOptimalConfiguration(
      baseline_tree = tree,
      trait_data = traits,
      min_descendant_tips = 3,
      num_cores = 1,
      shift_acceptance_threshold = -Inf,
      uncertaintyweights = TRUE,
      uncertaintyweights_par = FALSE,
      plot = FALSE,
      IC = "GIC",
      store_model_fit_history = FALSE,
      verbose = FALSE,
      progress = FALSE,
      method = "LL"
    ))
  )

  testthat::expect_length(messages, 0L)
  testthat::expect_identical(events$sessions, 1L)
  testthat::expect_length(events$created, 0L)
  testthat::expect_false(result$user_input$progress)
  testthat::expect_gt(length(result$shift_nodes_no_uncertainty), 0L)
  testthat::expect_equal(
    nrow(result$ic_weights),
    length(result$shift_nodes_no_uncertainty)
  )
})
