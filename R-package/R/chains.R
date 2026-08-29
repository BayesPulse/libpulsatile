#-------------------------------------------------------------------------------
# chains.R - Extracting and inventorying posterior MCMC chains
#-------------------------------------------------------------------------------

#' Posterior MCMC chains
#'
#' @description Functions for locating and extracting the various levels of
#'   posterior MCMC chains held in a model fit.
#'
#'   \code{chains()} prints an inventory of every chain in a fit, the
#'   expression that extracts it, and its dimensions. The single-subject model
#'   has accessor functions (\code{patient_chain()}, \code{pulse_chain()});
#'   the population and joint models store their chains as list elements, and
#'   the population per-pulse chains are nested two levels deep (subject, then
#'   saved draw), so \code{chains()} is the quickest way to see what a fit
#'   holds and how to reach it.
#'
#' @param fit A model fit from \code{\link{fit_pulse}},
#'   \code{\link{fit_pulse_population}}, or \code{\link{fit_pulse_joint}}.
#'   Accessors \code{patient_chain()} and \code{pulse_chain()} currently
#'   support single-subject (\code{pulse_fit}) objects only.
#' @return \code{chains()} prints the inventory and invisibly returns it as a
#'   data frame with one row per chain and columns \code{chain} (name),
#'   \code{access} (the expression that extracts it, written in terms of an
#'   object named \code{fit}), \code{rows}, and \code{columns}. Accessors in
#'   the \code{access} column that contain \code{[[i]]} or \code{[[j]]} are
#'   subscripted by subject and saved draw respectively.
#'
#'   \code{patient_chain()} returns the subject-level (common parameter) chain
#'   and \code{pulse_chain()} the per-pulse chain, each as a tibble.
#' @import tibble
#' @keywords pulse fit
#' @examples
#'
#' pulse <- simulate_pulse()
#' spec  <- pulse_spec()
#' fit   <- fit_pulse(data = pulse, iters = 1000, thin = 10,
#'                    burnin = 100, spec = spec)
#' chains(fit)
#' head(patient_chain(fit))
#' head(pulse_chain(fit))
#'
#' @export
chains <- function(fit) UseMethod("chains")


#' @export
chains.default <- function(fit) {
  stop(paste("chains() requires a model fit: a pulse_fit from fit_pulse(),",
             "a population_fit from fit_pulse_population(), or a joint_fit",
             "from fit_pulse_joint()."))
}


#' @export
chains.pulse_fit <- function(fit) {
  inventory <- rbind(
    chain_entry("patient", "patient_chain(fit)", fit$patient_chain),
    chain_entry("pulse",   "pulse_chain(fit)",   fit$pulse_chain))
  print_chain_inventory(inventory, "Single-subject fit (pulse_fit)")
}


#' @export
chains.population_fit <- function(fit) {
  inventory <- rbind(
    chain_entry("population", "fit$population_chain", fit$population_chain),
    chain_entry("subject", "fit$subject_chains[[i]]", fit$subject_chains[[1]]),
    chain_entry("pulse", "fit$pulse_chains[[i]][[j]]",
                fit$pulse_chains[[1]][[1]]))
  print_chain_inventory(
    inventory,
    paste0("Population fit (population_fit), ", fit$num_subjects, " subjects"),
    paste0("i = subject (1 to ", fit$num_subjects, "); ",
           "j = saved draw (1 to ", length(fit$pulse_chains[[1]]), "). ",
           "Subject and pulse rows report the dimensions of i = j = 1."))
}


#' @export
chains.joint_fit <- function(fit) {
  inventory <- rbind(
    chain_entry("driver",         "fit$driver_chain",      fit$driver_chain),
    chain_entry("response",       "fit$response_chain",    fit$response_chain),
    chain_entry("association",    "fit$association_chain", fit$association_chain),
    chain_entry("driver pulse",   "fit$driver_pulse_chain",
                fit$driver_pulse_chain),
    chain_entry("response pulse", "fit$response_pulse_chain",
                fit$response_pulse_chain))
  print_chain_inventory(
    inventory, "Joint driver-response fit (joint_fit)",
    paste("The association chain holds the coupling draws (rho, nu).",
          "Pulse chains are stacked across saved draws and carry no",
          "iteration column."))
}


# One inventory row for a single chain.
chain_entry <- function(chain, access, obj) {
  data.frame(chain = chain, access = access,
             rows = nrow(obj), columns = ncol(obj),
             stringsAsFactors = FALSE)
}


# Print an inventory as a left-aligned table, then return it invisibly.
print_chain_inventory <- function(inventory, header, footer = NULL) {

  cell <- rbind(c("chain", "access", "rows", "columns"),
                cbind(inventory$chain, inventory$access,
                      format(inventory$rows), format(inventory$columns)))
  widths <- apply(nchar(cell), 2, max)
  # Every column but the last is padded; padding the last would only leave
  # trailing whitespace on each line.
  for (k in seq_len(ncol(cell) - 1)) {
    cell[, k] <- formatC(cell[, k], width = -widths[k])
  }

  cat("\n", header, "\n\n", sep = "")
  for (i in seq_len(nrow(cell))) {
    cat("  ", paste(cell[i, ], collapse = "  "), "\n", sep = "")
  }
  if (!is.null(footer)) {
    cat("\n")
    cat(strwrap(footer, width = 76, indent = 2, exdent = 2), sep = "\n")
  }
  cat("\n")

  invisible(inventory)
}


#' @rdname chains
#' @export
pulse_chain <- function(fit) UseMethod("pulse_chain")

#' @export
pulse_chain.pulse_fit <- function(fit) {
  fit$pulse_chain
}

#' @rdname chains
#' @export
patient_chain <- function(fit) UseMethod("patient_chain")

#' @export
patient_chain.pulse_fit <- function(fit) {
  fit$patient_chain
}


#-------------------------------------------------------------------------------
# End of file
#-------------------------------------------------------------------------------
