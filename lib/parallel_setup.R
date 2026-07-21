# ══════════════════════════════════════════════════════════════
#  parallel_setup.R: Portable parallel backend for R
#
#  Works on Linux/macOS (fork via mcmapply) and Windows (PSOCK
#  via clusterMap). Source this file at the top of any script
#  that needs parallelism.
#
#  Primary API: par_map (parallel Map, multivariate)
#               par_mapply (par_map + simplify)
#               par_seeds (generate per-worker seeds)
#
#  lapply/sapply and their parallel wrappers are intentionally
#  excluded. Use Map/mapply (and par_map/par_mapply) exclusively.
#
#  Style: functional (no $, no side-effect assignment, no for).
#  All extraction via [["name"]]. All modification via
#  functional replacement forms (`[[<-`, `class<-`, etc.).
#
#  Usage:
#    source("lib/parallel_setup.R")
#    cfg <- setup_parallel(n_cores = 4)
#    seeds <- par_seeds(6489L, length(items))
#    results <- par_map(cfg, function(s, x) {
#        set.seed(s)
#        # ... work ...
#    }, seeds, items)
#    cleanup_parallel(cfg)
#
#  Worker initialisation (PSOCK + fork):
#    cfg <- setup_parallel(
#        n_cores = 4,
#        packages = c("lavaan", "semTools"),
#        init_worker = function() {
#            Sys.setenv(OPENBLAS_NUM_THREADS = "1")
#        }
#    )
#    # On PSOCK: init_worker runs via clusterCall after packages
#    # On fork:  init_worker runs once in the parent before forking
#    #           (forked children inherit the parent's environment)
# ══════════════════════════════════════════════════════════════

suppressPackageStartupMessages(library(parallel))

# ── Backend detection ─────────────────────────────────────────

.is_fork_safe <- function() {

    # fork works on Unix-like systems but not inside RStudio on
    # macOS (risk of GUI deadlocks). WSL counts as Linux here.
    os <- .Platform[["OS.type"]]
    if (os != "unix") return(FALSE)

    # Respect the user/system override
    env_val <- Sys.getenv("R_PARALLEL_BACKEND", unset = "")
    if (nzchar(env_val)) return(tolower(env_val) == "fork")

    # Disable fork inside RStudio on macOS to avoid crashes
    if (Sys.info()[["sysname"]] == "Darwin" &&
        nzchar(Sys.getenv("RSTUDIO", unset = ""))) {
        return(FALSE)
    }

    TRUE
}

# ── Seed generation ──────────────────────────────────────────

#' Generate n unique per-worker seeds from a master seed.
#'
#' Uses the default RNG (Mersenne-Twister). Saves and restores the
#' global RNG state so calling par_seeds has no side effect on the
#' caller's random stream.
#'
#' @param seed  Integer master seed.
#' @param n     Number of unique seeds to generate.
#' @return Integer vector of length n.
par_seeds <- function(seed, n) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
    } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(seed)
    sample.int(.Machine[["integer.max"]], size = n)
}

# ── Setup ─────────────────────────────────────────────────────

#' Initialise the parallel backend.
#'
#' @param n_cores     Number of worker cores. NULL = detectCores() - 1.
#' @param packages    Character vector of packages to load on workers.
#'                    PSOCK: loaded via library() on each worker.
#'                    Fork: loaded in the parent (children inherit).
#' @param export      Character vector of object names to export to PSOCK
#'                    workers. Looked up in the calling environment.
#'                    Ignored for fork (children share parent memory).
#' @param init_worker A function (no arguments) to run on each worker
#'                    for arbitrary initialisation (e.g., Sys.setenv,
#'                    options, .libPaths). PSOCK: runs via clusterCall.
#'                    Fork: runs once in the parent before forking.
#'                    NULL (default) = no extra initialisation.
#' @return A par_cfg list with [["type"]], [["n_cores"]], [["cl"]].
#'
#' @details Seeds are not managed by setup_parallel. Use par_seeds()
#'   to generate per-worker seeds, then pass them as an iterable to
#'   par_map. Each worker calls set.seed() with its own seed. This
#'   uses the default RNG (Mersenne-Twister) throughout, with no
#'   L'Ecuyer-CMRG.
setup_parallel <- function(n_cores     = NULL,
                           packages    = character(0),
                           export      = character(0),
                           init_worker = NULL) {

    max_cores <- parallel::detectCores(logical = TRUE)
    n_cores <- if (is.null(n_cores)) max(1L, max_cores - 1L) else n_cores
    n_cores <- n_cores |> as.integer() |> min(max_cores) |> max(1L)

    use_fork <- .is_fork_safe() && n_cores > 1L

    backend <- if (n_cores == 1L) "sequential"
               else if (use_fork) "fork"
               else "psock"

    cfg <- `class<-`(list(
        type    = backend,
        n_cores = n_cores,
        cl      = NULL
    ), "par_cfg")

    if (backend == "sequential") {
        message("[parallel] Running sequentially (1 core).")
        if (length(packages) > 0L) {
            invisible(Map(
                function(p) library(p, character.only = TRUE),
                packages
            ))
        }
        if (!is.null(init_worker)) init_worker()
        return(cfg)
    }

    if (backend == "fork") {
        message(sprintf("[parallel] Using fork backend with %d cores.",
                        n_cores))
        if (length(packages) > 0L) {
            invisible(Map(
                function(p) library(p, character.only = TRUE),
                packages
            ))
        }
        if (!is.null(init_worker)) init_worker()
        return(cfg)
    }

    # PSOCK path
    message(sprintf("[parallel] Using PSOCK backend with %d cores.",
                    n_cores))
    cl <- parallel::makeCluster(n_cores)

    # Load packages on workers
    if (length(packages) > 0L) {
        parallel::clusterCall(cl, function(pkgs) {
            invisible(Map(
                function(p) library(p, character.only = TRUE),
                pkgs
            ))
        }, pkgs = packages)
    }

    # Arbitrary worker initialisation (env vars, options, etc.)
    if (!is.null(init_worker)) {
        parallel::clusterCall(cl, init_worker)
    }

    # Export objects from the caller's environment
    if (length(export) > 0L) {
        envir <- parent.frame()
        parallel::clusterExport(cl, varlist = export, envir = envir)
    }

    `class<-`(`[[<-`(cfg, "cl", cl), "par_cfg")
}

# ── Cleanup ───────────────────────────────────────────────────

#' Stop the PSOCK cluster if one exists. Safe to call multiple times.
cleanup_parallel <- function(cfg) {
    if (!is.null(cfg[["cl"]])) {
        tryCatch(
            parallel::stopCluster(cfg[["cl"]]),
            error = function(e) NULL
        )
    }
    invisible(NULL)
}

# ── Internal: wrap fun for safe execution + progress ──────────

.make_safe_fun <- function(FUN) {
    function(...) {
        tryCatch(
            list(error = FALSE, value = FUN(...)),
            error = function(e) list(error = TRUE, message = conditionMessage(e))
        )
    }
}

.make_progress_fun <- function(FUN, prog_file, n) {
    safe_fun <- .make_safe_fun(FUN)
    function(...) {
        res <- safe_fun(...)
        tryCatch({
            cur <- prog_file |> readLines(warn = FALSE) |> as.integer()
            writeLines(as.character(cur + 1L), prog_file)
            if ((cur + 1L) %% max(1L, n %/% 20L) == 0L || (cur + 1L) == n) {
                message(sprintf("[parallel] %d / %d (%.0f%%)",
                                cur + 1L, n, 100 * (cur + 1L) / n))
            }
        }, error = function(e) NULL)
        res
    }
}

.prepare_fun <- function(FUN, n, .progress) {
    if (.progress) {
        prog_file <- tempfile("par_progress_")
        writeLines("0", prog_file)
        list(
            run_fun   = .make_progress_fun(FUN, prog_file, n),
            prog_file = prog_file
        )
    } else {
        list(
            run_fun   = .make_safe_fun(FUN),
            prog_file = NULL
        )
    }
}

# ── Primary dispatch: par_map (parallel Map) ──────────────────

#' Parallel Map, the primary parallel dispatch.
#'
#' Applies fun element-wise over one or more iterables, using
#' mcmapply (fork), clusterMap (PSOCK), or Map (sequential).
#'
#' For reproducibility, generate per-worker seeds with par_seeds()
#' and pass them as the first iterable:
#'
#'   seeds <- par_seeds(6489L, length(items))
#'   par_map(cfg, function(s, x) { set.seed(s); ... }, seeds, items)
#'
#' @param cfg       A par_cfg object from setup_parallel().
#' @param fun       Function to apply. Receives one element from each
#'                  iterable as positional arguments.
#' @param ...       One or more vectors/lists to iterate over in
#'                  parallel (passed as the dots of Map/mcmapply/clusterMap).
#' @param .progress Logical.  If TRUE, print a progress counter.
#' @return A list the same length as the first iterable.
par_map <- function(cfg, FUN, ..., .progress = FALSE) {

    dots <- list(...)
    n <- length(dots[[1L]])
    prep <- .prepare_fun(FUN, n, .progress)

    results <- switch(cfg[["type"]],
        fork = {
            parallel::mcmapply(
                prep[["run_fun"]], ...,
                mc.cores = cfg[["n_cores"]],
                mc.set.seed = FALSE,
                SIMPLIFY = FALSE
            )
        },
        psock = {
            parallel::clusterExport(
                cfg[["cl"]],
                varlist = c("FUN"),
                envir = environment()
            )
            if (.progress) {
                prog_file <- prep[["prog_file"]]
                n_val <- n
                parallel::clusterExport(
                    cfg[["cl"]],
                    varlist = c("prog_file", "n_val"),
                    envir = environment()
                )
            }
            parallel::clusterMap(
                cfg[["cl"]],
                prep[["run_fun"]], ...,
                SIMPLIFY = FALSE
            )
        },
        sequential = Map(prep[["run_fun"]], ...)
    )

    if (!is.null(prep[["prog_file"]])) unlink(prep[["prog_file"]])
    results
}

#' Parallel mapply, which is like par_map but simplifies the result.
par_mapply <- function(cfg, FUN, ..., .progress = FALSE) {
    res <- par_map(cfg, FUN, ..., .progress = .progress)
    vals <- Map(function(r) {
        if (isTRUE(r[["error"]])) NA else r[["value"]]
    }, res)
    tryCatch(simplify2array(vals), error = function(e) vals)
}

# ── Convenience wrappers ──────────────────────────────────────

#' Chunked parallel Map, which splits a data frame into n_cores
#' chunks, sends one chunk per worker, and Reduce(rbind)s the
#' results. Ideal for PSOCK where you want to minimise serialisation.
#'
#' fun receives a data.frame chunk and must return a data.frame.
par_map_chunked <- function(cfg, df, FUN, ..., .progress = FALSE) {
    n <- nrow(df)
    idx <- seq_len(n) |> cut(breaks = cfg[["n_cores"]], labels = FALSE)
    chunks <- split(df, idx)

    results <- par_map(cfg, function(chunk, ...) {
        FUN(chunk, ...)
    }, chunks, .progress = .progress)

    # Unwrap safe_fun structure and Reduce(rbind)
    results |>
        Map(f = function(r) {
            if (isTRUE(r[["error"]])) {
                warning("Chunk failed: ", r[["message"]])
                NULL
            } else {
                r[["value"]]
            }
        }) |>
        Filter(f = Negate(is.null)) |>
        Reduce(f = rbind)
}

# ── Utilities ─────────────────────────────────────────────────

#' Print method for par_cfg
print.par_cfg <- function(x, ...) {
    cat(sprintf("Parallel config: %s backend, %d core(s)\n",
                x[["type"]], x[["n_cores"]]))
    invisible(x)
}

#' Quick benchmark: compare sequential vs parallel on a test function.
#' Returns a data.frame of timings.
par_benchmark <- function(cfg, X, FUN, ..., times = 3L) {
    seq_cfg <- `class<-`(list(
        type = "sequential", n_cores = 1L, cl = NULL
    ), "par_cfg")

    timings <- Map(function(run) {
        t_seq <- system.time(par_map(seq_cfg, FUN, X, ...))[["elapsed"]]
        t_par <- system.time(par_map(cfg, FUN, X, ...))[["elapsed"]]
        data.frame(
            backend = c("sequential", cfg[["type"]]),
            run     = c(run, run),
            elapsed = c(t_seq, t_par)
        )
    }, seq_len(times)) |> Reduce(f = rbind)

    agg <- aggregate(elapsed ~ backend, data = timings, FUN = function(x) {
        c(mean = mean(x), sd = sd(x))
    })

    message("[benchmark] Mean elapsed times:")
    invisible(Map(function(i) {
        message(sprintf("  %-12s %.3fs (sd %.3fs)",
                        agg[["backend"]][i], agg[["elapsed"]][i, "mean"],
                        agg[["elapsed"]][i, "sd"]))
    }, seq_len(nrow(agg))))

    timings
}
