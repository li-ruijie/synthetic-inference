#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 00-shared.r -- Shared utility library                                  |
# +------------------------------------------------------------------------+
#
# Provides the functional programming primitives, GLM helpers, and
# environment management used by every script in the simulation pipeline.
# This file is the foundation layer: it clears the global environment,
# sets reproducibility options (seed 6489, JIT level 3), and exports a
# library of pure functions (Map, Reduce, Filter, Compose, Curry) that
# enforce the functional style used throughout.
#
# Role in the paper:
# The paper (Section 6, "Simulation study") states: "All simulations were
# conducted in R 4.6.1 ..., on Windows 11 and Linux (x86_64)." This file ensures
# a clean, reproducible environment for every script that implements the
# theoretical results (the paper's main theorems) and the six-method comparison.
# The GLM helpers (get.glm.dfres, get.glm.np, get.glm.var) extract the
# quantities that the paper's estimating equation and its synthetic
# counterpart operate on: residual df, parameter count p, and DV/IV
# variable names.
#
# Usage:
# source("lib/00-shared.r") at the top of each script. To prevent re-init
# (e.g. when sourcing from another sourced file), set
# .shared.loaded <- TRUE before sourcing.
if(exists(".shared.loaded") && isTRUE(.shared.loaded)) {
    # already loaded.  Skip re-init
} else {
list.ut <- (function()
    {
        # +------------------------------------------------------------+
        # | Utility functions                                        |
        # +------------------------------------------------------------+
        # Display, environment clearing, and formatting helpers.
        # These maintain the clean-room approach: every script starts
        # from a cleared globalenv, ensuring no state leaks between
        # the simulation scenarios (00-simulate.r) and the generators
        # (07-calibrated, 09-symmetric, 11-projection).
        ut.box <- function(v.contents, chr.head = "+", chr.sep = "-",
                           chr.tail = "+", chr.wall = "|", s.align = "l") {
            if(class(v.contents) |>
               match(c("integer", "numeric", "character"),
                     nomatch = 0) |>
               (any |> Negate())())
                stop("Incorrect contents format.")

            s.contents.lenmax <- nchar(v.contents) |> max()
            s.box.lenmax <- length(c(chr.head, chr.sep, chr.tail)) |> sum()
            if(s.contents.lenmax < s.box.lenmax)
                s.contents.lenmax <- s.box.lenmax
            s.contents.lenpadded <- s.contents.lenmax + 2
            s.contents.lenpadded <- (4 - (s.contents.lenpadded %% 4)) + s.contents.lenpadded

            s.contents.lenpadded.use <-
                switch(s.align,
                       l = -s.contents.lenpadded,
                       r = s.contents.lenpadded,
                       c = s.contents.lenpadded/2)
            v.messages <-
                switch((s.align == "c") + 1,
                       sprintf(paste0(chr.wall, " %*s ", chr.wall),
                               s.contents.lenpadded.use, v.contents),
                       {
                         l.contents <-
                           Map(function(x)
                               {
                                 s.len <- length(x) %/% 2
                                 Map(function(v) paste0(v, collapse = ""),
                                     list(head(x, s.len), tail(x, -s.len)))
                               }, strsplit(v.contents, ""))

                         mapply(function(x)
                                {
                                  paste0(sprintf(paste0(chr.wall, " %*s"),
                                                 s.contents.lenpadded.use, x[[1]]),
                                         sprintf(paste0("%*s ", chr.wall),
                                                 -s.contents.lenpadded.use, x[[2]]),
                                         collapse = "")
                                }, l.contents)
                       })
            v.sep <- paste0(c(chr.head,
                              rep(chr.sep,
                                  s.contents.lenpadded +
                                    ((nchar(c(" ", chr.wall)) * 2) |> sum()) -
                                    (nchar(c(chr.head, chr.tail))  |> sum())),
                              chr.tail),
                            collapse = "")
            mapply(message, c(v.sep, v.messages, v.sep)) |> invisible()
        }
        ut.clear.disp <- function() {
            if(!is.null(dev.list())) {
                invisible(mapply(dev.off, dev.list()))
            }
            switch(commandArgs()[1],
                   RStudio = cat("\014"),
                   system("clear"))
        }
        ut.clear.dispdev <- function() {
            if(!is.null(dev.list()))
                ut.silent(mapply(dev.off, dev.list()))
        }
        ut.clear.globalwspc <- function() {
            local(rm(list = ls()), envir = .GlobalEnv)
        }
        ut.get.all <- function() {
            env <- parent.frame(1)
            mget(ls(envir = env), envir = env)
        }
        ut.suppress <- function(...) {
            suppressMessages(...) |> suppressWarnings()
        }
        ut.silent <- function(...) {
            sink(nullfile())
            try(..., silent = TRUE, outFile = nullfile()) |> ut.suppress()
            sink()
        }
        ut.permargs <- function(...) {
            l.args <- list(...) |> rev()
            Map(function(l, v.ind) l[v.ind],
                l.args,
                mapply(function(v)
                       seq_along(v),
                       l.args,
                       SIMPLIFY = FALSE) |>
                expand.grid())           |>
                  rev()
        }
        ut.save.all <- function() {
            Map(function(n)
                save(list = n,
                     file = paste0(n, ".RData"),
                     envir = .GlobalEnv),
                grep("^df|^list",
                     ls(".GlobalEnv"), value = TRUE))
        }
        ut.time <- function() {
            mapply(function(x)
                   format(Sys.time(), x),
                   c("Date: %Y%m%d %z", "Time: %H%M%S")) |>
                  ut.box()
        }
        # +------------------------------------------------------------+
        # | Parallel computing helpers                               |
        # +------------------------------------------------------------+
        # Used by the Monte Carlo replication scripts (18-mc-projection.r
        # among them) which run 100 independent reps to
        # confirm Theorem 7's structural exactness (sd = 0.0 across
        # replications). PSOCK clusters are used on Windows.
        make.cl <- function(int.cores = parallel::detectCores(),
                            env = environment()) {
            parallel::clusterExport(cl <- parallel::makeCluster(int.cores),
                                    ls(envir = env), envir = env)
            cl
        }
        max.cores <- function() {
            parallel::detectCores()
        }
        set.env.name <- function(x) {
            cal.parent <- as.character(sys.call(-1))[1]
            env.parent <- sys.frame(-1)
            s.name <- paste0(cal.parent, "()")
            invisible(structure(env.parent,
                                env.name = s.name))
        }
        # +------------------------------------------------------------+
        # | Random number generation                                 |
        # +------------------------------------------------------------+
        # gen.rand.int: generates seeds for MC replications.
        # cor.m: generates a random vector with a specified correlation
        # to existing columns. This is the Cholesky-based approach
        # used in 00-simulate.r to produce the correlated latent
        # normals that underlie all 10 DGP scenarios. The known
        # correlation structure ensures the true regression coefficients
        # are significant, which is necessary to test whether generators
        # preserve them (the central question).
        gen.rand.int <- function(n) {
            s.int.max <- .Machine[["integer.max"]]
            sample(-s.int.max:s.int.max, n, replace = FALSE)
        }
        cor.m <- function(m, v.rho) {
            # https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
            s.n <- nrow(m)
            if(is.null(s.n)) {
                s.n <- length(m)
                if(is.null(s.n)) {
                    stop("Error")
                } else {
                    m <- matrix(m, ncol = 1)
                }
            }
            m.sc <- scale(m)
            l.m.svd <- svd(m.sc)
            svd.u <- l.m.svd[["u"]]
            svd.v <- l.m.svd[["v"]]
            svd.d <- l.m.svd[["d"]]
            v.svd.d.dual <- ifelse(svd.d > 0, 1/svd.d, 0)
            if(length(svd.d) == 1) {
                m.svd.d.dual <- matrix(v.svd.d.dual, ncol = 1, nrow = 1)
            } else {
                m.svd.d.dual <- diag(v.svd.d.dual)
            }
            m.sc.dual <- (s.n - 1) * svd.u %*% m.svd.d.dual %*% t(svd.v)
            v.res <- residuals(lm(rnorm(s.n) ~ m.sc))
            s.sigma.sq <- c((1 - (v.rho %*% cov(m.sc.dual) %*% v.rho)) / var(v.res))
            if(s.sigma.sq < 0) {
                stop("Correlations provided not possible.")
            }
            m.sc.dual %*% v.rho + sqrt(s.sigma.sq) * v.res
        }
        # +------------------------------------------------------------+
        # | Functional programming constructs                        |
        # +------------------------------------------------------------+
        # Compose, Curry, Map, Reduce, Filter -- the primitives the
        # pipeline style used in all scripts is built from. These enforce
        # the pure-function, no-mutation discipline (R coding convention #5).
        #
        # Map variants (TT/TF/FT/FF) control simplify and USE.names
        # to avoid the pitfalls of sapply. Map.FT (= Map) is the
        # default: list output with names preserved (convention #1).
        Compose <- function(...) {
            func.list <- list(...)
            if(!all(unlist(Map(is.function, func.list))))
                stop("Argument is not a function")
            function(...) {
                Reduce(function(x, f) f(x), func.list, ...)
            }
        }
        ComposeM <- function(...) {
            fs <- list(...)
            if(!all(unlist(Map(is.function, fs))))
                stop("Argument is not a function")
            fs.tail <- fs[-1]
            fs.head <- fs[[1]]
            function(...) {
                Reduce(function(x, f) f(x), fs.tail, fs.head(...))
            }
        }
        Curry <- function(f, ...) {
            .orig <- list(...)
            function(...) {
                do.call(f, c(.orig, list(...)))
            }
        }
        CurryL <- function(f, ...) {
            .curried <- as.list(match.call())[c(-1, -2)]
            function(...) {
                .args <- as.list(match.call())[-1]
                eval(substitute(do.call(f, c(.curried, .args))))
            }
        }
        Filter <- function(f, x, nomatch = NA_integer_) {
            `if`(any((ind <- (unlist(Map(f, x))) > 0)),
                 x[ind], nomatch)
        }
        Find <- function(f, x, right = FALSE, nomatch = NULL) {
            `if`(any((pos <- Position(f, x, right, nomatch = 0L))),
                 x[[pos]], nomatch)
        }
        Map.TT <- function(f, ...) {
            mapply(FUN = match.fun(f), ...,
                   SIMPLIFY  = TRUE,
                   USE.NAMES = TRUE)
        }
        Map.TF <- function(f, ...) {
            mapply(FUN = match.fun(f), ...,
                   SIMPLIFY  = TRUE,
                   USE.NAMES = FALSE)
        }
        Map.FT <- function(f, ...) {
            mapply(FUN = match.fun(f), ...,
                   SIMPLIFY  = FALSE,
                   USE.NAMES = TRUE)
        }
        Map.FF <- function(f, ...) {
            mapply(FUN = match.fun(f), ...,
                   SIMPLIFY  = FALSE,
                   USE.NAMES = FALSE)
        }
        Map <- Map.FT
        Negate.v <- function(...) {
            !(...)
        }
        Position <- function(f, x, right = FALSE,
                             nomatch = NA_integer_) {
            f <- match.fun(f)
            ind <- seq_along(x)
            if(right) ind <- rev(ind)
            `for`(i, ind, if(f(x[[i]])) return(i))
            nomatch
        }
        Swap <- function(f, pos = 2L) {
            f <- match.fun(f)
            pos.coerce <- as.integer(pos) |> suppressWarnings()
            s <- `if`(is.na(pos.coerce),
                      {
                        paste0("function(arg1, ...) f(", pos, " = arg1, ...)", collapse = "")
                      },
                      {
                        if(!all.equal(pos, pos.coerce))
                          stop("pos needs to be either an integer or character string.")
                        paste0("function(",
                               paste0("arg", c(pos, seq(pos.coerce - 1)), collapse = ", "),
                               ", ...)", " f(",
                               paste0("arg", seq(pos), collapse = ", "),
                               ", ...)", collapse = "")
                      })
            structure(parse(text = s,
                            keep.source = FALSE)[[1]] |>
                      eval(),
                    func.name = all.names(sys.call())[-1],
                    func = substitute(f))
        }
        # +------------------------------------------------------------+
        # | Naming and subsetting helpers                            |
        # +------------------------------------------------------------+
        # cnms/nms/rnms: functional replacement forms for colnames,
        # names, rownames (convention #7: no side-effect assignment).
        # ss.b/ss.bb: `[` and `[[` as proper functions for use in
        # Map/Reduce pipelines (convention #6: [["name"]] over $).
        cnms <- function(v, v.nm) {
            `colnames<-`(v, v.nm)
        }
        nms <- function(v, v.nm) {
            `names<-`(v, v.nm)
        }
        nms.sort <- function(obj, ...) {
            `if`((v.nm <- names(obj)) |> is.null(),
                 obj, obj[v.nm |> sort(...)])
        }
        rnms <- function(v, v.nm) {
            `rownames<-`(v, v.nm)
        }
        ## class checks
        is.0 <- function(v) {
            !is.null(v) & (length(v) == 0)
        }
        is.v <- function(obj, v.class, simplify = any) {
            # multivariate version of is
            Map(function(s.class)
                is(obj, s.class),
                v.class) |>
                  unlist()     |>
                  match.fun(simplify)()
        }
        ## flow control & comparisons
        identical.v <- function(v) {
            if(length(v) == 1)
                return(TRUE)
            bv.na <- is.na(v)
            if(all(bv.na)) {
                warning("Input contains NA.")
                return(TRUE)
            }
            if(any(bv.na)) {
                warning("Input contains NA.")
                return(FALSE)
            }
            all(v[1] == v)
        }
        if.t <- function(b.test, v.true, v.false) {
            # for use in pipes
            # if, binary and ternary forms
            if(missing(v.false)) {
                if(b.test) {
                    v.true
                }
            } else {
                if(b.test) {
                    v.true
                } else {
                    v.false
                }
            }
        }
        ## subsetting
        ss.b <- function(v, i) {
            `[`(v, i)
        }
        ss.bb <- function(v, i) {
            `[[`(v, i)
        }
        ss.e.nm <- function(obj, v.exclude.names) {
            # subset, excluding by names
            v.obj.names <- names(obj)
            if(is.null(v.obj.names)) {
                stop("Input object has no names.")
            } else {
                b.include <- match(v.obj.names,
                                   v.exclude.names,
                                   nomatch = 0) == 0
                b.exclude <- !b.include
                if(all(b.exclude))
                    stop("All items excluded.")
                if(any(b.include))
                    obj[b.include]
            }
        }
        ssa.b <- function(v, i, v.rep) {
            `[<-`(v, i, v.rep)
        }
        ssa.bb <- function(v, i, v.rep) {
            `[[<-`(v, i, v.rep)
        }
        ## attributes
        attr.exa <- function(obj) {
            attributes(obj)
        }
        attr.ex1 <- function(obj, attr.nm) {
            attr.exa(obj)[[attr.nm]]
        }
        attr.set <- function(obj, l.attr) {
            do.call(structure, c(obj, l.attr))
        }
        attr2env <- function(obj, env = parent.frame()) {
            l <- attributes(obj)
            b.include <- match(names(l), c("dim", "names"), nomatch = 0) == 0
            if(any(b.include))
                list2env(l[b.include], envir = env)
        }
        set.list.index <- function(l) {
            # create index
            l.index <- list(i = seq_along(l))
            v.names <- names(l)
            if(!is.null(v.names))
                l.index <- c(l.index, list(names = v.names))

            # index is now a list of lists of the indexes.
            # - this is a linked list and will need to be restructured
            # - restructure to be 1 list for each item present in l containing the attributes.
            l.index.rstruc <- do.call(Map,
                                      c(list(list),
                                        Map(as.list, l.index)))

            # attach the attributes
            mapply(function(i, l.attr)
                   structure(i, index = l.attr),
                   l, l.index.rstruc,
                   SIMPLIFY = FALSE,
                   USE.NAMES = TRUE)
        }
        # +------------------------------------------------------------+
        # | GLM extraction helpers                                   |
        # +------------------------------------------------------------+
        # These extract the key quantities from fitted GLMs that the
        # paper's framework (Section 2) operates on:
        #   - get.glm.dfres: residual degrees of freedom (n - p),
        #     used in dispersion estimation (Proposition 1(b)).
        #   - get.glm.np: sample size n and parameter count p,
        #     critical for the hat matrix H = X(X'X)^{-1}X' in the
        #     projection (Theorem 7, Corollary: Gaussian closed-form).
        #   - get.glm.var: dependent and independent variable names,
        #     needed to construct the model matrix X for the estimating
        #     equation U_n(beta) = 0.
        get.glm.dfres <- function(mo) {
            summary(mo)[["df"]][2]
        }
        get.glm.np <- function(mo) {
            v.df <- summary(mo)[["df"]]
            list(n = sum(v.df[2:3]), p = v.df[3])
        }
        get.glm.var <- function(mo) {
            v <- formula(mo) |> all.vars()
            list(dv = v[1], iv = v[-1])
        }
        # +------------------------------------------------------------+
        # | Summary statistics helpers                               |
        # +------------------------------------------------------------+
        # table.v/table.vec: tabulation used in the report scripts
        # (14-report-inferential.r among them) to
        # compute stratum sizes n_s and prevalences pi_s = n_s/n.
        # These are the quantities in Theorem 5 (sparse strata MSE)
        # and the shrinkage weight w_s = tau^2 / (tau^2 + sigma^2/n_s).
        table.v <- function(v, p = FALSE) {
            v <- table(v)
            v.n <- names(v)
            if(p) {
                v <- prop.table(v)
                v.n <- paste0(v.n, ".p")
            }
            `names<-`(as.vector(v), v.n)
        }
        table.vec <- function(v, prop = FALSE) {
            l.fn <- list(as.data.frame)
            l.fn <- c(`if`(prop, c(l.fn, list(prop.table)), l.fn), list(table))
            d <- Reduce(function(f, ...) f(...), l.fn, v, right = TRUE)
            `names<-`(d[,2], as.character(d[,1]))
        }
        # +------------------------------------------------------------+
        # | Initialisation: packages, options, seed                  |
        # +------------------------------------------------------------+
        # JIT level 3 for bytecode compilation of all R code.
        # Master seed 8 (default).  MC scripts override to 6489 per the
        # paper's reproducibility statement (Section 6).
        # contr.treatment for unordered factors ensures the model
        # matrix X uses dummy coding, matching the paper's notation
        # x_i = (1, age, female, race2, race3, income)'.
        library(compiler)
        # options
        enableJIT(3)
        set.seed(8)
        options(list(stringsAsFactors = FALSE,
                     save.defaults = list(compress = "xz",
                                          compression_level = 9),
                     contrasts = c(unordered = "contr.treatment",
                                   ordered = "contr.treatment"),
                     max.print = 80))
        # save functions
        list.ut <- ut.get.all()

        # clear workspace
        ut.clear.globalwspc()
        # turn off display devices
        ut.clear.dispdev()
        # import functions
        list2env(list.ut, envir = .GlobalEnv)

        # clear environment and show time
        ut.clear.globalwspc ()
        ut.clear.disp()
        ut.time()

        # return functions
        list.ut
    })()
} # end if(.shared.loaded) guard
