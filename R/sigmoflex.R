#------------------------------------------------------------------------------#
#' @import ggplot2
#------------------------------------------------------------------------------#
NULL

#------------------------------------------------------------------------------#
#' Sigmoflex
#'
#' Implementation of a proposal to quantify the differences between disease
#' progress curves of different data sets. In a nutshell, the method can be
#' divided into several steps. First, a reference data set must be chosen and a
#' model is fitted to it. Then for each other data set, the points are moved
#' along the x-axis (time axis), under some constraints, until the fitted model
#' match as well as possible the reference model. The relative movements of the
#' points are totally captured by two parameters, a and b, which correspond to
#' the translation and "streching" parameters. Use \code{sigmoflex_helper} for a
#' visual interpretation of a and b.
#'
#' @param data A data frame with 3 columns in the following order: the labels of
#' the data sets, the dates of the records, and the observation values (only on
#' a 0-1 scale at the moment).
#' @param reference The label of the reference data set.
#' @param ... Additional arguments to be passed to other methods.
#' @param threads Number of threads to perform the computations.
#'
#' @examples
#' # An example data frame to play with:
#' head(dummy_data)
#' str(dummy_data)
#'
#' res <- sigmoflex(dummy_data, "641")
#' res
#' plot(res)
#'
#' summary_res <- summary(res)
#' summary_res
#' plot(summary_res)
#'
#' sigmoflex_helper()
#'
#' @export
#------------------------------------------------------------------------------#
sigmoflex <- function(data, reference, ..., threads = 1) {
    stopifnot(is.data.frame(data))
    stopifnot(ncol(data) == 3)
    colnames(data) <- c("label", "time", "obs")
    stopifnot(is.numeric(data$time) && is.numeric(data$obs))
    data$label <- as.character(data$label)

    # Let's compute a glm model for the reference
    ref_x   <- data[data[["label"]] == reference, ]
    ref_glm <- glm(obs ~ time, family = quasibinomial, data = ref_x)

    # Fonction implementing the proposed formula:
    modify_time <- function(time, a, b) {
        stopifnot(length(time) >= 2) # Needs at least 2 observations over time.
        tlags   <- time[-1] - time[-length(time)]
        tlags   <- tlags * b
        time[1] <- time[1] - a
        sapply(2:length(time), function(i) time[i] <<- time[i - 1] + tlags[i - 1])
        time
    }

    # The objective is to minimize the sum of squared differences between the
    # ordinates of predicted values of the reference and the considered data set
    # (that we can modify playing with a and b). Predicted values are based on
    # general linear regression.
    # Here is the function to minimize:
    fun_to_mini <- function(time, obs, ref_model, a, b) {
        time     <- modify_time(time, a, b)
        new_glm  <- glm(obs ~ time, family = quasibinomial)
        new_pred <- predict(new_glm, newdata = data.frame(time = time))
        ref_pred <- predict(ref_model, newdata = data.frame(time = time))
        sum((new_pred - ref_pred)^2)
    }

    # Maximum likelihood estimation:
    epsilon <- 0.001 # To avoid perfect match issue for the reference.
    labels <- unique(data[["label"]])
    names(labels) <- labels # Trick to get a named list as lapply output.
    res_mle <- pbapply::pblapply(labels, function(lab) {
        dataset <- dplyr::filter(data, label == lab)
        res <- NULL
        # mle2(): wrapper over optim(), but other optimizer are possible.
        try(res <- bbmle::mle2(fun_to_mini,
                        data = list(time      = dataset[["time"]],
                                    obs       = dataset[["obs"]],
                                    ref_model = ref_glm),
                        start = list(a = 0 + epsilon, b = 1 + epsilon),
                        # a = 0 and b = 1 are values for the reference.
                        method = "L-BFGS-B"),
            silent = TRUE)
        res
    }, cl = threads)

    # Remove unsuccessful MLEs (i.e. return is null) if any:
    res_mle <- res_mle[!sapply(res_mle, is.null)]

    # Complete data:
    data$moved_time <- unlist(sapply(unique(data$label), function(name) {
        sub_data <- dplyr::filter(data, label == name)
        if (!is.null(res_mle[[name]])) { # A missing list element returns null.
            coefs <- bbmle::coef(res_mle[[name]])
            modify_time(sub_data$time, coefs[["a"]], coefs[["b"]])
        } else {
            rep(NA, nrow(sub_data))
        }
    }))
    data$is_ref <- data$label == reference

    # Theoretical data sets for curves:
    theo <- dplyr::bind_rows(lapply(unique(data$label), function(name) {
        time <- seq(min(c(data$time, data$moved_time)),
                    max(c(data$time, data$moved_time)),
                    by = 0.01)
        sub_data <- dplyr::filter(data, label == name)
        if (!is.null(res_mle[[name]])) { # A missing list element returns null.
            glm_premv  <- glm(obs ~ time, family = quasibinomial,
                              data = sub_data)
            glm_postmv <- glm(obs ~ moved_time, family = quasibinomial,
                              data = sub_data)
            obs_premv  <- predict(glm_premv,
                                  newdata = data.frame(time = time),
                                  type = "response")
            obs_postmv <- predict(glm_postmv,
                                  newdata = data.frame(moved_time = time),
                                  type = "response")
        } else {
            obs_premv <- obs_postmv <- rep(NA, nrow(sub_data))
            # Thus, we should keep the same colors for the points and the models
            # if there are unsuccessful MLEs. To check!
        }
        data.frame(label = name, time = time,
                   obs_premv = obs_premv, obs_postmv = obs_postmv,
                   stringsAsFactors = FALSE)
    }))

    # Returns:
    structure(list(mle  = res_mle,
                   data = data,
                   theo = theo),
              class = "sigmoflex")
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
print.sigmoflex <- function(x, ...) {
    cat("A sigmoflex object.\nUse summary() for more information.\n")
    invisible(x)
}

#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
plot.sigmoflex <- function(x, type = c("all", "before_moves", "after_moves"),
                           lines = FALSE, models = TRUE, ...) {
    type <- match.arg(type)
    if (type == "all" || type == "before_moves") {
        gg <- ggplot(x$data, aes(time, obs, col = label)) +
            geom_point()
        if (lines) {
            gg <- gg + geom_line(linetype = "dashed")
        }
        if (models) {
            gg <- gg + geom_line(data = x$theo, aes(time, obs_premv))
        }
        gg <- gg + scale_color_discrete("Data set") +
            labs(x = "Time", y = "Intensity") +
            geom_hline(yintercept = c(0, 1), linetype = "dashed",
                       col = "grey") +
            annotate("text",
                     x = min(c(x$data$time, x$data$moved_time)),
                     y = max(x$data$obs), label = "Before moves",
                     size = 8, hjust = "left", vjust = "top") +
            theme_bw()
        print(gg)
    }
    if (type == "all" || type == "after_moves") {
        gg <- ggplot(x$data, aes(moved_time, obs, col = label)) +
            geom_point()
        if (lines) {
            gg <- gg + geom_line(linetype = "dashed")
        }
        if (models) {
            gg <- gg + geom_line(data = x$theo, aes(time, obs_postmv))
        }
        gg <- gg + scale_color_discrete("Data set") +
            labs(x = "Time", y = "Intensity") +
            geom_hline(yintercept = c(0, 1), linetype = "dashed",
                       col = "grey") +
            annotate("text",
                     x = min(c(x$data$time, x$data$moved_time)),
                     y = max(x$data$obs), label = "After moves",
                     size = 8, hjust = "left", vjust = "top") +
            theme_bw()
        print(gg)
    }
    invisible(NULL)
}


#------------------------------------------------------------------------------#
#' @export
#------------------------------------------------------------------------------#
summary.sigmoflex <- function(object, ...) {
    # Let's regroup all the estimates of a and b in a data frame. It's easier to
    # handle.
    params <- do.call(rbind, lapply(seq_len(length(object$mle)), function(i) {
        coefs <- bbmle::coef(bbmle::summary(object$mle[[i]]))
        data.frame(label = names(object$mle)[i],
                   a = coefs["a", "Estimate"],
                   ase = coefs["a", "Std. Error"],
                   b = coefs["b", "Estimate"],
                   bse = coefs["b", "Std. Error"],
                   stringsAsFactors = FALSE)
    }))
    structure(list(params = params), class = "summary.sigmoflex")
}

#------------------------------------------------------------------------------#
#' @method print summary.sigmoflex
#' @export
#------------------------------------------------------------------------------#
print.summary.sigmoflex <- function(x, ...) {
    print(x$params)
    invisible(x)
}

#------------------------------------------------------------------------------#
#' @method plot summary.sigmoflex
#' @export
#------------------------------------------------------------------------------#
plot.summary.sigmoflex <- function(x, ...) {
    rect <- data.frame(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0)
    gg <- ggplot(x$params, aes(x = a, y = b, col = label)) +
        geom_rect(inherit.aes = FALSE, data = rect,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  color = "grey20", fill = "lightgrey", linetype = "dashed") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey20") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
        scale_colour_discrete(name = "Data set") +
        geom_point(size = 3) +
        geom_errorbarh(aes(xmin = (a - ase), xmax = (a + ase)), height = 0.1) +
        geom_errorbar( aes(ymin = (b - bse), ymax = (b + bse)), width = 0.1) +
        theme_bw()
    print(gg)
    invisible(NULL)
}

#------------------------------------------------------------------------------#
#' @rdname sigmoflex
#' @export
#------------------------------------------------------------------------------#
sigmoflex_helper <- function() {
    rect <- data.frame(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0)
    dummy <- data.frame(x=0,y=1)
    gg <- ggplot(data = dummy, aes(x, y)) +
        geom_rect(inherit.aes = FALSE, data = rect,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  color = "grey20", fill = "lightgrey", linetype = "dashed") +
        geom_hline(yintercept = 1, linetype = "dashed", color = "grey20") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
        geom_point(size = 5) +
        geom_segment(data = data.frame(x1 = c(0, 0, 1.5, 1.5),
                                       x2 = c(-0.6, 0.6, 1.5, 1.5),
                                       y1 = c(0.5, 0.5, 1.0, 1.0),
                                       y2 = c(0.5, 0.5, 1.4, 0.6)),
                     aes(x = x1, y = y1, xend = x2, yend = y2),
                     arrow = arrow(length = unit(0.4, "cm"))) +
        annotate(geom = "text",
                 x = c(0.30, 0.5, 1.2, 1.2, -0.3, 0.3),
                 y = c(1.15, -0.2, 1.2, 0.8, 0.5, 0.5),
                 label = c("Reference", "Impossible values",
                           "Faster", "Slower",
                           "Starts\nbefore", "Starts\nafter")) +
        xlim(-1, 2) +
        ylim(-0.5, 2) +
        labs(x = "a", y = "b") +
        theme_bw()
    print(gg)
}


