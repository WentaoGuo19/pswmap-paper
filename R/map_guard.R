#' @title PSW-MAP prior with a time limit and retries
#'
#' @description Runs pswmap::psw_map_prior() in a forked child process with a time limit. In a very small share of
#'    fits the Stan sampler stalls, typically when one source's weight rests on one or two subjects, so that its
#'    design-based variance is close to 0. Because R cannot interrupt compiled Stan code, the fit runs in a child
#'    process. If it has not finished within `timeout` seconds, the child and its chain processes are stopped and the
#'    fit is retried with new random numbers, up to `max_tries` times. The seed of each attempt is drawn from the
#'    calling process's random number stream, so a replicate stays reproducible. This guard keeps a long simulation
#'    running. It is not part of the method.
#'
#'    Forking is not available on Windows. There the prior is fitted directly, without a time limit.
#'
#' @param summaries,sigma_ref,tau_c,control,robust_weight passed to pswmap::psw_map_prior()
#' @param timeout seconds allowed for one attempt (a normal fit takes 1-2 seconds)
#' @param max_tries number of attempts
#'
#' @return the output of pswmap::psw_map_prior(), or NULL if every attempt failed or timed out
map_prior_guarded = function(summaries, sigma_ref, tau_c, control, robust_weight = 0.2, timeout = 60, max_tries = 3)
{
  fit = function() pswmap::psw_map_prior(summaries, sigma_ref = sigma_ref, tau_c = tau_c, control = control,
                                         robust_weight = robust_weight)
  if (.Platform$OS.type == "windows") return(fit())

  for (attempt in seq_len(max_tries)) {
    seed = sample.int(.Machine$integer.max, 1)
    job = parallel::mcparallel({ set.seed(seed); fit() }, mc.set.seed = FALSE, silent = TRUE)
    res = parallel::mccollect(job, wait = FALSE, timeout = timeout)
    if (!is.null(res) && !inherits(res[[1]], "try-error")) return(res[[1]])
    #time-out or error: stop the child's chain processes and the child itself, then collect it
    system(paste("pkill -9 -P", job$pid), ignore.stdout = TRUE, ignore.stderr = TRUE)
    tools::pskill(job$pid, tools::SIGKILL)
    invisible(suppressWarnings(parallel::mccollect(job, wait = TRUE)))
    message("map_prior_guarded: attempt ", attempt, " of ", max_tries,
            if (is.null(res)) paste0(" timed out after ", timeout, " s") else
              paste0(" failed: ", conditionMessage(attr(res[[1]], "condition"))))
  }
  NULL
}
