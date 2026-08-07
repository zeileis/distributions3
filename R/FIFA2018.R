#' Goals scored in all 2018 FIFA World Cup matches
#'
#' Data from all 64 matches in the 2018 FIFA World Cup along with predicted
#' ability differences based on bookmakers odds.
#'
#' To investigate the number of goals scored per match in the 2018 FIFA World Cup,
#' \code{FIFA2018} provides two rows, one for each team, for each of the matches
#' during the tournament. In addition some basic meta-information for the matches
#' (an ID, team name abbreviations, type of match, group vs. knockout stage),
#' information on the estimated log-ability for each team is provided. These
#' have been estimated by Zeileis et al. (2018) prior to the start of the
#' tournament (2018-05-20) based on quoted odds from 26 online bookmakers using
#' the bookmaker consensus model of Leitner et al. (2010). The difference in
#' log-ability between a team and its opponent is a useful predictor for the
#' number of goals scored.
#'
#' To model the data a basic Poisson regression model provides a good fit.
#' This treats the number of goals by the two teams as independent given the
#' ability difference which is a reasonable assumption in this data set.
#'
#' @usage data("FIFA2018", package = "distributions3")
#'
#' @format A data frame with 128 rows and 7 columns.
#' \describe{
#'   \item{goals}{integer. Number of goals scored in normal time (90 minutes), \
#'     i.e., excluding potential extra time or penalties in knockout matches.}
#'   \item{team}{character. 3-letter FIFA code for the team.}
#'   \item{match}{integer. Match ID ranging from 1 (opening match) to 64 (final).}
#'   \item{type}{factor. Type of match for groups A to H, round of 16 (R16), quarter final,
#'     semi-final, match for 3rd place, and final.}
#'   \item{stage}{factor. Group vs. knockout tournament stage.}
#'   \item{logability}{numeric. Estimated log-ability for each team based on
#'     bookmaker consensus model.}
#'   \item{difference}{numeric. Difference in estimated log-abilities between
#'     a team and its opponent in each match.}
#'   }
#'
#' @source The goals for each match have been obtained from Wikipedia
#'   (\url{https://en.wikipedia.org/wiki/2018_FIFA_World_Cup}) and the log-abilities
#'   from Zeileis et al. (2018) based on quoted odds from Oddschecker.com and Bwin.com.
#'
#' @references Leitner C, Zeileis A, Hornik K (2010).
#'   Forecasting Sports Tournaments by Ratings of (Prob)abilities: A Comparison for the EURO 2008.
#'   \emph{International Journal of Forecasting}, \bold{26}(3), 471-481.
#'   \doi{10.1016/j.ijforecast.2009.10.001}
#'
#' Zeileis A, Leitner C, Hornik K (2018).
#'   Probabilistic Forecasts for the 2018 FIFA World Cup Based on the Bookmaker Consensus Model.
#'   Working Paper 2018-09, Working Papers in Economics and Statistics,
#'   Research Platform Empirical and Experimental Economics, University of Innsbruck.
#'   \url{https://econpapers.repec.org/RePEc:inn:wpaper:2018-09}
#'
#' @examples
#' ## load data
#' data("FIFA2018", package = "distributions3")
#'
#' ## observed relative frequencies of goals in all matches
#' obsrvd <- prop.table(table(FIFA2018$goals))
#'
#' ## expected probabilities assuming a simple Poisson model,
#' ## using the average number of goals across all teams/matches
#' ## as the point estimate for the mean (lambda) of the distribution
#' p_const <- Poisson(lambda = mean(FIFA2018$goals))
#' p_const
#' expctd <- pdf(p_const, 0:6)
#'
#' ## comparison: observed vs. expected frequencies
#' ## frequencies for 3 and 4 goals are slightly overfitted
#' ## while 5 and 6 goals are slightly underfitted
#' cbind("observed" = obsrvd, "expected" = expctd)
#'
#' ## instead of fitting the same average Poisson model to all
#' ## teams/matches, take ability differences into account
#' m <- glm(goals ~ difference, data = FIFA2018, family = poisson)
#' summary(m)
#' ## when the ratio of abilities increases by 1 percent, the
#' ## expected number of goals increases by around 0.4 percent
#'
#' ## this yields a different predicted Poisson distribution for
#' ## each team/match
#' p_reg <- Poisson(lambda = fitted(m))
#' head(p_reg)
#'
#' ## as an illustration, the following goal distributions
#' ## were expected for the final (that France won 4-2 against Croatia)
#' p_final <- tail(p_reg, 2)
#' p_final
#' pdf(p_final, 0:6)
#' ## clearly France was expected to score more goals than Croatia
#' ## but both teams scored more goals than expected, albeit not unlikely many
#'
#' ## assuming independence of the number of goals scored, obtain
#' ## table of possible match results (after normal time), along with
#' ## overall probabilities of win/draw/lose
#' res <- outer(pdf(p_final[1], 0:6), pdf(p_final[2], 0:6))
#' sum(res[lower.tri(res)]) ## France wins
#' sum(diag(res))           ## draw
#' sum(res[upper.tri(res)]) ## France loses
#'
#' ## update expected frequencies table based on regression model
#' expctd <- pdf(p_reg, 0:6)
#' head(expctd)
#' expctd <- colMeans(expctd)
#' cbind("observed" = obsrvd, "expected" = expctd)
"FIFA2018"
