library(ggplot2)
out <- commandArgs(TRUE)[1]
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
erfinv <- function(y) qnorm((y + 1) / 2) / sqrt(2)
sc <- list(
  exp = function(x) rep(1, length(x)),
  linear = function(x) 1 / abs(x),
  log = function(x) log(abs(x)) * abs(x),
  gamma = function(x) 1 / digamma(x),
  factorial = function(x) 1 / digamma(x + 1),
  logit = function(x, a, b) -1 * (x - a)^2 * (-1 + (b - a) / (x - a)) * log(abs(-1 + (b - a) / (x - a))) / (b - a),
  expit = function(x, a, b) exp(x) * (1 + exp(-x))^2 * (a + (b - a) / (1 + exp(-x))) / (b - a),
  probitInv = function(x, a, b) sqrt(2) * exp(0.5 * x^2) * sqrt(pi) * (a + 0.5 * (b - a) * (1 + erf(x / sqrt(2)))) / (b - a),
  probit = function(x, a, b) { e <- erfinv(-1 + 2 * (x - a) / (b - a)); sqrt(2) * (b - a) * e / sqrt(pi) / 2 * exp(e^2) }
)
Nf <- list(
  logit = function(x, a, b) (x - a) * (b - x) / (b - a),
  probit = function(x, a, b) (x - a) * (b - x) / (b - a),
  expit = function(x, a, b) (a + (b - a) / (1 + exp(-x))) / (b - a),
  probitInv = function(x, a, b) (a + 0.5 * (b - a) * (1 + erf(x / sqrt(2)))) / (b - a)
)
frac <- list(logit = c(4e-4, 40), probit = c(4e-4, 40), expit = c(1, 100), probitInv = c(1, 100))
mk <- function(kind, xr, a = NA, b = NA, logx = FALSE, band = NULL, sing = NA, label) {
  x <- if (logx) exp(seq(log(xr[1]), log(xr[2]), length.out = 400)) else seq(xr[1], xr[2], length.out = 400)
  s <- if (is.na(a)) sc[[kind]](x) else sc[[kind]](x, a, b)
  if (!is.null(Nf[[kind]])) {
    n <- Nf[[kind]](x, a, b); lo <- n * frac[[kind]][1]; hi <- n * frac[[kind]][2]
  } else if (!is.null(band)) { lo <- band[1]; hi <- band[2] } else { lo <- NA; hi <- NA }
  inb <- is.finite(s) & s > 0 & (is.na(lo) | (s >= lo & s <= hi))
  data.frame(panel = label, x = x, y = pmin(pmax(abs(s), 1e-5), 1e5), lo = lo, hi = hi,
             region = ifelse(inb, "in band: formula kept", "out of band / <= 0 / wrong sign: fallback"),
             sing = sing, logx = logx)
}
plotSet <- function(d, file, title, ncol = 2, logx = FALSE) {
  d$panel <- factor(d$panel, levels = unique(d$panel))
  d$grp <- cumsum(c(1, diff(as.integer(factor(d$region))) != 0 | diff(as.integer(d$panel)) != 0))
  s <- unique(d[!is.na(d$sing), c("panel", "sing")])
  p <- ggplot(d, aes(x, y)) +
    geom_ribbon(aes(ymin = pmax(lo, 1e-5), ymax = pmin(hi, 1e5)), fill = "#2f9e6f", alpha = 0.15, na.rm = TRUE) +
    geom_vline(data = s, aes(xintercept = sing), linetype = 3, colour = "grey40") +
    geom_line(aes(colour = region, group = grp), linewidth = 1) +
    scale_y_log10() +
    scale_colour_manual(values = c("#2f9e6f", "#d9480f"), name = NULL) +
    facet_wrap(~panel, scales = "free", ncol = ncol) +
    labs(x = "initial estimate (theta)", y = "|scaleC| (log scale)", title = title,
         subtitle = "shaded: guard band; dotted: singularity of the transform's derivative-based scaleC") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  if (logx) p <- p + scale_x_log10()
  ggsave(file, p, width = 9, height = 3 * ceiling(length(levels(d$panel)) / ncol) + 1.2, dpi = 110)
}
bounded <- rbind(
  mk("logit", c(0.02, 0.98), 0, 1, sing = 0.5, label = "logit(theta, 0, 1)"),
  mk("logit", c(2, 99), 1, 100, sing = 50.5, label = "logit(theta, 1, 100)"),
  mk("expit", c(-6, 6), 0, 1, label = "expit(theta, 0, 1)"),
  mk("expit", c(-6, 6), 1000, 1001, label = "expit(theta, 1000, 1001)"),
  mk("probit", c(0.03, 0.97), 0, 1, sing = 0.5, label = "probit(theta, 0, 1)"),
  mk("probitInv", c(-2.6, 2.6), 0, 1, label = "probitInv(theta, 0, 1)")
)
plotSet(bounded, file.path(out, "scalec-bounded.png"), "Bounded transforms: band N * [lo, hi] built from each parameter's own bounds")
fixed <- rbind(
  mk("exp", c(0.01, 100), logx = TRUE, label = "exp() / log-normal: scaleC = 1"),
  mk("linear", c(1e-3, 1e3), logx = TRUE, band = c(0.1, 10), label = "linear: 1/|init|, band scaleCband = [0.1, 10]"),
  mk("log", c(0.1, 120), logx = TRUE, band = c(0.1, 10), sing = 1, label = "log(): log|init| * |init|"),
  mk("gamma", c(0.2, 12), logx = TRUE, band = c(0.1, 10), sing = 1.4616, label = "gamma: 1/digamma(init)"),
  mk("factorial", c(0.05, 8), logx = TRUE, band = c(0.1, 10), sing = 0.4616, label = "factorial: 1/digamma(init + 1)")
)
fixed$x2 <- fixed$x
plotSet(fixed, file.path(out, "scalec-fixed.png"), "Unbounded / structural transforms: fixed bands", logx = TRUE)
# linear: the scaleC actually used and bobyqa's first step (rhobeg = 0.2) as % of init
x <- exp(seq(log(1e-2), log(1e2), length.out = 600))
inv <- ifelse(1 / x >= 0.1 & 1 / x <= 10, 1 / x, x)
lin <- rbind(
  data.frame(x = x, step = 100 * 0.2 * inv / x, rule = 'scaleLinear = "inverse" (default)'),
  data.frame(x = x, step = 100 * 0.2 * x / x , rule = 'scaleLinear = "abs" (|init|, NONMEM)')
)
lin$grp <- paste(lin$rule, cumsum(c(1, abs(diff(log10(lin$step))) > 0.5)))
p <- ggplot(lin, aes(x, step, colour = rule, group = grp)) +
  annotate("rect", xmin = 0.1, xmax = 10, ymin = 0, ymax = Inf, fill = "#2f9e6f", alpha = 0.12) +
  geom_line(aes(linetype = rule), linewidth = 1) + scale_x_log10() + scale_y_log10() + scale_linetype_manual(values = c("22", "solid"), name = NULL) +
  scale_colour_manual(values = c("#0b7285", "#d9480f"), name = NULL) +
  labs(x = "initial estimate of a linear theta", y = "first bobyqa step, % of init",
       title = "Linear theta: first optimizer step relative to its value",
       subtitle = "shaded: scaleCband = [0.1, 10], where the default keeps 1/|init|") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(out, "scalec-linear-step.png"), p, width = 8, height = 4.5, dpi = 110)
