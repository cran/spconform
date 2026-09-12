## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5,
  fig.align = "center",
  eval = requireNamespace("sp", quietly = TRUE)
)

## ----setup--------------------------------------------------------------------
library(spconform)

## -----------------------------------------------------------------------------
library(sp)
data(meuse)

s <- as.matrix(meuse[, c("x", "y")])
y <- log(meuse$zinc)

## ----fig-layout---------------------------------------------------------------
plot(meuse$x, meuse$y, col = rgb(0.2, 0.4, 0.8, 0.5), pch = 19,
     xlab = "X coordinate", ylab = "Y coordinate",
     main = "Meuse Sampling Locations")

## -----------------------------------------------------------------------------
pred_fun <- function(s_train, y_train, s_new) {
  fit <- lm(y_train ~ s_train[, 1] + s_train[, 2] +
              I(s_train[, 1]^2) + I(s_train[, 2]^2))
  cbind(1, s_new[, 1], s_new[, 2], s_new[, 1]^2, s_new[, 2]^2) %*% coef(fit)
}

## -----------------------------------------------------------------------------
set.seed(1)
n <- nrow(s)
idx <- sample(n, floor(0.7 * n))

s_train <- s[idx, ]; y_train <- y[idx]
s_test  <- s[-idx, ]; y_test  <- y[-idx]

out <- scp_geostatistical(s_train, y_train, s_test, pred_fun,
                          alpha = 0.1, seed = 1)
print(out)

## -----------------------------------------------------------------------------
coverage_report(out, y_test)

## ----fig-intervals------------------------------------------------------------
plot(out, y_true = y_test)

## ----fig-diagnostics, fig.width = 7, fig.height = 5.5-------------------------
diag <- diagnose(out, y_true = y_test, s_test = s_test, plot = TRUE)
print(diag)

## -----------------------------------------------------------------------------
set.seed(123)
coverages <- numeric(50)
widths    <- numeric(50)

for (i in 1:50) {
  idx_i <- sample(n, floor(0.7 * n))
  s_tr <- s[idx_i, ]; y_tr <- y[idx_i]
  s_te <- s[-idx_i, ]; y_te <- y[-idx_i]

  out_i <- scp_geostatistical(s_tr, y_tr, s_te, pred_fun,
                              alpha = 0.1, seed = i)

  rep_i <- coverage_report(out_i, y_te)
  coverages[i] <- rep_i$coverage
  widths[i]    <- rep_i$mean_width
}

mean(coverages)
sd(coverages)
mean(widths)

## ----fig-coverage-hist--------------------------------------------------------
hist(coverages, breaks = 15, col = "lightblue", border = "white",
     main = "Empirical Coverage Across 50 Random Splits",
     xlab = "Empirical Coverage", xlim = c(0.7, 1))
abline(v = 0.90, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = "Nominal target (0.90)",
       col = "red", lty = 2, bty = "n")

## ----fig-spatial-width--------------------------------------------------------
plot_df <- data.frame(
  x = s_test[, 1],
  y = s_test[, 2],
  width = out$upper - out$lower
)

plot(plot_df$x, plot_df$y,
     cex = plot_df$width, pch = 19,
     col = rgb(0.2, 0.4, 0.8, 0.5),
     xlab = "X coordinate", ylab = "Y coordinate",
     main = "Spatial Distribution of Interval Width")

## -----------------------------------------------------------------------------
xbreaks <- seq(min(meuse$x), max(meuse$x), length.out = 7)
ybreaks <- seq(min(meuse$y), max(meuse$y), length.out = 7)

meuse$cell_x  <- cut(meuse$x, xbreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_y  <- cut(meuse$y, ybreaks, include.lowest = TRUE, labels = FALSE)
meuse$cell_id <- (meuse$cell_y - 1) * 6 + meuse$cell_x

agg <- aggregate(log(zinc) ~ cell_id, data = meuse, FUN = mean)
names(agg) <- c("cell_id", "y")

cell_coords <- unique(meuse[, c("cell_id", "cell_x", "cell_y")])
agg <- merge(agg, cell_coords, by = "cell_id")
agg <- agg[order(agg$cell_id), ]

n_cells <- nrow(agg)
adj <- matrix(0, n_cells, n_cells)
for (i in 1:n_cells) {
  for (j in 1:n_cells) {
    if (i != j) {
      dx <- abs(agg$cell_x[i] - agg$cell_x[j])
      dy <- abs(agg$cell_y[i] - agg$cell_y[j])
      if (dx <= 1 && dy <= 1) adj[i, j] <- 1
    }
  }
}

## -----------------------------------------------------------------------------
out2 <- scp_areal(agg$y, adjacency = adj, alpha = 0.2)
print(out2)
summary(out2)
coverage_report(out2, agg$y)

## ----fig-areal-intervals------------------------------------------------------
plot(out2, y_true = agg$y)

## ----fig-width-comparison-----------------------------------------------------
boxplot(list(Geostatistical = out$upper - out$lower,
             Areal          = out2$upper - out2$lower),
        main = "Interval Width Comparison",
        ylab = "Interval Width",
        col  = c("lightblue", "lightgreen"))

