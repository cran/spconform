## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6.5,
  fig.height = 5,
  fig.align = "center",
  eval = requireNamespace("sp", quietly = TRUE)
)

## ----setup--------------------------------------------------------------------
library(spconform)

## -----------------------------------------------------------------------------
library(sp)
data(meuse)

coords <- as.matrix(meuse[, c("x", "y")])
coords_scaled <- scale(coords)
y <- log(meuse$zinc)

## ----fig-layout---------------------------------------------------------------
plot(meuse$x, meuse$y, col = rgb(0.2, 0.4, 0.8, 0.6), pch = 19,
     xlab = "Easting (m)", ylab = "Northing (m)",
     main = "Meuse River Topsoil Sampling Locations (n = 155)")

## -----------------------------------------------------------------------------
pred_fun <- function(s_train, y_train, s_new) {
  df_tr <- data.frame(y = y_train, x1 = s_train[, 1], x2 = s_train[, 2])
  df_new <- data.frame(x1 = s_new[, 1], x2 = s_new[, 2])
  fit <- lm(y ~ x1 + x2 + I(x1^2) + I(x2^2) + I(x1 * x2), data = df_tr)
  as.numeric(predict(fit, newdata = df_new))
}

## -----------------------------------------------------------------------------
set.seed(42)
n <- nrow(coords_scaled)
train_idx <- sample(n, floor(0.70 * n))
test_idx  <- setdiff(seq_len(n), train_idx)

s_train <- coords_scaled[train_idx, ]; y_train <- y[train_idx]
s_test  <- coords_scaled[test_idx, ];  y_test  <- y[test_idx]

out <- scp_geostatistical(
  s_train = s_train,
  y_train = y_train,
  s0 = s_test,
  pred_fun = pred_fun,
  alpha = 0.10,
  split = 0.50,
  seed = 123
)

print(out)
summary(out)

# Standard S3 methods for seamless integration with R workflows:
head(predict(out, interval = "prediction"))
head(residuals(out, y_true = y_test, type = "abs"))
head(as.data.frame(out))

## ----fig-intervals------------------------------------------------------------
plot(out, y_true = y_test)

## ----fig-diagnostics, fig.width = 7.5, fig.height = 6-------------------------
# Generate structured diagnostic object
diag_res <- diagnose(out, y_true = y_test, s_test = s_test, plot = FALSE)

# S3 print method displays text summary including Moran's I
print(diag_res)

# S3 plot method produces multi-panel spatial diagnostic layout
plot(diag_res)

## -----------------------------------------------------------------------------
# Aggregate Meuse observations onto a 6x6 spatial grid
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

# Construct Rook/Queen graph adjacency matrix
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

# Run 80% Areal Conformal Prediction
out_areal <- scp_areal(y = agg$y, adjacency = adj, alpha = 0.20, decay = 1.0)
print(out_areal)
coverage_report(out_areal, agg$y)

## ----fig-areal-intervals------------------------------------------------------
plot(out_areal, y_true = agg$y)

