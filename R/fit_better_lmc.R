library(gstat)

#' Make matrix positive definite
posdef = function(X) {
  q = eigen(X)
  d = q$values
  d[d < 0] = 0
  q$vectors %*% diag(d, nrow = length(d)) %*% t(q$vectors)
}

fit_better_lmc = function (v, g, model, ...)
{
  if (!inherits(v, "gstatVariogram"))
    stop("v should be of class gstatVariogram")
  if (!inherits(g, "gstat"))
    stop("g should be of class gstat")
  if (!missing(model)) {
    if (!inherits(model, "variogramModel"))
      stop("model should be of class variogramModel")
  }
  n = names(g$data)
  for (i in 1:length(n)) {
    for (j in i:length(n)) {
      name = ifelse(i == j, n[i], cross.name(n[i], n[j]))
      x = v[v$id == name, ]
      if (nrow(x) == 0)
        stop(paste("gstatVariogram", name, "not present"))
      m = g$model[[name]]
      if (!missing(model))
        m = model
      g$model[[name]] = fit.variogram(x, m, fit.ranges = fit.ranges,
                                      ...)
    }
  }

  m = g$model[[n[1]]]
  for (k in 1:nrow(m)) {
    psill = matrix(NA, nrow = length(n), ncol = length(n))
    for (i in 1:length(n)) {
      for (j in i:length(n)) {
        name = ifelse(i == j, n[i], cross.name(n[i],
                                               n[j]))
        psill[i, j] = psill[j, i] = g$model[[name]][k,
                                                    "psill"]
      }
    }
    psill = posdef(psill)
    diag(psill) = diag(psill) * correct.diagonal
    for (i in 1:length(n)) {
      for (j in i:length(n)) {
        name = ifelse(i == j, n[i], cross.name(n[i],
                                               n[j]))
        g$model[[name]][k, "psill"] = psill[i, j]
      }
    }
  }
  g
}
