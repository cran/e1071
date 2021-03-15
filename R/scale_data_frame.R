scale_data_frame <-
  function(x, center = TRUE, scale = TRUE)
  {
    if (is.numeric(x)) return (scale(x, center, scale))
    i <- sapply(x, is.numeric)
    if (ncol(x[, i, drop = FALSE])) {
      x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
      if(center || !is.logical(center))
        attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
      if(scale || !is.logical(scale))
        attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
    }
    x
  }

