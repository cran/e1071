scale_data_frame <-
  function(x, center = TRUE, scale = TRUE)
  {
    if (isFALSE(center) && isFALSE(scale)) return(x)
    if (!is.data.frame(x)) return (scale(x, center, scale))
    i <- vapply(x, is.numeric, NA) | vapply(x, is.logical, NA)
    if (any(i)) {
      x[, i] <- tmp <- scale.default(x[, i, drop = FALSE], na.omit(center), na.omit(scale))
      if(center || !is.logical(center))
        attr(x, "scaled:center")[i] <- attr(tmp, "scaled:center")
      if(scale || !is.logical(scale))
        attr(x, "scaled:scale")[i]  <- attr(tmp, "scaled:scale")
    }
    x
  }

