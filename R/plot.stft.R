plot.stft <- function (Y, col = gray (63:0/63))
  {
    x <- Y$values
    image(x=1:dim(x)[1], y=1:dim(x)[2], z=x, col=col)
}
