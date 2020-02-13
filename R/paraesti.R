


paraesti <- function(rdd_dat_indh, j, y, df) {
  reg_para <- rdd_reg_lm(rdd_dat_indh, order = 1)

  if (y == 0) {
    rmin <- min(award[, "Patient"]) - 20
    rmax <- max(award[, "Patient"]) + 10
    if (j == 0) {
     pl<- plot(
        reg_para,
        xlab = c("Period"),
        ylab = c("Number of patients"),
        ylim = range(rmin, rmax),
        col = "black"
      )
     par(new = TRUE)
    } else{
      plot(
        reg_para,
        xlab = c("Period"),
        ylab = c("Number of patients"),
        ylim = range(rmin, rmax),
        col = "red"
      )
    }
 #   par(new = TRUE)
  } else{
    rmin <- min(award[, "Views"]) - 3000
    rmax <- max(award[, "Views"]) + 2000
    if (j == 0) {
      pl<- plot(
        reg_para,
        xlab = c("Period"),
        ylab = c("Number of views"),
        ylim = range(rmin, rmax)
      )
      par(new = TRUE)
    } else{
        plot(
        reg_para,
        xlab = c("Period"),
        ylab = c("Number of views"),
        ylim = range(rmin, rmax),
        col = "red"
      )
    }
    #use par and plot on the same graph but different axis.
  }
  legend(
    "bottomleft",
    legend = c("Before awards", "After awards"),
    pch = c(16, 16),
    col = c("black", "red"),
    horiz = T
  )#, title="Group"
  return(pl)
}
