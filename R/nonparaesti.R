

nonparaesti<-function(rdd_dat_indh,j,y,df){
#  for (j in seq(0, 1)) {
    bw_ik <- rdd_bw_ik(rdd_dat_indh)
    reg_nonpara <-
      rdd_reg_np(rdd_object = rdd_dat_indh, bw = bw_ik)

    if (y == 0) {
      rmin <- min(award[, "Patient"]) - 20
      rmax <- max(award[, "Patient"]) + 10
      if (df == award1) {
        plot(
          reg_nonpara,
          xlab = c("Period"),
          ylab = c("Number of patients"),
          ylim = range(rmin, rmax)
        )
      } else{
        plot(
          reg_nonpara,
          xlab = c("Period"),
          ylab = c("Number of patients"),
          ylim = range(rmin, rmax),
          col = "red"
        )
      }
      par(new = TRUE)
    } else{
      rmin <- min(award[, "Views"]) - 3000
      rmax <- max(award[, "Views"]) + 2000
      if (df == award1) {
        plot(
          reg_nonpara,
          xlab = c("Period"),
          ylab = c("Number of views"),
          ylim = range(rmin, rmax)
        )
      } else{
        plot(
          reg_nonpara,
          xlab = c("Period"),
          ylab = c("Number of views"),
          ylim = range(rmin, rmax),
          col = "red"
        )
      }
      #use par and plot on the same graph but different axis.
      par(new = TRUE)
    }
    legend(
      "bottomleft",
      legend = c("Before awards", "After awards"),
      pch = c(16, 16),
      col = c("black", "red"),
      horiz = T
    )#, title="Group"
#  }
}
