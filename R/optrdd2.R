

library("rdd")
library("rddtools")

library(readxl)
#award <- read_excel("~/Downloads/Rpkg/optrdd2/data/sum_avg_award_paper.xlsx",sheet = "df")
#head(award)
#award<-award[1:4152,]
award <- data_avg_treat
award1 <- award[award$D == 0, ]
award2 <- award[award$D == 1, ]
#2 choose 1
##
#y <- 1#'Patient'
#y <- 0 #'Views'

for (y in seq(0, 1)) {
  for (j in seq(0, 1)) {
    if (j == 0) {
      df <- award1
    } else{
      df <- award2
    }

    if (y == 0) {
      rdd_dat_indh <- rdd_data(
        y = Patient,
        x = period,
        data = df,
        cutpoint = 6
      )
    } else{
      rdd_dat_indh <- rdd_data(
        y = Views,
        x = period,
        data = df,
        cutpoint = 6
      )
    }

    str(rdd_dat_indh)
    summary(rdd_dat_indh)
    head(rdd_dat_indh)
    #award<-sum_avg_award_paper
    #plot(rdd_dat_indh)
    paraesti(rdd_dat_indh,j,y,df)
  }
}







for (y in seq(0, 1)) {
  for (j in seq(0, 1)) {
    if (j == 0) {
      df <- award1
    } else{
      df <- award2
    }

    if (y == 0) {
      rdd_dat_indh <- rdd_data(
        y = Patient,
        x = period,
        data = df,
        cutpoint = 6
      )
    } else{
      rdd_dat_indh <- rdd_data(
        y = Views,
        x = period,
        data = df,
        cutpoint = 6
      )
    }

    str(rdd_dat_indh)
    summary(rdd_dat_indh)
    head(rdd_dat_indh)
    #award<-sum_avg_award_paper
    #plot(rdd_dat_indh)
    nonparaesti(rdd_dat_indh,j,y,df)
  }
}
