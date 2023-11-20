library(survivalROC)
library(riskRegression)
# load your dataset
your_data <- read.csv("sample_data.csv")

# load two Fine and Gray models
model_VTE <- readRDS("VTE Risk Assessment Model for Lymphoma.rds")
model_PEDVT <- readRDS("DVT Risk Assessment Model for Lymphoma.rds")

# get predicted VTE or PE/LE-DVT risk for each patient at month "n"
n <- 6 ## could adjust to any time you'd like to predict
your_data$predited_VTE_incidence<- predict(model_VTE, newdata =your_data, type="risk",times = n)
your_data$predited_PEDVT_incidence<- predict(model_PEDVT, newdata =your_data, type="risk",times = n)

# prepare the binary outcomes for Time dependent AUC at month "n"
your_data$VTE_indicator <- ifelse(your_data$status_VTE=="1",1,0)
your_data$PEDVT_indicator <- ifelse(your_data$status_PEDVT=="1",1,0)

# TD_AUC using fgr predicted incidence as marker
make_survival_roc <- function(resampled_data, time_col, status_col, marker_col){
  survivalROC(resampled_data[[time_col]], resampled_data[[status_col]],
              marker = resampled_data[[marker_col]], predict.time = n, method = "NNE",
              span = 0.0055)
}

# make_bootstrap_ci calls this function
make_bootstrap_auc <- function(data, indices, time_col, status_col, marker_col) {
  resampled_data <- data[indices,]
  out1 <- make_survival_roc(resampled_data, time_col, status_col, marker_col)
  out1$AUC
}

# Use when you want a TD-ROC value and 95% CI that is bootstrapped
make_bootstrap_ci <- function(data, time_col, status_col, marker_col,
                              R, round_num = 3){
  boot_auc <- boot::boot(data = data, statistic = make_bootstrap_auc,
                         time_col = time_col, status_col = status_col,
                         marker_col = marker_col, R = R)
  est <- round(boot_auc$t0, round_num)
  ci <- boot::boot.ci(boot_auc, conf = 0.95, type = c("perc"))
  lower <- round(ci$percent[4], round_num)
  upper <- round(ci$percent[5], round_num)
  data.frame(Time = time_col,
             Status = status_col,
             Marker = marker_col,
             R = R,
             Est = est,
             Lower = lower,
             Upper = upper,
             Str = glue::glue("{est} ({lower}-{upper})"))
}

R <-10 # could be 1000 if your data size is large enough
TD_AUC_VTE <- make_bootstrap_ci(your_data, "time_VTE_month","VTE_indicator","predited_VTE_incidence",R)
TD_AUC_PEDVT <- make_bootstrap_ci(your_data, "time_PEDVT_month","PEDVT_indicator","predited_PEDVT_incidence",R)
