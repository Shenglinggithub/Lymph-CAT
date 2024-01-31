library(survivalROC)
library(riskRegression)
# load your dataset
your_data <- read.csv("sample_data.csv")

# load two Fine and Gray models
model_VTE <- readRDS("VTE Risk Assessment Model for Lymphoma.rds")
model_PEDVT <- readRDS("PEDVT Risk Assessment Model for Lymphoma.rds")



# make sure "Lymphoma_histology" is factor type with c("0","1", "2","3") levels.
your_data$Lymphoma_histology <-  factor(your_data$Lymphoma_histology, levels = c("0","1", "2","3"))
# get predicted VTE or PE/LE-DVT risk for each patient at 6 month
your_data$predited_VTE_incidence<- predict(model_VTE, newdata =your_data, type="risk",times = 6)
your_data$predited_PEDVT_incidence<- predict(model_PEDVT, newdata =your_data, type="risk",times = 6)

# prepare the binary outcomes for Time dependent AUC
your_data$VTE_indicator <- ifelse(your_data$status_VTE=="1",1,0)
your_data$PEDVT_indicator <- ifelse(your_data$status_PEDVT=="1",1,0)
unique(your_data$VTE_indicator)

## Use when you want a TD-ROC value and 95% CI that is bootstrapped
make_survival_roc <- function(resampled_data, time_col, status_col, marker_col){
  survivalROC(resampled_data[[time_col]], resampled_data[[status_col]],
              marker = resampled_data[[marker_col]], predict.time = 6, method = "NNE",
              # span = 0.04*nobs^(-0.2)
              span = 0.0055)
}

## make_bootstrap_ci calls this function
make_bootstrap_auc <- function(data, indices, time_col, status_col, marker_col) {
  resampled_data <- data[indices,]
  out1 <- make_survival_roc(resampled_data, time_col, status_col, marker_col)
  out1$AUC
}

## Use when you want a TD-ROC value and 95% CI that is bootstrapped
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
str(your_data)
R <-10 # could be 1000 if your sample size is large enough
TD_AUC_VTE <- make_bootstrap_ci(your_data, "time_VTE_month","VTE_indicator","predited_VTE_incidence",R)
TD_AUC_PEDVT <- make_bootstrap_ci(your_data, "time_PEDVT_month","PEDVT_indicator","predited_PEDVT_incidence",R)
TD_AUC_VTE
TD_AUC_PEDVT
