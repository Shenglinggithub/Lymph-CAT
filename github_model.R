model <- readRDS("VTE Risk Assessment Model for Lymphoma.rds")
#summary(model)
your_data <- read.csv("sample_data.csv")
your_data$predited_incidence<- predict(model, newdata =your_data, type="risk",times = 6)
