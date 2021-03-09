

# Maj Beldring Henningsen, majbh@sund.ku.dk

# NLS Curves using Bootstrapping


#-------------------------------------------------------
# Packages and settings:

library(tidyverse)
library(gridExtra)
library(data.table)
library(plotly)
library(GGally)
library(tidymodels)
library(nlstools) # for bootstrapping
library(nlme) # for nlslist
library(nlshelper) # for tidy(fit)
#library(ggExtra)
#library(ggalluvial)
Sys.setlocale("LC_ALL","English") # data formatting
memory.size()            # Checking your memory size
memory.limit()           # Checking the set limit
memory.limit(size=56000) # suggest for 64 bit
options(stringsAsFactors = FALSE) # prevent factorizing caracters


#-------------------------------------------------------
# Loading data and preparing data:

load("M:/PCR_data/PCR_merge2.RData") 
rm(df_pcr, df_curve); gc() 

# pre-preparing data;remove dates and DYR_ID:
df <- df_model %>% 
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB, TEAT_TREAT)

# ECO, No other AB treatments, Holstein, Parity 2:
df <- df %>% 
  filter(OTHER_AB == 0) %>%
  filter(BREED == 1) %>%
  filter(PARITY == 2) %>%
  filter(HERDTYPE == 0) %>%
  filter(DIM < 306) %>%
  dplyr::select(BES_ID, HERDTYPE, PARITY, BREED, DIM, logSCC, 
                MILK, IMI, DRY_TREAT)

#----------------------------------------------------------------------
# ECO, drytreated vs non drytreated

### TREATED
# only treated with max 200 observations:
df_treat <- df %>%
  filter(DRY_TREAT == 1) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)
df_treat$BES_ID <- factor(df_treat$BES_ID) # keep only used levels by resetting the variable


# Boot fitting treated:
nls_treat <- nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_treat, 
             start = list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
             control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                            printEval = FALSE, warnOnly = TRUE), 
             algorithm = "port", trace = FALSE)

boo_treat <- nlsBoot(nls_treat, niter = 400)

# boxplot of parameters
plot(boo_treat, type = "boxplot", ask = FALSE)

# Matrix with the bootstrapped parameter estimates
Theta_treat <- boo_treat$coefboot

# Model: logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2
fun <- function(DIM, theta) theta["a"] + theta["b"]*DIM + exp(-(exp(theta["k"]) * DIM)*theta["d"])

# Points where to evaluate the model
x_eval <- seq(min(df_treat$DIM), max(df_treat$DIM), length.out = 100)

# Matrix with the predictions
Pred_treat <- apply(Theta_treat, 1, function(theta) fun(x_eval, theta))

# Pack the estimates for plotting
Estims_treat <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_treat, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)

p_treat <- ggplot(data = Estims_treat, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(0.5), colour = "darkgreen") + 
  ylim(3.95, 5.0) +
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + labs(title = "Wilmink for: Treated, parity 2, organic, Holstein", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_treat




#--------------------------------------------------------------------------
### NON TREATED

# only non treated with max 200 observations:
df_no <- df %>%
  filter(DRY_TREAT == 0) %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n()) %>%
  filter(count > 200) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)

df_no$BES_ID <- factor(df_no$BES_ID) # keep only used levels by resetting the variable

# Boot fitting NON treated:
nls_no <- nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_no, 
                 start = list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
                 control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                                printEval = FALSE, warnOnly = TRUE), 
                 algorithm = "port", trace = FALSE)

boo_no <- nlsBoot(nls_no, niter = 400)

# boxplot of parameters
plot(boo_no, type = "boxplot", ask = FALSE)

# Matrix with the bootstrapped parameter estimates
Theta_no <- boo_no$coefboot

# Model: logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2
fun <- function(DIM, theta) theta["a"] + theta["b"]*DIM + exp(-(exp(theta["k"]) * DIM)*theta["d"])

# Points where to evaluate the model: same as for treated so not changing this
x_eval <- seq(min(df_treat$DIM), max(df_treat$DIM), length.out = 100)

# Matrix with the predictions
Pred_no <- apply(Theta_no, 1, function(theta) fun(x_eval, theta))

# Pack the estimates for plotting
Estims_no <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_no, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)

p_no <- ggplot(data = Estims_no, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(0.5), colour = "darkred") + 
  ylim(3.95, 5.0) +
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + 
  labs(title = "Wilmink for: NON Treated, parity 2, organic, Holstein", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)
p_no




