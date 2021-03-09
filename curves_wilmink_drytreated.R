

# Maj Beldring Henningsen, majbh@sund.ku.dk

# Curves: Drytreated
## drytreated only organic
## drytreated all farms
## ONLY HOLSTEIN
## PARITY DIVIDED (parity 2 for poster)

# exclude:
## Other AB treat

# create organic curves and all curves


#-------------------------------------------------------
# Packages and settings:

library(tidyverse)
library(gridExtra)
library(data.table)
library(plotly)
library(GGally)
library(tidymodels)
library(nlstools) # for bootstrapping
library('minpack.lm') # nlsLM (not using this for now)
library('nls2') # nls2 (not using this for now)
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
# Loading data:

load("M:/PCR_data/PCR_merge2.RData") 
rm(df_pcr, df_curve); gc() 

#-----------------------------------------------------------
# pre-preparing data;remove dates and DYR_ID:

df <- df_model %>% 
  dplyr::select(BES_ID, DYR_ID, PARITY, BREED, HERDTYPE, DIM, logSCC, 
                MILK, IMI, DRY_TREAT, PCR_TEST, RES_MAJOR, OTHER_AB, TEAT_TREAT)

# create coloumn counting herd occurences
df <- df %>%
  group_by(BES_ID, PARITY) %>%
  mutate(count = n())



#--------------------------------------------------------------
# 1: ECO CURVES = df1: choosing only organic farms and with no other treatments than dryoff

df1 <- df %>% 
  filter(HERDTYPE == 0) %>%
  filter(OTHER_AB == 0) %>%
  filter(BREED == 1) %>%
  filter(DIM < 306) %>%
  filter(count > 200) %>%
  dplyr::select(BES_ID, PARITY, BREED, DIM, logSCC, 
                MILK, IMI, DRY_TREAT)
df1$BES_ID <- factor(df1$BES_ID) # keep only used levels by resetting the variable


## TREATED, ORGANIC, NO OTHER TREATMENTS
# parity 2
df2 <- df1 %>%
  filter(PARITY == 2) %>% 
  filter(DRY_TREAT == 1)
# parity 3:
df3 <- df1 %>%
  filter(PARITY == 3) %>% 
  filter(DRY_TREAT == 1)
# parity 3+:
df4 <- df1 %>%
  filter(PARITY == 4) %>% 
  filter(DRY_TREAT == 1)

## NOT TREATED, ORGANIC, NO OTHER TREATMENTS
# parity 2:
df2_0 <- df1 %>%
  filter(PARITY == 2) %>% 
  filter(DRY_TREAT == 0)
# parity 2:
df3_0 <- df1 %>%
  filter(PARITY == 3) %>% 
  filter(DRY_TREAT == 0)
# parity 3+:
df4_0 <- df1 %>%
  filter(PARITY == 4) %>% 
  filter(DRY_TREAT == 0)



## STARTLIST for NLS, from equation 8 from Græsbøll
st <- list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6)
f_nls <- logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d | BES_ID


## RUN NLS for: Eco, treated, Holstein, no other AB, Parity stratified:
m_2 <- nlsList(f_nls, df2, start = sapply(st, mean),
               control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024,  printEval = FALSE, warnOnly = TRUE))
out_2 <- coef(m_2) %>%
  drop_na()

m_3 <- nlsList(f_nls, df3, start = sapply(st, mean),
               control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = TRUE))
out_3 <- coef(m_3) %>%
  drop_na()

m_4 <- nlsList(f_nls, df4, start = sapply(st, mean),
               control = list(maxiter = 400, tol = 1e-05, minFactor = 1/1024, printEval = FALSE, warnOnly = TRUE))
out_4 <- coef(m_4) %>%
  drop_na()


## Histograms of output. Filter must be applied to avoid crazy outliers

out_2 <- out_2 %>%
  drop_na() %>%
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5)

m2_a <- qplot(out_2$a,
              geom="histogram",
              binwidth = 0.1,  
              #main = "a: con, parity 2, holstein", 
              xlab = "Parity 2, a values",  
              fill=I("white"), col=I("green4"))
m2_b <- qplot(out_2$b,
              geom="histogram",
              binwidth = 0.0001,  
              #main = "b: con, parity 2, holstein", 
              xlab = "Parity 2, b values",  
              fill=I("white"), col=I("dodgerblue4"))
m2_k <- qplot(out_2$k,
              geom="histogram",
              binwidth = 0.1,  
              #main = "k: con, parity 2, holstein", 
              xlab = "Parity 2, k values",  
              fill=I("white"), col=I("darkorange"))
m2_d <- qplot(out_2$d,
              geom="histogram",
              binwidth = 0.1,  
              #main = "d: con, parity 2, holstein", 
              xlab = "Parity 2, d values",  
              fill=I("white"), col=I("firebrick4"))

out_3 <- out_3 %>%
  drop_na() %>%
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5)

m3_a <- qplot(out_3$a,
              geom="histogram",
              binwidth = 0.1,  
              #main = "a: con, parity 2, holstein", 
              xlab = "Parity 3, a values",  
              fill=I("white"), col=I("green4"))
m3_b <- qplot(out_3$b,
              geom="histogram",
              binwidth = 0.0001,  
              #main = "b: con, parity 2, holstein", 
              xlab = "Parity 3, b values",  
              fill=I("white"), col=I("dodgerblue4"))
m3_k <- qplot(out_3$k,
              geom="histogram",
              binwidth = 0.1,  
              #main = "k: con, parity 2, holstein", 
              xlab = "Parity 3, k values",  
              fill=I("white"), col=I("darkorange"))
m3_d <- qplot(out_3$d,
              geom="histogram",
              binwidth = 0.1,  
              #main = "d: con, parity 2, holstein", 
              xlab = "Parity 3, d values",  
              fill=I("white"), col=I("firebrick4"))

out_4 <- out_4 %>%
  drop_na() %>%
  filter(a>1, a<6) %>%
  filter(b>0.001, b<0.005) %>%
  filter(k > (-3), k < (-1)) %>%
  filter(d>0.9, d<3.5)
  
m4_a <- qplot(out_4$a,
              geom="histogram",
              binwidth = 0.1,  
              #main = "a: con, parity 2, holstein", 
              xlab = "Parity >3, a values",  
              fill=I("white"), col=I("green4"))
m4_b <- qplot(out_4$b,
              geom="histogram",
              binwidth = 0.0001,  
              #main = "b: con, parity 2, holstein", 
              xlab = "Parity >3, b values",  
              fill=I("white"), col=I("dodgerblue4"))
m4_k <- qplot(out_4$k,
              geom="histogram",
              binwidth = 0.1,  
              #main = "k: con, parity 2, holstein", 
              xlab = "Parity >3, k values",  
              fill=I("white"), col=I("darkorange"))
m4_d <- qplot(out_4$d,
              geom="histogram",
              binwidth = 0.1,  
              #main = "d: con, parity 2, holstein", 
              xlab = "Parity >3, d values",  
              fill=I("white"), col=I("firebrick4"))



# Histogram of Parameters: plotting all parity 2 histograms:
grid.arrange(m2_a, m2_b, m2_k, m2_d, 
             m3_a, m3_b, m3_k, m3_d,
             m4_a, m4_b, m4_k, m4_d,
             ncol=4, nrow=3)


## BOOTstrapping:

df_boo2 <- df1 %>%
  filter(PARITY == 2) %>% 
  filter(DRY_TREAT == 0) %>%
  ungroup() %>%
  dplyr::select(BES_ID, DIM, logSCC)

# df_boo2 <- head(df_boo2, 1000)

nls_boo2 <- nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2, 
             start = list(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
             control = list(maxiter = 100, tol = 1e-05, minFactor = 1/1024, 
                            printEval = FALSE, warnOnly = TRUE), 
             algorithm = "port", trace = FALSE)


boo_2 <- nlsBoot(nls_boo2, niter = 400)


# Matrix with the bootstrapped parameter estimates
Theta_mat <- boo_2$coefboot

# Model
# logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d, data = df_boo2
fun <- function(DIM, theta) theta["a"] + theta["b"]*DIM + exp(-(exp(theta["k"]) * DIM)*theta["d"])

# Points where to evaluate the model
x_eval <- seq(min(df_boo2$DIM), max(df_boo2$DIM), length.out = 100)

# Matrix with the predictions
Pred_mat <- apply(Theta_mat, 1, function(theta) fun(x_eval, theta))

# Pack the estimates for plotting
Estims_plot <- cbind(
  x = x_eval, 
  as.data.frame(t(apply(Pred_mat, 1, function(y_est) c(
    median_est = median(y_est), 
    ci_lower_est = quantile(y_est, probs = 0.025, names = FALSE), 
    ci_upper_est = quantile(y_est, probs = 0.975, names = FALSE)
  ))))
)

p_boo <- ggplot(data = Estims_plot, aes(x = x, y = median_est, ymin = ci_lower_est, ymax = ci_upper_est)) + 
  geom_ribbon(alpha = 0.7, fill = "grey") + 
  geom_line(size = rel(1.5), colour = "black") + 
  #geom_point(data = df_boo2, aes(x = DIM, y = logSCC), size = rel(1), colour = "darkgreen", inherit.aes = FALSE) + 
  theme_bw() + labs(title = "Wilmink curve test 1 with Bootstrap", x = "DIM", y = "logSCC")
#ggsave("bootpstrap_results.pdf", height = 5, width = 9)

p_boo



#-------------------------------------------------------------------------
# ggpairs tidy model: intercorrelation of wilmink curve parametres
# Only working with test data (map function and singular gradient issues w. full data)
# 
# df_tidy <- df_nls %>%
#   filter(BREED == "1") %>%
#   filter (HERDTYPE == "1") %>%
#   dplyr::select(BES_ID, DIM, logSCC)
# 
# df_tidy <- tail(df_tidy, 5000)

df_tidy <- df_2

tidy_nls <- df_tidy %>%
  select(c("BES_ID", "DIM", "logSCC")) %>%  # select needed columns
  mutate(BES_ID = BES_ID %>% fct_drop()) %>%
  nest(yields = c(logSCC, DIM)) %>%
  mutate(model = map(
    yields, safely(~ nls(formula = logSCC ~ a + b * DIM + exp(-(exp(k)) * DIM)*d,
                  start = c(a = 3.9, b = 0.0027, k = -1.94, d = 2.6),
                  data = .x, 
                  control = list(
                    maxiter = 400, tol = 1e-05, minFactor = 1/1024, 
                    printEval = FALSE, warnOnly = TRUE),
                  algorithm = "port", 
                  lower=c(0.0001,0.000001,-10,0.0001),upper=c(20,0.7,-0.001,10),
                  trace = FALSE)))) # data must be vector, not dataframe





# model tidying + adjusting p values:
slopes_nls <- tidy_nls %>%
  mutate(coefs = map(model, tidy)) %>%
  unnest(coefs) 

tidy_nls
slopes_nls

pairs_nls <- slopes_nls %>% 
  select(BES_ID, term, estimate) %>% 
  pivot_wider(names_from = term, values_from = estimate, id_cols = BES_ID  )

ggpairs(pairs_nls, columns = c("a","b","k","d"))



#-----------------------------------------------------------------------------------
# Create wilmink curve with parameters:
# median, CI 75 and CI 95:

save.image("M:/PCR_data/PCR_nls.RData") # save data and parameters for PCR_curves


# END OF CODE..

