

# This script is provided in relation to the data analysis described in Section 2.5 of the Methods section of the article:
# Title: The effect of the discontinued use of zinc oxide on antimicrobial usage in Danish pig farms
#
# The primary objective of the analysis, described in the article,
# was to quantify the effect of the nationwide ban on zinc oxide usage, 
# as well as the influence of individual farms and veterinarians, on farm-level AMU.
#
# Using pig-farm and prescription data from 2018 to 2023 from Danish national databases, 
# we fitted a linear mixed-effect model to a three-level nested dataset, consisting of 
# monthly average standardized AMU (Defined Animal Daily Doses per pig-day), on a farm (n=4,020), overseen by a veterinarian (n=146).
#
# The script is shared to facilitate transparency, ensure reproducibility, and make the methodology accessible.
#
# The data provided with this script has been anonymized to ensure compliance with the General Data Protection Regulation (GDPR).
# Specifically, all farm characteristics have been removed, and farm and vet ID numbers have been altered to protect confidentiality.
#
# While the results from this version may not match exactly with those in the article, 
# due to the removal of structural fixed effects from the linear mixed-effects model, 
# the script remains consistent with the methodology outlined in the article.
#
# Generative AI tools have been applied in the development of this script.
#
# For any questions, please contact Josefine Ostenfeld Nielsen (josos@dtu.dk).








################################################################################
################################################################################
################################################################################

# Reading packages and libraries 

local({
  
  # Function to check and install a package if not present
  check_and_install <- function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(package_name)
    }
    library(package_name, character.only = TRUE)
  }
  
  # List of required packages
  packages <- c(
    # LMM 
    "nlme", 
    "MuMIn",
    # Plot 
    "ggplot2", 
    "gridExtra", 
    # Data
    "dplyr", 
    "tidyr", 
    "forcats"
  )
  
  # Install (if necessary) each package
  lapply(packages, check_and_install)
  
})








################################################################################
################################################################################
################################################################################

### Loading and preparing data

# In session, set working directory to source file location to load data
# Reading data 
data <- read.table("data_restricted.txt", sep = "\t", header = TRUE)

# Preparing data
data <- data %>% 
  mutate_at(vars(
    sin_month,            # calendar month (sin curve)
    cos_month,            # calendar month (cos curve)
    study_month,          # study month (1-72)
    ADD,                  # number of ADD used (defined animal daily doses) 
    log_ADD               # log transformed ADD
  ), as.numeric) %>% 
  mutate_at(vars(
    farm_id,              # farm id (anonymized)
    vet_id,               # veterinarian id (anonymized)
    rearing_stage,        # rearing stage og pigs (weaners/finishers)
    zinc_status,          # zinc use in weaners before ban (non-using/using)
    zinc,                 # zinc temporal variable (0 = before, 1 = 1-5 months after, and 2 = >5 months after ban/discontinuation)
    zinc_using,           # zinc temporal variable only for zinc using weaners (else 0)
    zinc_non_using        # zinc temporal variable only for zinc non-using weaners (else 0)
  ), as.factor) %>% 
  
  # Defining reference levels of categorical explanatory variables:
  mutate(                 
    zinc_status = fct_relevel(zinc_status, "using"),
    zinc = fct_relevel(zinc, "0"),
    zinc_using = fct_relevel(zinc_using, "0"),
    zinc_non_using = fct_relevel(zinc_non_using, "0")
  )


# Splitting date into a weaner and finisher dataframe
df_weaner <- data %>% filter(rearing_stage=="weaner")
df_finisher <- data %>% filter(rearing_stage=="finisher") 








################################################################################
################################################################################
################################################################################

### Determining autocorrelation structure in model 0

# An autoregressive residual covariance structure was used to account for temporal dependencies between study months within farms. 
# Different autoregressive terms (AR 1-3) was included in Model 0 (no fixed effects, only the random veterinarian and farm effects)
# ACF plots of normalized residuals were visually inspected to determine the appropriate number of autoregressive terms.


# Defining data, replace to test autocorrelation for finisher data
df <- df_weaner  


local({
  
  # Model 0 with no autocorrelation structure 
  m_none <- lme(
    fixed = log_ADD ~ 1,
    random = list(vet_id = ~ 1, farm_id = ~ 1 | vet_id),
    method = "ML",
    data = df
  )
  
  # Model 0 with a residual covariance structure of autoregressive order 1 (AR1)
  m_ar1 <- lme(
    fixed = log_ADD ~ 1,
    random = list(vet_id = ~ 1, farm_id = ~ 1 | vet_id),
    correlation = corARMA(p = 1, q = 0, form = ~ study_month | vet_id/farm_id),
    method = "ML",
    data = df
  )
  
  # Model 0 with a residual covariance structure of autoregressive order 2 (AR2)
  m_ar2 <- lme(
    fixed = log_ADD ~ 1,
    random = list(vet_id = ~ 1, farm_id = ~ 1 | vet_id),
    correlation = corARMA(p = 2, q = 0, form = ~ study_month | vet_id/farm_id),
    method = "ML",
    data = df
  )
  
  # Model 0 with a residual covariance structure of autoregressive order 3 (AR3)
  m_ar3 <- lme(
    fixed = log_ADD ~ 1,
    random = list(vet_id = ~ 1, farm_id = ~ 1 | vet_id),
    correlation = corARMA(p = 3, q = 0, form = ~ study_month | vet_id/farm_id),
    method = "ML",
    data = df
  )
  
  # Saving models in list
  models <- list(
    "none" = m_none, 
    "AR(1)" = m_ar1, 
    "AR(2)" = m_ar2, 
    "AR(3)" = m_ar3
  )
  
  # Printing phi values for the autocorrelation structure of each model in the list
  for (model_name in names(models)) {
    model <- models[[model_name]]
    print(model_name)
    print(model[["modelStruct"]][["corStruct"]])
    cat("________________________________________________\n\n")
  }
  
  # Empty list to store ACF plots
  plots <- list()
  
  # Generate ACF plots for models in list
  for (i in seq_along(models)) {
    model_name <- names(models)[i]  
    acf_result <- ACF(models[[i]], maxLag = 12, resType = "normalized")
    plots[[i]] <- ggplot(acf_result, aes(x = lag, y = ACF)) + 
      geom_bar(stat = "identity", position = "dodge", fill = "black", width = 0.25) +
      geom_hline(yintercept = 0) + 
      ggtitle(paste("Autocorrelation structure:", model_name)) +
      labs(x = "", y = "") +
      coord_cartesian(ylim=c(-0.1,1), xlim=c(0,12))+
      scale_x_continuous(breaks = 0:12)+
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
      ) 
  }
  
  # Arranging plots in a 2-column grid
  grid_arrange_plots <- grid.arrange(grobs = plots, ncol = 2)
  
  # Adding common x and y axis labels
  grid.text("Lag [months]", y = unit(0.02, "npc"), gp = gpar(fontsize = 12))  
  grid.text("Autocorrelation", x = unit(0.02, "npc"), rot = 90, gp = gpar(fontsize = 12))  
})








################################################################################
################################################################################
################################################################################

### Estimating R2 and ICC in model 1

# Model 1 is an extension of model 0, including temporal fixed effects.
# This included:
    # the long-term time trend (continuous linear effect of study month 1-72),
    # seasonality (sine-cosine curve of the calendar months),
    # the categorical zinc oxide ban/discontinuation variable (before, after 1-5 months and after >5 months), 
    # and in the weaner model two interaction terms wtih the zinc oxide status on weaner farms, 
      # to estimate long-term time trend and the effect of the zinc ban/discontinuation seperatly for weaners on zinc-using and non-zinc-using farms.
# Structural fixed effects are not included in the publicly available script.
# R-squared and ICC are estimated before random slopes are added in the final model. 


# Defining data, replace to estimate R2 and ICC for finisher data
df <- df_weaner


# Execute all lines (ctrl+enter here)
local({
  if (identical(df, df_weaner)) {
    
    # Model 1 for weaners
    m1 <- lme(
      fixed = log_ADD ~ 
        study_month*zinc_status + 
        sin_month + 
        cos_month + 
        zinc*zinc_status,
      random = list(vet_id = pdDiag(~ 1), farm_id = pdDiag(~ 1)),
      correlation = corARMA(p = 1, q = 0, form = ~ study_month | vet_id/farm_id),
      method = "ML",
      data = df
    )
    r2 <- data.frame(r.squaredGLMM(m1))
  } else if (identical(df, df_finisher)) {
    
    # Model 1 for finishers
    m1 <- lme(
      fixed = log_ADD ~ 
        study_month + 
        sin_month + 
        cos_month + 
        zinc,
      random = list(vet_id = pdDiag(~ 1), farm_id = pdDiag(~ 1)),
      correlation = corARMA(p = 1, q = 0, form = ~ study_month | vet_id/farm_id),
      method = "ML",
      data = df
    )
    r2 <- data.frame(r.squaredGLMM(m1))
  }
  # model summary
  # round(summary(m1)$tTable, 3)
  
  
  # R2 (coefficients of determination) = variance explained by model (in total and amount explained by fixed and random effects)
  cat(
    "variance explained by model (total):", round(r2$R2c * 100, 1), "%\n",
    "variance explained by fixed effects:", round(r2$R2m * 100, 1), "%\n",
    "variance explained by random effects:", round((r2$R2c-r2$R2m) * 100, 1), "%\n\n"
  )
  
  # ICC (intra-class correlation coefficient) = variance explained by random effects, attributed to unknown farm-level or veterinarian level factors
  cat(
    "variance explained by random effects, attributed to unknown veterinarian-level factors:", 
    round(as.numeric(VarCorr(m1)[2, "Variance"]) / 
            (as.numeric(VarCorr(m1)[2, "Variance"]) +
               as.numeric(VarCorr(m1)[4, "Variance"])) * 100, 1), "%\n",
    
    "variance explained by random effects, attributed to unknown farm-level factors:", 
    round(as.numeric(VarCorr(m1)[4, "Variance"]) /
            (as.numeric(VarCorr(m1)[2, "Variance"]) +
               as.numeric(VarCorr(m1)[4, "Variance"])) * 100, 1), "%\n"
  )
})









################################################################################
################################################################################
################################################################################

### Estimating final log-linear mixed-effect model (model 3)
# Model 3 is an extension of model 1, including radom slopes of the effect the zinc ban/discontinuation on the individual farms.


# Model 3 for weaners 
m3_weaner <- lme(
  fixed = log_ADD ~ 
    study_month*zinc_status + 
    sin_month + 
    cos_month + 
    zinc*zinc_status,
  random = list(
    vet_id = pdDiag(~ 1), 
    farm_id = pdDiag(~ 1 + zinc_using + zinc_non_using)),
  correlation = corARMA(p = 1, q = 0, form = ~ study_month | vet_id/farm_id),
  method = "ML",
  data = df_weaner
)
# model summary
# round(summary(m3_weaner)$tTable, 3)

# Model 3 for finishers 
m3_finisher <- lme(
  fixed = log_ADD ~ 
    study_month + 
    sin_month + 
    cos_month + 
    zinc,
  random = list(
    vet_id = pdDiag(~ 1), 
    farm_id = pdDiag(~ 1 + zinc)),
  correlation = corARMA(p = 1, q = 0, form = ~ study_month | vet_id/farm_id),
  method = "ML",
  data = df_finisher
)
# model summary
# round(summary(m3_finisher)$tTable, 3)









################################################################################
################################################################################
################################################################################

### Transforming the derived log-scale fixed effects and 95% CI to obtain the percentage change in AMU (ADDkg/pig-day) 
 # of to a one unit increase in the continuous predictors or compared to the reference level of the categorical predictors

local({
  
  # Transfoming weaner model 3 output
  df <- as.data.frame(intervals(m3_weaner, which = "fixed")[["fixed"]]) %>% rename("x"="est.", "l"="lower", "u"="upper")
  df_m3_weaner <- df %>% 
    dplyr::select(x, l, u) %>% 
    mutate(
      x = case_when(row_number() ==1 ~ exp(x), 
                    rownames(df)=="zinc_statusnon_using:zinc1" ~ (exp(x+df$x[rownames(df)=="zinc1"]) - 1) * 100, 
                    rownames(df)=="zinc_statusnon_using:zinc2" ~ (exp(x+df$x[rownames(df)=="zinc2"]) - 1) * 100, 
                    rownames(df)=="study_month:zinc_statusnon_using" ~ (exp(x+df$x[rownames(df)=="study_month"]) - 1) * 100,
                    TRUE ~ (exp(x) - 1) * 100),
      l = case_when(row_number() ==1 ~ exp(l), 
                    rownames(df)=="zinc_statusnon_using:zinc1" ~ (exp(l+df$x[rownames(df)=="zinc1"]) - 1) * 100, 
                    rownames(df)=="zinc_statusnon_using:zinc2" ~ (exp(l+df$x[rownames(df)=="zinc2"]) - 1) * 100, 
                    rownames(df)=="study_month:zinc_statusnon_using" ~ (exp(l+df$x[rownames(df)=="study_month"]) - 1) * 100,
                    TRUE ~ (exp(l) - 1) * 100),
      u = case_when(row_number() ==1 ~ exp(u), 
                    rownames(df)=="zinc_statusnon_using:zinc1" ~ (exp(u+df$x[rownames(df)=="zinc1"]) - 1) * 100, 
                    rownames(df)=="zinc_statusnon_using:zinc2" ~ (exp(u+df$x[rownames(df)=="zinc2"]) - 1) * 100, 
                    rownames(df)=="study_month:zinc_statusnon_using" ~ (exp(u+df$x[rownames(df)=="study_month"]) - 1) * 100,
                    TRUE ~ (exp(u) - 1) * 100),
    ) 
  df_m3_weaner <- df_m3_weaner[!(rownames(df_m3_weaner) %in% c("sin_month", "cos_month")), ]
  
  # Transfoming finisher model 3 output
  df <- as.data.frame(intervals(m3_finisher, which = "fixed")[["fixed"]]) %>% rename("x"="est.", "l"="lower", "u"="upper")
  df_m3_finisher <- df %>% 
    dplyr::select(x, l, u) %>% 
    mutate(
      x = case_when(row_number() ==1 ~ exp(x), TRUE ~ (exp(x) - 1) * 100),
      l = case_when(row_number() ==1 ~ exp(l), TRUE ~ (exp(l) - 1) * 100),
      u = case_when(row_number() ==1 ~ exp(u), TRUE ~ (exp(u) - 1) * 100),
    ) 
  df_m3_finisher <- df_m3_finisher[!(rownames(df_m3_finisher) %in% c("sin_month", "cos_month")), ]
  
  # Printing results 
  cat("estimated fixed effects and 95% confidence intervals expressed as percentage change in AMU\n")
  cat("\nweaners:\n")
  print(round(as.matrix(df_m3_weaner), 1))
  cat("\nfinishers:\n")
  print(round(as.matrix(df_m3_finisher), 1))
})









################################################################################
################################################################################
################################################################################

### Predicted random effect (unknown factors) of the individual farms and veterinarians on AMU

# The best linear unbiased predictions (BLUPs) represent the predicted deviations of the AMU (log ADDkg/pig-day) from the 
  # average population (weaner or finisher) AMU (log ADDkg/pig-day) attributed to the individual farms or veterinarians. 


local({
  m3_models <- list(
    "weaner" = m3_weaner, 
    "finisher" = m3_finisher
  )
  
  df <- NULL
  
  for (i in seq_along(m3_models)) {
    m3_rearing_stage <- names(m3_models)[i]  
    m3 <- m3_models[[i]]
    
    # Predicted random effect of the individual farms
    df_farm <- coef(m3, level = 2)
    df_farm$names <- rownames(df_farm)
    rownames(df_farm) <- NULL
    df_farm <- df_farm %>% 
      rename("int_farm_vet"="(Intercept)") %>% 
      separate(names, into = c("vet_id", "farm_id"), sep = "/") %>%
      dplyr::select(vet_id, farm_id, int_farm_vet) %>% 
      mutate(
        rearing_stage = m3_rearing_stage,
        int_farm_vet = int_farm_vet - fixef(m3)["(Intercept)"]
      )
    
    # Predicted random effect of the individual vets
    df_vet <- coef(m3, level = 1) 
    df_vet$vet_id <- rownames(df_vet)
    rownames(df_vet) <- NULL
    df_vet <- df_vet %>% 
      rename("int_vet"="(Intercept)") %>% 
      dplyr::select(vet_id, int_vet) %>% 
      mutate(
        rearing_stage = m3_rearing_stage,
        int_vet = int_vet - fixef(m3)["(Intercept)"]
      ) %>% 
      left_join(df_farm, by=c("vet_id", "rearing_stage")) %>% 
      mutate(int_farm = int_farm_vet - int_vet) 
    
    # Combined df for random effects 
    df <- bind_rows(df, df_vet)
  }
  
  df <- df %>% mutate(rearing_stage = factor(rearing_stage, levels = c("weaner", "finisher")))
  
  
  # Plot of random effects
  ggplot(df,aes(x=reorder_within(vet_id, int_vet, rearing_stage), y=int_farm))+
    facet_grid(
      cols = vars(rearing_stage), scales = "free_x", 
      labeller = labeller(rearing_stage = c("weaner" = "Weaners", "finisher" = "Finishers"))
    )+
    geom_hline(yintercept = 0, color="azure4", size=1)+
    geom_point(size=0.75, shape=2, color="gray")+
    geom_point(aes(y=int_vet), color="black", size=1)+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_blank()
    )+
    labs(
      x="Veterinarians (with farms clustered by veterinarian)", 
      y="Random effect (BLUP)" 
    )
})









################################################################################
################################################################################
################################################################################

### Predicted random slopes (effect of zinc ban/discontinuation) on the individual farms

# For weaners, short-term and long-term effects of the zinc oxide ban/discontinuation were 
  # predicted seperatly for weaners on zinc-using and non-zinc using farms 



local({
  m3_models <- list(
    "weaner" = m3_weaner, 
    "finisher" = m3_finisher
  )
  
  effects <- c(
    "short_term", # short-term effect = 1-5 months after ban/discontinuation
    "long_term"   # long-term effect = >5 months after ban/discontinuation
  )
  
  df_effect <- NULL
  
  for (i in seq_along(m3_models)) {
    m3_rearing_stage <- names(m3_models)[i]  
    m3 <- m3_models[[i]]
    
    df_m3 <- coef(m3) 
    df_m3$names <- rownames(df_m3)
      
    for (j in effects) {
      if (j=="short_term") {
        x <- 2
        w1 <- "zinc_using1"
        w2 <- "zinc_non_using1"
        f <- "zinc1"
        fixed1 <- fixef(m3)["zinc1"]
        fixed2 <- fixef(m3)["zinc1"] + fixef(m3)["zinc_statusnon_using:zinc1"]
        
      } else if (j=="long_term") {
        x <- 3
        w1 <- "zinc_using2"
        w2 <- "zinc_non_using2"
        f <- "zinc2"
        fixed1 <- fixef(m3)["zinc2"]
        fixed2 <- fixef(m3)["zinc2"] + fixef(m3)["zinc_statusnon_using:zinc2"]
      }
      
      if (m3_rearing_stage == "weaner") {
        
        # Weaner farms with observations before and after ban/discontinuation and zinc status
        df_farm_weaner <- df_weaner %>% 
          dplyr::select(farm_id, zinc_status, zinc_using, zinc_non_using) %>% 
          group_by(farm_id) %>% 
          filter(n_distinct(zinc_using)>=x | n_distinct(zinc_non_using)>=x) %>% 
          ungroup() %>% 
          dplyr::select(farm_id, zinc_status) %>% 
          distinct()
        
        # Model output + weaner farm data
        df <- df_m3 %>% 
          separate(names, into = c("vet_id", "farm_id"), sep = "/") %>% 
          inner_join(df_farm_weaner, by = "farm_id")
        
        # Predicted effects in zinc using weaner farms 
        df_w_using <- df %>% 
          filter(zinc_status=="using") %>% 
          dplyr::select(farm_id, zinc_status, all_of(w1)) %>% 
          rename("slope" = w1) %>% 
          mutate(slope = slope + fixed1) %>% 
          mutate(slope = (exp(slope)-1)*100) %>% 
          mutate(effect = j)
        
        # Predicted effects in zinc non-using weaner farms 
        df_w_non_using <- df %>% 
          filter(zinc_status=="non_using") %>% 
          dplyr::select(farm_id, zinc_status, all_of(w2)) %>% 
          rename("slope" = w2) %>% 
          mutate(slope = slope + fixed2) %>% 
          mutate(slope = (exp(slope)-1)*100) %>% 
          mutate(effect = j)
        
        df_effect <- bind_rows(df_effect, df_w_using, df_w_non_using)
        
      } else if (m3_rearing_stage == "finisher") {
      
        # Finisher farms with observations before and after ban
        df_farm_finisher <- df_finisher %>%
          dplyr::select(farm_id, zinc) %>%
          group_by(farm_id) %>%
          filter(n_distinct(zinc)>=2) %>% 
          ungroup() %>% 
          dplyr::select(farm_id) %>% 
          distinct()
        
        # Predicted effects in finisher farms 
        df_f <- df_m3 %>%
          separate(names, into = c("vet_id", "farm_id"), sep = "/") %>%
          inner_join(df_farm_finisher, by = "farm_id") %>%
          dplyr::select(farm_id, all_of(f)) %>%
          rename("slope" = f) %>%
          mutate(slope = (exp(slope)-1)*100) %>% 
          mutate(zinc_status="f", effect=j)
      
        df_effect <- bind_rows(df_effect, df_f)
      }
    }
  }
  
  df_effect <- df_effect %>% mutate(effect = factor(effect, levels = c("short_term", "long_term")))
  
  # Plot of predicted effect of the zinc ban in individual farms 
  ggplot(df_effect, aes(y=slope, x=zinc_status)) +
    facet_grid(cols = vars(effect), labeller = as_labeller(c(
      "short_term"="1-5 months after the zinc oxide ban/discontinuation",
      "long_term"=">5 months after the zinc oxide ban/discontinuation"
    )))+
    geom_hline(yintercept = 0, color="black", linetype="dashed")+
    geom_jitter(aes(y = slope), width = 0.2, height = 0, alpha = 0.3, size=0.5, color="darkgrey")+
    geom_boxplot(outlier.shape = NA, alpha=0) + 
    labs(y="% change in AMU (ADDkg/pig-day)\ncompared to before the ban/biscontinuation", x="")+
    scale_x_discrete(labels = c(
      "non_using" = "Weaners\n(non-zinc-using)",
      "using" = "Weaners\n(zinc-using)",
      "f" = "Finishers"
    ), limits = c("non_using", "using", "f"))+
    theme_bw()+
    theme(legend.position = "none", 
          text = element_text(family = "serif", size = 10.5),
          strip.text.x = element_text(size = 10.5),
          axis.text = element_text(size = 10.5),
          panel.grid.major.x = element_blank()
    ) 
})










################################################################################
################################################################################
################################################################################


### model diagnostics 


# defining model, replace to get model diagnostics for finisher model
m3 <- m3_weaner

# residual vs fitted
local({std_residuals <- residuals(m3, type = "pearson")
plot(fitted(m3), std_residuals, main = "Residuals vs Fitted",
     xlab = "Fitted Values", ylab = "Standardized Residuals")})

# autocorrelation function plot
local({residuals_normalized <- residuals(m3, type = "normalized")
acf(residuals_normalized, lag.max = 12, main = "")
title(main = "ACF (Normalized Residuals)", line = 0.5)})

# scale-location plot
local({std_residuals <- residuals(m3, type = "pearson")
sqrt_abs_std_residuals <- sqrt(abs(std_residuals))
plot(fitted(m3), sqrt_abs_std_residuals,
     xlab = "Fitted Values", ylab = "Sqrt(|Standardized Residuals|)",
     main = "Scale-Location")})

# normal Q-Q plot (residuals)
local({residuals <- residuals(m3)
qqnorm(residuals, main = "Normal Q-Q Plot (Residuals)")
qqline(residuals, col = "red")})

# normal Q-Q plot (veterinarian random effect)
local({random_effect <-  ranef(m3)
qqnorm(random_effect[["vet_id"]][["(Intercept)"]], main = "Normal Q-Q Plot (Veterinarian random effect)")
qqline(random_effect[["vet_id"]][["(Intercept)"]], col = "red")})

# normal Q-Q plot (farm random effect)
local({random_effect <-  ranef(m3)
qqnorm(random_effect[["farm_id"]][["(Intercept)"]], main = "Normal Q-Q Plot (Farm random effect)")
qqline(random_effect[["farm_id"]][["(Intercept)"]], col = "red")})










