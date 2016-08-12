## ----setup, warning=FALSE, message=FALSE---------------------------------
library(tcpl)
library(toxboot)
library(data.table)
library(RMySQL)
library(DBI)
library(magrittr)
library(ggplot2)
library(pander)

## ----erl3data_code, eval = FALSE-----------------------------------------
#  tcplConf(db = "prod_external_invitrodb_v2")
#  assay_names <- c("NVS_NR_bER",
#                   "OT_ER_ERaERa_1440",
#                   "ATG_ERa_TRANS_up",
#                   "TOX21_ERa_LUC_BG1_Agonist",
#                   "ACEA_T47D_80hr_Positive")
#  
#  aeid_table_full <- tcplLoadAeid()
#  aeid_table <- aeid_table_full[aenm %in% assay_names]
#  aeids <- aeid_table[,aeid]
#  
#  dat <- toxbootQueryToxCast(aeids = aeids)
#  
#  set.seed(12345)
#  m4ids <- sample(unique(dat[, m4id]), size = 200)
#  erl3data <- dat[m4id %in% m4ids]

## ----memory_toxboot, echo=TRUE, eval = TRUE------------------------------
dat <- toxbootmc(dat = erl3data, 
                 boot_method = "smooth", 
                 m4ids = tail(erl5data[hitc == 1L, m4id], 10),
                 cores = 1, 
                 destination = "memory", 
                 replicates = 10)
dim(dat)

## ----file_toxboot, echo=TRUE, eval = FALSE-------------------------------
#  toxbootmc(dat = erl3data,
#            boot_method = "smooth",
#            cores = 8,
#            destination = "file",
#            replicates = 10)

## ----mongo_toxboot, eval = FALSE-----------------------------------------
#  toxbootConf(mongo_host = "123.45.67.89",
#              db = "bootstrap",
#              DBNS = "bootstrap.prod_external_invitrodb_v2",
#              user = "username",
#              pass = "password")
#  
#  toxbootmc(dat = erl3data,
#            boot_method = "smooth",
#            cores = 8,
#            destination = "mongo",
#            replicates = 10)

## ----mongo_toxboot_query, eval = FALSE-----------------------------------
#  m4ids <- unique(erl3data[, m4id])
#  fields <- c("m4id", "max_med", "hill_ga", "hill_gw", "hill_tp", "hill_aic",
#              "gnls_ga", "gnls_gw", "gnls_tp",  "gnls_la", "gnls_lw", "gnls_aic",
#              "cnst_aic")
#  dat_boot <- toxbootGetMongoFields(m4id = m4ids, fields = fields)

## ----mysql_make_toxboot, echo = TRUE, eval = FALSE-----------------------
#  toxbootMysqlCreateTable()

## ----mysql_toxboot, echo=TRUE, eval = FALSE------------------------------
#  toxbootmc(dat = erl3data,
#            boot_method = "smooth",
#            cores = 32,
#            destination = "mysql",
#            replicates = 10)
#  
#  dat_boot <- toxbootGetMySQLFields()

## ----erl5data_command, eval = FALSE--------------------------------------
#  m4ids <- unique(erl3data[, m4id])
#  erl5data <- tcplLoadData(5, fld = "m4id", val = m4ids, type = "mc")

## ----modl_hit------------------------------------------------------------
dat_tb <- toxbootHitParamCI(dat, erl5data)

## ----hit_pct_plot--------------------------------------------------------
dat_sum <- dat_tb[, .(hit_pct = sum(boot_hitc)/10), by = m4id]
dat_sum

ggplot(dat_sum, 
       aes(x = hit_pct)) + 
  geom_histogram(binwidth = 0.1) + 
  theme_bw()

## ----hit_pct_ecdf--------------------------------------------------------
ggplot(dat_sum, 
       aes(x = hit_pct)) + 
  stat_ecdf() + 
  theme_bw()

## ----parameter_tables----------------------------------------------------
pander(erl5data[m4id == 9057756, .(modl, 
                                   hill_ga, 
                                   hill_gw, 
                                   hill_tp, 
                                   gnls_ga, 
                                   gnls_gw, 
                                   gnls_tp, 
                                   gnls_la, 
                                   gnls_lw)], 
       split.table = Inf)

## ----pipeline_plot-------------------------------------------------------
  hill_ga <- erl5data[m4id == 9057756, hill_ga]
  hill_gw <- erl5data[m4id == 9057756, hill_gw]
  hill_tp <- erl5data[m4id == 9057756, hill_tp]
  gnls_ga <- erl5data[m4id == 9057756, gnls_ga]
  gnls_gw <- erl5data[m4id == 9057756, gnls_gw]
  gnls_tp <- erl5data[m4id == 9057756, gnls_tp]
  gnls_la <- erl5data[m4id == 9057756, gnls_la]
  gnls_lw <- erl5data[m4id == 9057756, gnls_lw]
  
ggplot(erl3data[m4id == 9057756],
       aes(x=logc, 
           y=resp)) +
  stat_function(fun = hill_curve, 
                args=list(hill_tp = hill_tp, 
                          hill_ga = hill_ga, 
                          hill_gw = hill_gw),
                alpha = 1,
                color = "red", 
                size = 1) +
  stat_function(fun = gnls_curve, 
                args=list(top = gnls_tp, 
                          ga = gnls_ga, 
                          gw = gnls_gw, 
                          la = gnls_la, 
                          lw = gnls_lw),
                alpha = 1, 
                color = "blue", 
                size = 1,
                linetype = 2) +
  theme_bw() +
  geom_point(size=5,alpha=1) +
  theme(legend.position="none", legend.title=element_blank()) +
  ylab("Percent Activity") +
  xlab("Log Concentration")

## ----9057756-------------------------------------------------------------
ggplot(dat_tb[m4id == 9057756], 
       aes(x = modl_ga)) + 
  stat_ecdf() + 
  theme_minimal()

## ----9057756_1000_10000, eval = FALSE------------------------------------
#  dat1000 <- toxbootmc(dat = erl3data,
#                       m4ids = rep(9057756, 8),
#                       boot_method = "smooth",
#                       cores = 8,
#                       destination = "memory",
#                       replicates = 125) %>%
#    toxbootHitParamCI(erl5data)
#  
#  dat10000 <- toxbootmc(dat = erl3data,
#                        m4ids = rep(9057756, 8),
#                        boot_method = "smooth",
#                        cores = 8,
#                        destination = "memory",
#                        replicates = 1250) %>%
#    toxbootHitParamCI(erl5data)

## ----read_1000_10000-----------------------------------------------------
dim(dat1000)
dim(dat10000)

## ----plot_1000_10000-----------------------------------------------------
ggplot(dat10000, 
       aes(x = modl_ga)) + 
  stat_ecdf() + 
  stat_ecdf(data = dat1000,
            color = "blue",
            linetype = 2) +
  stat_ecdf(data = dat_tb[m4id == 9057756],
            color = "red",
            linetype = "dotdash") +
  theme_bw()

## ----boot_fits-----------------------------------------------------------
rep_num <- 1000
xmin <- min(erl3data[m4id == 9057756, logc])
xmax <- max(erl3data[m4id == 9057756, logc])

dat_boot_curve <- expand.grid(replicate = 1:rep_num,
                              lconc = seq(xmin,
                                          xmax,
                                          length.out = 100)) %>%
  data.table()
dat_result <- copy(dat1000)
dat_result[, repnum := 1:.N]
dat_boot_curve <- merge(dat_boot_curve,
                        dat_result,
                        by.x = "replicate",
                        by.y = "repnum")
dat_boot_curve[modl == "hill", 
               resp := hill_curve(hill_tp = hill_tp, 
                                  hill_ga = hill_ga, 
                                  hill_gw = hill_gw, 
                                  lconc)]
dat_boot_curve[modl == "gnls", 
               resp := gnls_curve(top = gnls_tp, 
                                  ga = gnls_ga, 
                                  gw = gnls_gw, 
                                  la = gnls_la, 
                                  lw = gnls_lw, 
                                  lconc)]

hill_ga <- erl5data[m4id == 9057756, hill_ga]
hill_gw <- erl5data[m4id == 9057756, hill_gw]
hill_tp <- erl5data[m4id == 9057756, hill_tp]
gnls_ga <- erl5data[m4id == 9057756, gnls_ga]
gnls_gw <- erl5data[m4id == 9057756, gnls_gw]
gnls_tp <- erl5data[m4id == 9057756, gnls_tp]
gnls_la <- erl5data[m4id == 9057756, gnls_la]
gnls_lw <- erl5data[m4id == 9057756, gnls_lw]

ggplot(dat_boot_curve, 
                      aes(x=lconc, 
                          y=resp,
                          color = modl)) +
  geom_line(size = 2,
            alpha = 0.01,
            aes(group = replicate)) +
  geom_point(data = erl3data[m4id == 9057756],
             aes(x = logc,
                 y = resp),
             alpha = 1,
             size = 5,
             color = 'black', 
             fill = 'cyan', 
             shape=21) +
  stat_function(fun = hill_curve, 
                args=list(hill_tp = hill_tp, 
                          hill_ga = hill_ga, 
                          hill_gw = hill_gw),
                alpha = 1,
                color = "cyan", 
                size = 1) +
  scale_color_manual(values = c("hill" = "red", "gnls" = "blue")) +
  ylab("Percent Activity") +
  xlab("Log Concentration (uM)") +
  expand_limits(y = c(120, -40)) +
  theme_bw() +
  guides(color=FALSE)

