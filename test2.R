.libPaths( c( .libPaths(), "C:/Users/danau/R/Rlibs") )
library(tidyverse)
library(lubridate)
rm(list = ls())
load("C:/Users/danau/Dropbox/Dan's Playground/data_needed.Rdata")
data_reg <- data_needed %>%
    select(permno, date, adj_ret, size, b2m, e2p, exchcd, vwretd, at, be, ib)

som <- function(date) {
    floor_date(date, "month")
}

eom <- function(date) {
    ceiling_date(date, "month") - days(1)
}

tmp_beta <- data_reg %>%
    select(permno, date, adj_ret) %>%
    na.omit(adj_ret) %>%
    group_by(permno) %>%
    arrange(permno, date) %>%
    mutate(date = som(date)) %>%
    mutate(past5 = date %m-% months(60)) %>%
    mutate(num = row_number()) %>%
    mutate(index = date) %>%
    complete(index = seq(min(past5), max(date), by = "month")) %>%
    arrange(permno, index) %>%
    fill(num, .direction = "up")

tmp_beta1 <- tmp_beta
tmp_beta2 <- left_join(tmp_beta, tmp_beta1, by = c("permno", "index" = "past5")) %>%
    mutate(past5.y = index) %>%
    arrange(permno, date.y) %>%
    mutate(num = num.y - num.x) %>%
    mutate(date = eom(date.y)) %>%
    mutate(adj_ret = adj_ret.y) %>%
    select(permno, date, adj_ret, num)

data_beta <- tmp_beta2 %>%
    left_join(select(data_reg, permno, date, size, exchcd, vwretd, at, be, ib), by = c("permno", "date")) %>%
    drop_na(date) %>% 
    data.frame()

pre_beta <- function(df){
    df <- df %>%
        mutate(beta = NA) %>%
        mutate(beta1 = NA) %>%
        mutate(beta0 = NA) %>%
        mutate(vwretd_l1 = lag(vwretd)) %>%
        mutate(reg = ifelse((month(date) == 6) & (year(date) >= 1962) & (num >=24) & 
                                !is.na(at) & !is.na(be) & !is.na(ib) & (at != 0) & (be != 0), 1, NA))
    for (i in 1:nrow(df)) {
        if (!is.na(df$reg[i])) {
            df_tmp <- df[(i - df$num[i]):i, ]
            coeff <- lm(adj_ret ~ vwretd + vwretd_l1, data = df_tmp)$coef
            # print(coeff)
            df$beta0[i] <- coeff[2]
            df$beta1[i] <- coeff[3]
            df$beta[i] <- sum(coeff[2:3])
        }
    }
    df <- filter(df, !is.na(beta))
    df
}

preBetaList <- lapply(split(data_beta, data_beta$permno), pre_beta)
preBeta <- do.call(rbind, preBetaList)
nyse_preBeta <- preBeta %>% filter(exchcd %in% c(1, 31))

breakpoints <- function(df) {
    bps <- list()
    bps$sizes <- quantile(df$size, seq(0.1, 0.9, 0.1), na.rm = TRUE)
    df <- df %>% 
        mutate(sg = ntile(size, 10))
    dfls <- split(data.frame(df), df$sg)
    for (i in 1:length(dfls)) {
        bps[[ paste("beta_", as.character(dfls[[i]]$sg[1]), sep = "") ]] <- quantile(dfls[[i]]$beta, seq(0.1, 0.9, 0.1), na.rm = TRUE)
    }
    bps
}
nyse_breakpoints <- lapply(split(nyse_preBeta, nyse_preBeta$date), breakpoints)

preBetaGroup <- function(df, bps) {
    bps <- bps[[ as.character(df$date[1]) ]]
    gps <- df %>%
        mutate(size_g = findInterval(size, c(-Inf, bps$size, Inf)))
    gpl <- lapply(split(gps, gps$size_g), function(df) {
        df <- df %>% mutate(beta_g = findInterval(beta, c(-Inf, bps[[ paste("beta_", size_g[1], sep="") ]], Inf)))
    })
    do.call(rbind, gpl)
}
beta_groups <- do.call(rbind, lapply(split(preBeta, preBeta$date), preBetaGroup, nyse_breakpoints))

report_groups <- beta_groups %>% group_by(date, size_g, beta_g) %>% summarise(count = n())
report_groups2 <- beta_groups %>% group_by(size_g, beta_g) %>% summarise(count = n())

data_tmp <- data_beta %>% 
    # should be 6 months, there's an error in front
    mutate(testDate = floor_date(date %m-% months(6), "year") + months(6) - days(1)) %>%
    left_join(select(beta_groups, permno, date, size_g, beta_g), by = c("permno", "testDate" = "date"))
    
data_tmp2 <- data_tmp %>% 
    group_by(date, size_g, beta_g) %>% 
    summarize(avg_rtn = mean(adj_ret, na.rm = TRUE)) %>% 
    na.omit()

get_data <- function(query) {
    res <- dbSendQuery(wrds, query)
    data <- dbFetch(res, n=-1)
    dbClearResult(res)
    data
}

user = 'dan94'
password = 'TNafphzQBiM923W'

wrds <- dbConnect(Postgres(),
                  host='wrds-pgdata.wharton.upenn.edu',
                  port=9737,
                  dbname='wrds',
                  sslmode='require',
                  user = user,
                  password = password)

q_idx <- "select date, vwretd
          from crsp.msi"
get_data <- function(query) {
    res <- dbSendQuery(wrds, query)
    data <- dbFetch(res, n=-1)
    dbClearResult(res)
    data
}
data_idx <- get_data(q_idx)
adjust_date <- function(df) {
    df <- df %>% mutate(date = ceiling_date(date, "month") - days(1))
    df
}
data_idx <- adjust_date(data_idx) %>% mutate(vwretd_l1 = lag(vwretd))

data_tmp3 <- data_tmp2 %>% 
    left_join(data_idx, by = "date") %>% 
    mutate(g = paste(size_g, beta_g, sep = "_"))

data_post <- data_tmp3 %>%
    filter((date >= "1963-07-01") & (date < "1991-01-01"))

post_beta <- function(df){
    df_tmp <- df %>% na.omit()
    coeff <- lm(avg_rtn ~ vwretd + vwretd_l1, data = df_tmp)$coef
    # df <- filter(df, !is.na(beta))
    coeff[2:3]
}
    
postBeta <- as.data.frame(t(sapply(split(data_post, data_post$g), post_beta)))
postBeta$g <- rownames(postBeta)
postBeta <- postBeta %>%
    separate(g, c("size_g", "beta_g"), "_") %>% 
    mutate(beta = vwretd + vwretd_l1) %>% 
    select(size_g, beta_g, beta) %>% 
    spread(beta_g, beta) %>% 
    select(c(1, 2, 4:11, 3)) %>% 
    arrange(as.numeric(size_g))

t1pa <- data_tmp2 %>%
    filter((date >= "1963-07-01") & (date < "1991-01-01")) %>% 
    group_by(size_g, beta_g) %>% 
    summarize(avg_rtn = 100 * mean(avg_rtn)) %>% 
    spread(beta_g, avg_rtn)

t1pb <- postBeta

t1pc <- data_tmp %>% 
    filter((date >= "1963-07-01") & (date < "1991-01-01")) %>% 
    drop_na(size, size_g, beta_g) %>%
    group_by(size_g, beta_g) %>% 
    summarize(avg_size = mean(size)) %>% 
    spread(beta_g, avg_size)

test <- data_tmp %>% 
    select(date, permno, adj_ret, size_g, beta_g)

test3 <- as.data.frame(t(sapply(split(data_post, data_post$g), post_beta)))
test3$g <- rownames(test3)
test3 <- test3 %>%
    separate(g, c("size_g", "beta_g"), "_") %>% 
    mutate(beta = vwretd + vwretd_l1) %>% 
    select(size_g, beta_g, beta) %>% 
    mutate(size_g = as.numeric(size_g)) %>% 
    mutate(beta_g = as.numeric(beta_g))

test2 <- left_join(test, test3, by= c("size_g", "beta_g")) %>% 
    filter(date >= "1963-07-01" & date <= "1990-12-31")
beta_months <- split(test2, test2$date)

get_beta_slope <- function(df) {
    df <- drop_na(df, beta, adj_ret)
    coeff <- coef(lm(adj_ret *100 ~ beta, data=df))
    coeff
}

ttdf <- t(sapply(beta_months, get_beta_slope))
betas <- ttdf[,2]
mean(betas)
sd(betas)/sqrt(sum(length(betas)))

for (i in 1:738) {
    get_beta_slope(beta_months[[i]])
}

# testLs = lapply(nyse_prebeta, testFunc)
# 
# df = test
# bps = nyse_breakpoints
# bps <- bps[[ as.character(df$date[1]) ]]
# gps <- df %>%
#     mutate(size_g = findInterval(size, c(-Inf, bps$size, Inf)))
# gpl <- split(gps, gps$size_g)
# lapply(gpl, function(df) {
#     df <- df %>% mutate(beta_g = findInterval(beta, c(-Inf, bps[[ paste("beta_", size_g[1], sep="") ]], Inf)))
# })
# do.call(rbind, gpl)
# 
# test2 = gpl[[1]]
# test2 = test2 %>% select(permno, date, size_g, beta)
# View(test2 %>% mutate(beta_g = findInterval(beta, c(-Inf, bps[[ paste("beta_", size_g[[1]], sep="") ]], Inf))) %>% arrange(beta))
# 
# 
# preBetaBps <- lapply(split(), breakpoints)
# 
# nyse_preBeta <- do.call(rbind, nysePreBetaList)
# an_preBeta <- do.call(rbind, anPreBetaList)
# 
# preBetaGroups <- preBeta %>%
#     select(permno, date, size, beta) %>% 
#     split(preBeta$date)
# 
# 
# a = preBetaGroup(test, nyse_breakpoints)
# a %>% mutate(beta_g = paste("beta_", size_g, sep = ""))
# 
# test1 <- test %>% mutate(sg = ntile(size, 10))
# test2 <- test1 %>% mutate(sg2 = cut)
# test2 <- test1 %>% mutate(sg2 = findInterval(size, c(-Inf, bps$size, Inf)))
# bps = nyse_breakpoints[[1]]
# 
# factor(findInterval(test1$size, c(-Inf, quantile(test1$size, seq(0.1, 0.9, 0.1), na.rm = TRUE), Inf)), labels=1:10)
# 
# #
# # test <- data_beta %>% filter(permno == 10147)
# # 
# # test <- test %>%
# #     mutate(beta = NA) %>%
# #     mutate(beta1 = NA) %>%
# #     mutate(beta0 = NA) %>%
# #     mutate(vwretd_l1 = lag(vwretd)) %>%
# #     mutate(reg = ifelse((exchcd %in% c(1, 31)) & (month(date) == 7) & (year(date) >= 1962) & (num >=24) & !is.na(at) & !is.na(be) & !is.na(ib) & (at != 0) & (be != 0), 1, NA))
# # for (i in 1:nrow(test)) {
# #     if (!is.na(test$reg[i])) {
# #         df_tmp <- test[(i - test$num[i]):i, ]
# #         coeff <- lm(adj_ret ~ vwretd + vwretd_l1, data = df_tmp)$coef
# #         print(coeff)
# #         test$beta0[i] <- coeff[2]
# #         test$beta1[i] <- coeff[3]
# #         test$beta[i] <- sum(coeff[2:3])
# #     }
# # }
# # 
# # View(test %>% filter(!is.na(beta)))
# 
# # save(data_beta, file = "C:/Users/danau/data_beta.Rdata")
# # # load("C:/Users/danau/data_beta.Rdata")
# 
# nyseNames <- unique(filter(data_beta, exchcd %in% c(1, 31))$permno)
# nyse <- data_beta %>%
#     filter(permno %in% nyseNames) %>%
#     data.frame()
# 
# amexNasNames <- unique(filter(data_beta, !(exchcd %in% c(1, 31)))$permno)
# amexNas <- data_beta %>%
#     filter(permno %in% amexNasNames) %>%
#     data.frame()
# 
# 
# 
