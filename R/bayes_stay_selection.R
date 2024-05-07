#' Model selection for the residence time model based on WAIC for REST/REST-RAD models
#'
#' @param formula_stay Model formula for the residence time within a focal area. Grammar follows lme4::glmer function. e.g. Stay ~ 1 + (1|Group)
#' @param station_data A data frame containing information for each camera station in each row
#' @param stay_data A data frame containing the time spent within a focal area
#' @param family The probability distribution of the time spent. Specify the probability distribution selected by the bayes_stay_selection function.
#' @param local_stay Specify whether to use the expected value of the time spent for each camera station (local_stay = TRUE) or the global expected value (FALSE) for density estimation.
#' @param plot If TRUE, plots the expected values of residence times.
#' @param cores The number of cores to use for parallel computation. The default value is 1.
#' @param iter The number of iterations. The default value is 2000.
#' @param warmup The number of warmup iterations. The default value is NULL (half of the iterations will be used for warmup).
#' @param chains The number of chains. The default value is 2.
#' @param thin The thinning interval. The default value is 1.
#' @param all_comb If TRUE, models with all combinations of covariates are compared. If FALSE, the WAIC value of the designated model is shown.
#' @return WAIC values for each model
#' @import dplyr ggplot2
#' @export
#' @examples
#' stay_data <- format_stay(
#'   detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens",
#'   target_species = "A"
#' )
#' bayes_stay_selection(formula_stay = Stay ~ 1 + x1,
#'                     station_data,
#'                     stay_data,
#'                     family = "lognormal",
#'                     local_stay = TRUE,
#'                     plot = TRUE,
#'                     cores = 3,
#'                     iter = 2000,
#'                     warmup = NULL,
#'                     chains = 2,
#'                     thin = 1,
#'                     all_comb = FALSE)

bayes_stay_selection <- function(formula_stay,
                                 station_data,
                                 stay_data,
                                 family = "lognormal",
                                 local_stay = FALSE,
                                 plot = TRUE,
                                 cores = 3,
                                 iter = 2000,
                                 warmup = NULL,
                                 chains = 2,
                                 thin = 1,
                                 all_comb = TRUE){

  #######################
  ###private functions###
  #######################
  centering <- function(data,xname_stay){
    for(i in 1:length(xname_stay)){
      temp <- regexpr(":",xname_stay[i],fixed=TRUE)
      if(temp[1]>0){
        for(j in 1:length(data)){
          if(names(data[j])==xname_stay[i] && names(data[j]) != "(intercept)"){
            data[j] <- data[j]-mean(as.matrix(data[j]),na.rm=T)
          }
        }
      }
    }
    return(data)
  }
  nobars <- function(term) {
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
      nb <- nobars(term[[2]])
      if (is.null(nb)) return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }

  findbars <- function(term) {
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
  }
  subbars <- function(term){
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
      term[[2]] <- subbars(term[[2]])
      return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
      term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
  }
  dvcalc <- function(data,yname_stay){
    yname2 <- unlist(strsplit(yname_stay,"+",fixed=TRUE))
    yname2 <- unlist(strsplit(yname2,"-",fixed=TRUE))
    yname2 <- unlist(strsplit(yname_stay,"*",fixed=TRUE))
    yname2 <- unlist(strsplit(yname2,"/",fixed=TRUE))
    if(length(yname2)>1){
      lyn2 <- length(yname2)-1
      count <- 0
      for(i in 1:lyn2){
        temp <- substr(yname_stay,count+nchar(yname2[i])+1,count+nchar(yname2[i])+1)
        if(temp=="+"){
          data[yname2[1]] <- data[yname2[1]] + data[yname2[i+1]]
        }else if(temp=="-"){
          data[yname2[1]] <- data[yname2[1]] - data[yname2[i+1]]
        }else if(temp=="*"){
          data[yname2[1]] <- data[yname2[1]] * data[yname2[i+1]]
        }else if(temp=="/"){
          data[yname2[1]] <- data[yname2[1]] / data[yname2[i+1]]
        }
        count <- count + 1 + nchar(yname2[i])
      }
    }
    return(data)
  }

  bit.test <- function(number, n) {
    (number %/% (2^n)) %% 2
  }

  full_terms <- function(x) {
    lapply(1:2^length(x),
           function(i) {
             r <- ""
             for (j in 1:length(x)) {
               r <- paste(r,
                          ifelse(bit.test((i - 1), (j - 1)),
                                 paste(x[j], " + ", sep=""), ""),
                          sep="")
             }
             r <- paste("Stay ~ 1 + ", r, "1", sep = "") #切片がベータの場合
             r <- strsplit(r, " \\+ 1$")[[1]][1]
           }
    )
  }

  #################
  ###option ?###
  #################
  center = TRUE
  cauchy = 2.5
  lkj_corr = 2
  #################
  ###check input###
  #################

  if(family=="lognormal" || family=="gamma" || family=="weibull" || family=="exponential"){

  }else{
    stop(paste0("Input family type(",family,") is incorrect."))
  }

  if(is.null(warmup)==TRUE){
    warmup = floor(iter/2)
  }

  if(is.null(cores)==TRUE){
    cores <- chains
  }
  if(cores>getOption("mc.cores",detectCores())){
    cores <- getOption("mc.cores",detectCores())
  }
  #####################
  ### convert original dataframe
  #####################

  stay_data <- data.frame(stay_data)
  stay_data_join <- stay_data %>%
    left_join(station_data, by = "Station")

  stay_data_join_temp <- stay_data_join

  ####################
  ### Creating full model
  ###################
  formula_stay <- as.formula(formula_stay)
  formula_stay <- paste0(formula_stay[2],formula_stay[1],formula_stay[3])

  # 文字列を空白で分割する

  if(all_comb == TRUE){
    terms <- unlist(strsplit(formula_stay, "~|\\+"))
    terms <- trimws(terms)
    terms <- terms[terms != "1"][-1]
    t <- full_terms(terms)
  } else {
    t <- formula_stay
  }

  fitstan <- list(0); waic2 <- lppd <- p_waic <- numeric(0)

  for(h in 1:length(t)){
    #####################
    ###variable names for stay
    #####################
    stay_data_join <- stay_data_join_temp

    formula_stay <- as.formula(t[[h]])
    formula_stay_paste <- paste0(formula_stay[2],formula_stay[1],formula_stay[3])
    xformula_stay <- nobars(formula_stay)
    if (is.name(formula_stay) && length(formula_stay) == 1) {
      formula_stay <- reformulate(formula_stay, response = NULL)
    }
    xname_stay <- colnames(model.matrix(xformula_stay, stay_data_join))

    # 文字列を空白で分割する
    terms <- unlist(strsplit(formula_stay_paste, "~|\\+"))
    terms <- trimws(terms)
    terms <- terms[terms != "1"][-1]

    yname_stay <- deparse(xformula_stay[[2]] )
    stay <- model.frame(xformula_stay, stay_data_join)[,1]


    zname_stay <- list()
    if ( nchar(formula_stay[3]) != nchar(nobars(formula_stay)[3]) ) {
      zformula <- findbars(formula_stay)
      for ( i in 1:length(zformula) ) {
        idname_stay <- all.vars( zformula[[i]][[3]])
        v <- zformula[[i]][[2]]
        if (is.numeric(v)) {
          zname_stay[[idname_stay]] <- "(Intercept)"
        } else {
          tempv <- gsub(" ", "", deparse( v ),fixed=TRUE)
          tempv <- unlist(strsplit(tempv,"+",fixed=TRUE))
          if(tempv[1]=="1"){
            f <- as.formula( paste( "~" , deparse( v ) , sep="" ) )
            zname_stay[[idname_stay]] <- colnames( model.matrix( f , stay_data_join ) )
          }else{
            f <- as.formula( paste( "~" , deparse( v ) , sep="" ) )
            f <- as.data.frame(model.matrix( f , stay_data_join ))[-1]
            zname_stay[[idname_stay]] <- colnames( f )
          }
        }
      }
    }

    idname_stay <- list()
    R_stay <- length(zname_stay)
    P_stay <- length(xname_stay)
    idname_stay <- c()
    if(R_stay>0) for(i in 1:R_stay) idname_stay[i] <- attr(zname_stay,"name")[i]

    ######################
    ###creating dataset###
    ######################
    # 滞在時間
    dat2 <- stay_data_join
    dat2["(Intercept)"] <- 1
    if(center == TRUE) dat2 <- centering(dat2, xname_stay)
    datname <- c(yname_stay, xname_stay)
    if(R_stay > 0) for(i in 1:R_stay) datname <- c(datname, idname_stay[i], zname_stay[[i]])

    x_stay_expand <- model.matrix(xformula_stay, dat2)
    dat3 <- subset(cbind(x_stay_expand, dat2), select = unique(datname)) # ランダム効果も含むモデル式内の変数

    x_meanstay_covs <- sapply(data.frame(x_stay_expand), mean)
    names(x_meanstay_covs) <- NULL
    #　サイトレベルでの滞在時間の共変量
    station_data["Stay"] <- 0
    x_stay_station <- data.frame(model.matrix(xformula_stay, station_data))
    names(x_stay_station)[names(x_stay_station) == "X.Intercept."] <- "(Intercept)"

    # 応答変数（滞在時間）
    stay <- as.numeric(dat3[ , yname_stay])
    censored <- stay_data_join[ , "Cens"]
    x_stay <- subset(dat3, select = xname_stay)
    N_station <- nrow(station_data)
    # 滞在時間のランダム効果
    z_stay <- list()
    id_stay <- list()
    Q_stay <- array()
    G_stay <- array()
    if(R_stay>0){
      for(i in 1:R_stay){
        Q_stay[i] <- length(zname_stay[[i]])
        if(Q_stay[i]==1){
          z_stay[[i]] <- subset(dat3,select = zname_stay[[i]])[,]
        }else {
          z_stay[[i]] <- subset(dat3,select = zname_stay[[i]])
        }
        id_stay[[i]] <- as.numeric(as.factor(as.numeric(as.factor(dat3[,idname_stay[i]]))))
        G_stay[i] <- length(unique(id_stay[[i]]))
      }
    }else{
      Q_stay <- c(0,0,0,0)
      G_stay <- c(0,0,0,0)
    }
    N_stay <- nrow(dat3)

    z2name_stay <- c()
    id2name_stay <- c()
    for(i in 1:R_stay){
      z2name_stay[i] <- paste("z_stay", i, sep = "")
      id2name_stay[i] <- paste("id_stay", i, sep = "")
    }

    # datalistの作成
    if (R_stay > 0) {
      datastan <- list(N_stay = N_stay, P_stay = P_stay, R_stay = R_stay, Q_stay = Q_stay, G_stay = G_stay, stay = stay, x_stay = x_stay, censored = censored, N_station = N_station)
    } else {
      datastan <- list(N_stay = N_stay, P_stay = P_stay, stay = stay, x_stay = x_stay, censored = censored, N_station = N_station)
    }

    if(R_stay>0){
      for(i in 1:R_stay){
        datastan[z2name_stay[i]] <- list(z_stay[[i]])
        datastan[id2name_stay[i]] <- list(id_stay[[i]])
      }
    }

    if(local_stay == TRUE){
      datastan["x_stay_station"] <- list(x_stay_station)
    }

    if(local_stay == FALSE){
      datastan["x_meanstay_covs"] <- list(x_meanstay_covs)
    }

    ##########################
    ####creating stan code####
    ##########################

    ###data
    data_code <- "data{\n"
    temp1 <- "\tint<lower=1> N_stay;\n"
    temp1 <- paste0(temp1, "\tvector[N_stay] stay;\n")
    temp1 <- paste0(temp1, "\tarray[N_stay] int<lower=0, upper = 1> censored;\n")

    temp1 <- paste0(temp1, "\tint<lower=1> P_stay;\n")
    temp1 <- paste0(temp1, "\trow_vector[P_stay] x_stay[N_stay];\n")

    if(R_stay > 0){
      temp1 <- paste0(temp1,"\tint<lower=1> R_stay;\n\tint<lower=1> G_stay[R_stay];\n\tint<lower=1> Q_stay[R_stay];\n" )
    }

    if(R_stay > 0){
      for(i in 1:R_stay){
        temp1 <- paste0(temp1,"\t","int<lower=1> ", "id_stay",i,"[N_stay];","\n")
      }
    }

    temp2 <- ""
    if(R_stay > 0){
      for(i in 1:R_stay){
        if(Q_stay[i]==1){
          temp2 <- paste0(temp2,"\t","real", " z_stay",i,"[N_stay];","\n")
        }else{
          temp2 <- paste0(temp2,"\t","row_vector[Q_stay[",i, "]] z_stay",i,"[N_stay];","\n")
        }
      }
    }

    #  要追加項目
    temp3 <- "\tint<lower=1>  N_station;\n"

    if(local_stay == TRUE){
      temp3 <- paste0(temp3, "\trow_vector[P_stay] x_stay_station[N_station];\n")
    }

    if(local_stay == FALSE){
      if(length(x_meanstay_covs) > 1) temp3 <- paste0(temp3, "\trow_vector[P_stay] x_meanstay_covs;\n");
    }

    data_code <- paste0(data_code,temp1,temp2,temp3,"}\n")

    ###transformed data
    td_code <- "transformed data{\n"

    if(R_stay>0){
      temp1 <- ""
      for(i in 1:R_stay){
        temp1 <- paste0(temp1,"\t","vector[Q_stay[",i, "]] mu_stay",i,";\n")
      }
      temp2 <- ""
      for(i in 1:R_stay){
        temp2 <- paste0(temp2,"\t","for(q in 1:Q_stay[",i,"]) mu_stay",i,"[q] = 0;\n")
      }
      td_code <- paste0(td_code,temp1,temp2,"}\n")
    }else{
      td_code <- paste0(td_code,"}\n")
    }

    ###parameters
    para_code <-'parameters{\n'
    temp1 <- "\tvector[P_stay] beta;\n"

    if(R_stay>0){
      for(i in 1:R_stay){
        if(Q_stay[i]==1){
          temp1 <- paste0(temp1,"\t","real r",i,"[G_stay[",i,"]];\n")
        }else{
          temp1 <- paste0(temp1,"\t","vector[Q_stay[",i, "]] r",i,"[G_stay[",i,"]];\n")
        }
      }
    }
    temp2 <- ''
    if(R_stay>0){
      for(i in 1:R_stay){
        if(Q_stay[i]==1){
          temp2 <- paste0(temp2,"\t","real<lower=0> tau_sd",i,";\n")
        }else{
          temp2 <- paste0(temp2,"\t","vector<lower=0>[Q_stay[",i, "]] tau_sd",i,";\n")
        }
      }
      for(i in 1:R_stay){
        if(Q_stay[i]>1){
          temp2 <- paste0(temp2,"\t","corr_matrix[Q_stay[",i, "]] tau_corr",i,";\n")
        }
      }
    }
    temp3 <- ''
    if(family == "lognormal" || family == "gamma" | family == "weibull"){
      temp3 <- paste0("\t","real<lower=0> s;\n")
    }
    para_code <- paste0(para_code,temp1,temp2,temp3,"}\n")

    ###transformed parameters
    tp_code <-'transformed parameters{\n'

    temp1 <- "\treal lambda[N_stay];\n"
    temp1 <- paste0(temp1, "\treal predict[N_stay];\n")

    if(family == "lognormal" || family == "gamma" | family == "weibull"){
      temp1 <- paste0(temp1,"\t","real<lower=0> scale;\n")
    }

    temp2 <- ''
    if(R_stay>0){
      for(i in 1:R_stay){
        if(Q_stay[i]==1){
          temp2 <- paste0(temp2,"\t","real<lower=0> tau",i,";\n")
        }else{
          temp2 <- paste0(temp2,"\t","cov_matrix[Q_stay[",i, "]] tau",i,";\n")
        }
      }
    }

    temp3 <- ''
    if(family == "gamma"){
      temp3 <- paste0(temp3,"\tscale = 1/s;\n")
    }else if(family == "lognormal"){
      temp3 <- paste0(temp3,"\tscale = s^2;\n")
    }else if(family == "weibull"){
      temp3 <- paste0(temp3,"\tscale = 1/s;\n")
    }
    temp4 <- ''
    if(R_stay>0){
      for(i in 1:R_stay){
        if(Q_stay[i]==1){
          temp4 <- paste0(temp4,"\t","tau",i," = tau_sd",i,"^2;\n")
        }else{
          temp4 <- paste0(temp4,"\ttau",i," <- quad_form_diag(tau_corr",i,",tau_sd",i,");\n")
        }
      }
    }

    # temp1
    temp5 <- '\tfor(n in 1:N_stay){\n'
    temp5 <- paste(temp5, "\t\tpredict[n] = x_stay[n]*beta")
    if(R_stay>0){
      for(i in 1:R_stay){
        temp5 <- paste0(temp5,"+z_stay",i,"[n]*r",i,"[id_stay",i,"[n]]")
      }
    }
    temp5 <-paste0(temp5,";\n")

    if(family=="exponential"){
      temp5 <- paste0(temp5,"\t\tlambda[n] = 1 / exp(predict[n]);\n") # 左辺のpredictはパラメータ
    }else if(family=="gamma"){
      temp5 <- paste0(temp5,"\t\tlambda[n] = s / exp(predict[n]);\n") # 分母が期待値
    }else if(family=="lognormal"){
      temp5 <- paste0(temp5,"\t\tlambda[n] = log(exp(predict[n])) - (s*s/2);\n") # 分母が期待値
    }else if(family=="weibull"){
      temp5 <- paste0(temp5,"\t\tlambda[n] = 1 / (log(s) - predict[n]);\n")
    }
    temp5 <-paste0(temp5,"\t}\n")

    # mean_stay
    if(local_stay == TRUE){
      temp5 <- paste0(temp5, "\tvector[N_station] mean_stay;\n")
      temp5 <- paste0(temp5, "\tfor(j in 1:N_station){\n")
      temp5 <- paste0(temp5,   "\t\tmean_stay[j] = exp(x_stay_station[j] * beta)")
      if(R_stay>0){
        for(i in 1:R_stay){
          temp5 <- paste0(temp5," * exp(z_stay",i,"[j] * r",i,"[id_stay",i,"[j]])")
        }
      }
      temp5 <- paste0(temp5, ";\n")
      temp5 <- paste0(temp5, "\t}\n")
    }
    if(local_stay == FALSE){
      temp5 <- paste0(temp5, "\treal mean_stay;\n");
      if(length(x_meanstay_covs) > 1) temp5 <- paste0(temp5, "\tmean_stay = exp(x_meanstay_covs * beta);\n")
      if(length(x_meanstay_covs) == 1) temp5 <- paste0(temp5, "\tmean_stay = exp(beta[1]);\n")
    }
    tp_code <- paste0(tp_code,temp1,temp2,temp3,temp4,temp5,"}\n")


    ###model
    model_code <- 'model{\n'

    temp3 <-"\tbeta ~ normal(0,100);\n"
    if(R_stay>0){
      for(i in 1:R_stay){
        if(Q_stay[i]==1){
          temp3 <- paste0(temp3,"\t","for(g in 1:G_stay[",i,"]) r",i,"[g] ~ normal(0, tau_sd",i,");\n")
        }else{
          temp3 <- paste0(temp3,"\t","r",i," ~ multi_normal(mu",i,",tau",i,");\n")
        }
      }
      for(i in 1:R_stay){
        if(Q_stay[i]>1){
          temp3 <- paste0(temp3,"\t","tau_corr",i," ~ lkj_corr(",lkj_corr,");\n")
        }
      }
    }
    if(family=="lognormal" | family=="gamma" | family == "weibull"){ # これは自分で足した（Gammaについて）
      temp3 <- paste0(temp3,"\t","s ~ cauchy(0,",cauchy,");\n")
    }

    temp2 <-'\t'
    if(family == "exponential"){
      temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
      temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
      temp2 <- paste0(temp2, "\t\ttarget += exponential_lpdf(stay[n] | lambda[n]);\n")
      temp2 <- paste0(temp2,"\t} else {\n")
      temp2 <- paste0(temp2, "\t\ttarget += exponential_lccdf(stay[n] | lambda[n]);\n")
    }else if(family=="gamma"){
      temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
      temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
      temp2 <- paste0(temp2, "\t\ttarget += gamma_lpdf(stay[n] | s, lambda[n]);\n")
      temp2 <- paste0(temp2,"\t} else {\n")
      temp2 <- paste0(temp2, "\t\ttarget += gamma_lccdf(stay[n] | s, lambda[n]);\n")
    }else if(family=="weibull"){
      temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
      temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
      temp2 <- paste0(temp2, "\t\ttarget += weibull_lpdf(stay[n] |lambda[n], s);\n")
      temp2 <- paste0(temp2,"\t} else {\n")
      temp2 <- paste0(temp2, "\t\ttarget += weibull_lccdf(stay[n] |lambda[n], s);\n")
    }else if(family=="lognormal"){
      temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
      temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
      temp2 <- paste0(temp2, "\t\ttarget += lognormal_lpdf(stay[n] |lambda[n], s);\n")
      temp2 <- paste0(temp2, "\t} else {\n")
      temp2 <- paste0(temp2, "\t\ttarget += lognormal_lccdf(stay[n] |lambda[n], s);\n")
    }
    temp2 <- paste0(temp2, "\t\t}\n\t}\n")
    model_code <- paste0(model_code,temp3,temp2,"}\n")

    ###generated quantities

    gq_code <-'generated quantities{\n'

    temp2 <- ''
    temp2 <- paste0(temp2, "\n\treal log_lik_stay[N_stay];\n")
    temp2 <- paste0(temp2, '\tfor(n in 1:N_stay){\n')

    if(family == "exponential"){
      temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = exponential_lpdf(stay[n] | lambda[n]);\n")
      temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = exponential_lccdf(stay[n] | lambda[n]);\n")
    }else if(family=="gamma"){
      temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = gamma_lpdf(stay[n] | s, lambda[n]);\n")
      temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = gamma_lccdf(stay[n] | s, lambda[n]);\n")
    }else if(family=="lognormal"){
      temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = lognormal_lpdf(stay[n] | lambda[n], s);\n")
      temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = lognormal_lccdf(stay[n] | lambda[n], s);\n")
    }else if(family=="weibull"){
      temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = weibull_lpdf(stay[n] | lambda[n], s);\n")
      temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = weibull_lccdf(stay[n] | lambda[n], s);\n")
    }
    temp2 <- paste0(temp2,"\t\t}\n\t}\n")

    gq_code <- paste0(gq_code, temp2, "}\n")

    codestan <- paste(data_code,td_code,para_code,tp_code,model_code,gq_code,"\n")

    ###################
    ###running rstan###
    ###################

    rstan_options(auto_write=TRUE)
    options(mc.cores=parallel::detectCores())


    options(warn = -1)

    modelname <- paste0(formula_stay_paste, " [", family, "]")
    cat(paste(modelname, "\nNow compiling!\n"))
    stanmodel <- rstan::stan_model(model_name=modelname,model_code=codestan)
    cat("\nMCMC sampling start.\n")

    fitstan[[h]] <- rstan::sampling(stanmodel, data=datastan,
                                    iter=iter, warmup = warmup,chains=chains,thin=thin,cores=cores)


    ################################################
    ###calculating parameters###
    ############################

    ###calculating WAIC
    loglik_stay <- rstan::extract(fitstan[[h]],"log_lik_stay")$log_lik_stay
    lppd[h] <- sum(log(colMeans(exp(loglik_stay))))
    p_waic[h] <- sum(apply(loglik_stay,2,var))
    waic <- -lppd[h]/N_stay + p_waic[h]/N_stay
    waic2[h] <- waic * (2*N_stay)
  }
  result <- data.frame(Model = unlist(t),
                    Family = rep(family, length(waic2)),
                    lppd = round(lppd, 2),
                    p_waic = round(p_waic, 2),
                    WAIC = round(waic2, 2)
                    ) %>%
    arrange(WAIC)

  fit <- fitstan[[which(t== result$Model[1])]]
  mstay <- summary(fit)$summary[grep("mean_stay", rownames(summary(fit)$summary)), ]

  if(local_stay == FALSE) {
    mstay_temp <- rep(mstay, nrow(station_data))
    mstay_temp2 <- matrix(mstay_temp, byrow = TRUE, ncol = length(mstay))
    colnames(mstay_temp2) <- names(mstay)
    mstay <- mstay_temp2
  }
  stay_expected_local <- station_data %>%
    bind_cols(mstay)

  if(plot == TRUE){
    g <- ggplot(data = stay_expected_local) +
      geom_pointrange(aes(x = Station, y = mean,
                          ymin = `2.5%`, ymax = `97.5%`),
                      size = 0.75, linewidth = 0.75) +
      ylim(0, max(stay_expected_local$mean) * 2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))

    result <- list(result = result, plot = g)
    return(result)
  } else{
    return(result)
  }
}
WAIC <- '2.5%' <- '97.5%' <- NULL
