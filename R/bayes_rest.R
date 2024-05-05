#' Bayesian parameter estimation of REST/RAD-REST model based on MCMC samplings with Rstan
#'
#' @param formula_stay Model formula for the staying time within a focal area. Grammar follows lme4::glmer function. e.g. Stay ~ 1 + (1|Group)
#' @param formula_density Model formula for animal density. e.g. ~ 1 + x1
#' @param station_data_complete A data frame containing information for each camera station in each row
#' @param stay_data A data frame containing the time spent within a focal area
#' @param activity_data A vector of detection times transformed into radians
#' @param activity_model The method for estimating the proportion of active animals. If "kernel", the fixed kernel method by Rowcliffe et al. () is used. If "mixed", a von Mises mixture distribution is used. The default value is "kernel".
#' @param K The number of mixtures when activity_model = "mixed". The default value is 5.
#' @param bw_adj Bandwidth adjustment when activity_model = "mixed". The default value is 1.5. See Rowcliffe et al. () for details.
#' @param stay_family The probability distribution of the time spent. Specify the probability distribution selected by the bayes_stay_selection function.
#' @param local_stay Specify whether to use the expected value of the time spent for each camera station (local_stay = TRUE) or the global expected value (FALSE) for density estimation.
#' @param focal_area The area of a focal area
#' @param cores The number of cores to use for parallel computation. The default value is 1.
#' @param iter The length of iterations. The default value is 2000.
#' @param warmup The length of warmup. The default value is NULL (half of the iterations will be used for warmup).
#' @param chains The number of chains. The default value is 2.
#' @param thin The thinning interval. The default value is 1.
#' @param model Specify the model to use. Choose either "REST" or "RAD-REST".
#'
#' @return An object of class "stanfit" representing
#' the fitted Bayesian model. It contains samples from the posterior
#' distribution of the model parameters, which can be accessed using
#' the rstan extraction functions like rstan::extract. The object also
#' contains other useful quantities like the MCMC chains, sampling metadata etc.
#' @export
#' @import rstan activity parallel
#' @importFrom rstan extract
#' @importFrom stats as.formula formula model.frame model.matrix quantile sd var
#' @examples
#' station_data_rest <- format_station_data(
#'   detection_data = detection_data,
#'   station_data = station_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_y = "y",
#'   model = "REST",
#'   target_species = "A"
#' )
#'
#' station_data_complete <- add_effort(
#'   detection_data = detection_data,
#'   station_data_formatted = station_data_rest,
#'   col_name_station = "Station",
#'   col_name_term = "Term",
#'   col_name_datetime = "DateTime",
#'   plot = TRUE,
#'   font_size = 5
#' )
#'
#' stay_data <- format_stay(
#'   detection_data = detection_data,
#'   col_name_station = "Station",
#'   col_name_species = "Species",
#'   col_name_stay = "Stay",
#'   col_name_cens = "Cens",
#'   target_species = "A"
#' )
#'
#' activity_data <- format_activity(
#'   detection_data = detection_data,
#'   col_name_species = "Species",
#'   col_name_datetime = "DateTime",
#'   target_species = "A"
#' )
#'
#' fitstan <- bayes_rest(formula_stay = Stay ~ 1,
#'   formula_density = ~ 1,
#'   station_data_complete,
#'   stay_data,
#'   activity_data,
#'   activity_model = "kernel",
#'   local_stay = FALSE,
#'   stay_family = "exponential",
#'   focal_area = 2.0,
#'   model = "REST")
bayes_rest <- function(formula_stay,
                       formula_density,
                       station_data_complete,
                       stay_data,
                       activity_data,
                       activity_model = "kernel",
                       K = 5,
                       bw_adj = 1.5,
                       stay_family = "lognormal",
                       local_stay = FALSE,
                       focal_area,
                       cores = 2, iter = 2000, warmup = NULL, chains = 2, thin = 1,
                       model = "REST"){
  #######################
  ###define functions###
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

  #################
  ###option ?###
  #################
  center = TRUE
  cauchy = 2.5
  lkj_corr = 2
  #################
  ###check input###
  #################

  if(stay_family=="lognormal" || stay_family=="gamma" || stay_family=="weibull" || stay_family=="exponential"){

  }else{
    stop(paste0("Input stay_family type(",stay_family,") is incorrect."))
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

  station_data <- data.frame(station_data_complete)
  station_data <- station_data %>%
    filter(Effort != 0)

  #####################
  ###variable names for stay
  #####################

  formula_stay <- as.formula(formula_stay)
  formula_stay_paste <- paste0(formula_stay[2],formula_stay[1],formula_stay[3])
  xformula_stay <- nobars(formula_stay)
  if (inherits(xformula_stay, "name") & length(xformula_stay) == 1) {
    xformula_stay <- nobars(as.formula(paste(deparse(formula_stay), "+ 1")))
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
        # just intercept
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

  #####################
  ###variable names for density
  #####################
  # 密度にはランダム効果を入れない

  formula_density <- as.formula(formula_density)
  formula_density_paste <- paste0(formula_density[1],formula_density[2])
  xformula_density <- nobars(formula_density)
  if (inherits(xformula_density, "name") & length(xformula_density) == 1) {
    xformula_density <- nobars(as.formula(paste(deparse(formula_density), "+ 1")))
  }
  xname_density <- colnames(model.matrix(xformula_density, station_data))
  P_density <- length(xname_density)

  ######################
  ###creating dataset###
  ######################

  # 活動時間データ準備
  activity_data <- as.vector(unlist(activity_data))
  N_activity <- length(activity_data)

  # 努力量(0は削除済み)
  eff <- station_data[ , "Effort"]

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
  # # カテゴリカル共変量の処理（ランダム効果が入った場合はどうなるか要確認）
  # is_cat_stay <- sapply(model.frame(xformula_stay, dat2), is.character)
  # n_cat_stay <- sum(unlist(is_cat_stay))
  #
  # if(n_cat_stay > 0) {
  #   # # 展開前での処理
  #   # df_cat_stay <- data.frame(dat2[, names(which(is_cat_stay))])
  #   # tb_cat_stay <- list(0)
  #   # for(i in 1:n_cat_stay) {
  #   #   tb_cat_stay[[i]] <- table(df_cat_stay[ , i])
  #   # }
  #
  #   # 展開後の処理
  #   # binary_cols <- sapply(data.frame(x_stay_expand), function(col) all(col %in% c(0, 1)))
  #   # names(binary_cols) <- NULL
  #   # cat_covs_stay <- which(binary_cols)[-1]
  #   x_meanstay_covs <- sapply(data.frame(x_stay_expand), mean)
  # }
#
#   cat_stay_name <- c()
#   for(i in 1:n_cat_stay){
#     cat_stay_name[i] <- paste("cat_stay", i, sep = "")
#   }


  #　サイトレベルでの滞在時間の共変量
  station_data["Stay"] <- 0
  x_stay_station <- data.frame(model.matrix(xformula_stay, station_data))
  names(x_stay_station)[names(x_stay_station) == "X.Intercept."] <- "(Intercept)"

  # 応答変数（滞在時間）
  stay <- as.numeric(dat3[ , yname_stay])
  censored <- stay_data_join[ , "Cens"]
  x_stay <- subset(dat3, select = xname_stay)

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

  # 密度
  dat2 <- station_data
  dat2["(Intercept)"] <- 1
  if(center==TRUE) dat2 <- centering(dat2, xname_density)
  dat3 <- model.matrix(formula_density, dat2)
  x_density_expand <- model.matrix(xformula_density, dat2)

  x_meandensity_covs <- sapply(data.frame(x_density_expand), mean)
  # categorical_covs_dens <- sum(sapply(data.frame(x_station_temp), is.character)) > 0
  # # カテゴリカル共変量の処理（ランダム効果が入った場合はどうなるか要確認）
  # is_cat_density <- sapply(model.frame(xformula_density, dat2), is.character)
  # n_cat_density <- sum(unlist(is_cat_density))
  #
  # if(n_cat_density > 0) {
  #   # 展開前での処理
  #   df_cat_density <- data.frame(dat2[, names(which(is_cat_density))])
  #   tb_cat_density <- list(0)
  #   for(i in 1:n_cat_density) {
  #     tb_cat_density[[i]] <- table(df_cat_density[ , i])
  #   }
  #   # 展開後の処理
  #   binary_cols <- sapply(data.frame(x_density_expand), function(col) all(col %in% c(0, 1)))
  #   names(binary_cols) <- NULL
  #   cat_covs_density <- which(binary_cols)[-1]
  #   x_meandensity_covs <- sapply(data.frame(x_density_expand[, cat_covs_density]), mean)
  #   names(x_meandensity_covs) <- NULL
  #   x_meandensity_covs <- c(1, x_meandensity_covs)
  # }
  # cat_density_name <- c()
  # for(i in 1:n_cat_stay){
  #   cat_density_name[i] <- paste("cat_density", i, sep = "")
  # }

  # データの定義
  x_density <- subset(dat3, select = xname_density)
  N_station <- nrow(x_density)

  if(model == "REST"){
    Y <- station_data$Y
  } else {
    N <- station_data[, "N"]
    y <- station_data[, grep("y", colnames(station_data))]
    N_class <- ncol(y)
  }

  if(activity_model == "kernel"){
    model_active <- fitact(activity_data, bw = bwcalc(activity_data, K = 3), adj = 1.5, reps=1)
    actv <- model_active@act[1]
  }

  # datalistの作成
  if (R_stay > 0) {
    datastan <- list(N_stay = N_stay, P_stay = P_stay, R_stay = R_stay, Q_stay = Q_stay, G_stay = G_stay, stay = stay, x_stay = x_stay, censored = censored,
                     N_station = N_station, P_density = P_density, x_density = x_density,
                     activity_data = activity_data, N_activity = N_activity, K = K, eff = eff, S = focal_area)
  } else {
    datastan <- list(N_stay = N_stay, P_stay = P_stay, stay = stay, x_stay = x_stay, censored = censored,
                     N_station = N_station, P_density = P_density, x_density = x_density,
                     activity_data = activity_data, N_activity = N_activity, K = K, eff = eff, S = focal_area)
  }

  if(R_stay>0){
    for(i in 1:R_stay){
      datastan[z2name_stay[i]] <- list(z_stay[[i]])
      datastan[id2name_stay[i]] <- list(id_stay[[i]])
    }
  }

  if(activity_model == "kernel"){
    for(i in 1:R_stay){
      datastan["actv"] <- list(actv)
    }
  }

  if(model == "REST"){
    for(i in 1:R_stay){
      datastan["Y"] <- list(Y)
    }
  } else {
    datastan["N"] <- list(N)
    datastan["y"] <- list(y)
    datastan["N_class"] <- N_class
  }

  if(local_stay == TRUE){
    datastan["x_stay_station"] <- list(x_stay_station)
  }


  # if(n_cat_stay > 0){
  #   for(i in 1:n_cat_stay){
  #     datastan[cat_stay_name[i]] <- list(tb_cat_stay[[i]])
  #   }
  # }
  # if(n_cat_density > 0){
  #   for(i in 1:n_cat_density){
  #     datastan[cat_density_name[i]] <- list(tb_cat_density[[i]])
  #   }
  # }

  if(local_stay == FALSE){
    #atastan["cat_covs_stay"] <- list(c(1, cat_covs_stay))
    #datastan["N_cat_covs_stay"] <- list(length(cat_covs_stay) + 1)
    datastan["x_meanstay_covs"] <- list(x_meanstay_covs)
  }

  # if(n_cat_density > 0){
    # datastan["cat_covs_dens"] <- list(c(1, cat_covs_dens))
    # datastan["N_cat_covs_dens"] <- list(length(cat_covs_dens) + 1)
    datastan["x_meandensity_covs"] <- list(x_meandensity_covs)
  # }

  # library(cmdstanr)
  # model_0 <- cmdstan_model("rest_test.stan")
  # fit_0 <- model_0$sample(datastan, iter_warmup = 1000, iter_sampling = 2000,
  #                         thin = 10, chains = 3, parallel_chains = 3, refresh = 200, show_messages = T,
  #                         adapt_delta = 0.8)

  #rstan::stan(file = "./../RAD-REST/rest_test.stan", data = datastan, chains = 1 )

  #stanmodel <- rstan::stan(file = "rest_test.stan", data = datastan, chains = 1 )
  #summary(stanmodel, pars="global_density")$summary
  #summary(stanmodel, pars="beta")$summary
  #summary(stanmodel, pars="actv")$summary

  ##########################
  ####creating stan code####
  ##########################

  ###function
  fun_code <- ""
  if(activity_model == "mixed"){
    fun_code <- 'functions{\n'
    temp1 <- "\tvector seq_fun(real start, real end, int N_by){\n"

    temp1 <- paste0(temp1, "\t\treal h;\n")
    temp1 <- paste0(temp1, "\t\tvector[N_by] out;\n")
    temp1 <- paste0(temp1, "\t\th=(end-start)/(N_by-1);\n")

    temp1 <- paste0(temp1, "\t\tfor(i in 1:N_by){\n")
    temp1 <- paste0(temp1, "\t\t\tout[i]=start+(i-1)*h;\n")
    temp1 <- paste0(temp1, "\t\t}\n")
    temp1 <- paste0(temp1, "\t\treturn(out);\n")
    temp1 <- paste0(temp1, "\t}\n")

    temp2 <- "\tvector VMM_prod(vector Mu, vector kappa, vector eta, int K){\n"
    temp2 <- paste0(temp2, "\t\tvector[627] t = seq_fun(0.01, 2*pi()-0.01, 627);\n")
    temp2 <- paste0(temp2, "\t\tvector[627] f_val;\n")
    temp2 <- paste0(temp2, "\t\treal ps_[627, K];\n")
    temp2 <- paste0(temp2, "\t\tvector[K] log_eta = log(eta);\n")
    temp2 <- paste0(temp2, "\t\tfor(i in 1:627){\n")
    temp2 <- paste0(temp2, "\t\t\tfor(k in 1:K){\n")
    temp2 <- paste0(temp2, "\t\t\t\tps_[i,k] = log_eta[k] + von_mises_lpdf(t[i]|Mu[k],kappa[k]);\n")
    temp2 <- paste0(temp2, "\t\t\t}\n")
    temp2 <- paste0(temp2, "\t\t\tf_val[i] = log_sum_exp(ps_[i,:]);\n")
    temp2 <- paste0(temp2, "\t\t}\n")
    temp2 <- paste0(temp2, "\t\treturn(exp(f_val));\n")
    temp2 <- paste0(temp2, "\t}\n")

    temp3 <- "\treal activity_calc(vector f_val){\n"
    temp3 <- paste0(temp3, "\t\treal f_max;\n")
    temp3 <- paste0(temp3, "\t\treal activity;\n")
    temp3 <- paste0(temp3, "\t\tf_max = max(f_val);\n")
    temp3 <- paste0(temp3, "\t\tactivity = 1.0/(2*pi()*f_max);\n")
    temp3 <- paste0(temp3, "\t\treturn(activity);\n")
    temp3 <- paste0(temp3, "\t}\n")

    fun_code <- paste0(fun_code,temp1,temp2, temp3,"}\n")
    }
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

    # 活動時間
    if(activity_model == "mixed"){
      temp1 <- paste0(temp1, "\tint N_activity;\n")
      temp1 <- paste0(temp1, "\tvector<lower=0>[N_activity] activity_data;\n")
      temp1 <- paste0(temp1, "\tint<lower=1> K;\n")
    } else {
      temp1 <- paste0(temp1, "\treal actv;\n")
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
    temp3 <- paste0(temp3, "\tint<lower=1> P_density;\n")
    temp3 <- paste0(temp3, "\trow_vector[P_density] x_density[N_station];\n")
    temp3 <- paste0(temp3, "\treal<lower=0.0> S;\n")
    temp3 <- paste0(temp3, "\tvector<lower=0.0>[N_station] eff;\n")
    if(model == "REST"){
      temp3 <- paste0(temp3, "\tarray[N_station] int Y;\n")
    } else {
      temp3 <- paste0(temp3, "\tint <lower=0> N_class;\n")
      temp3 <- paste0(temp3, "\tarray[N_station] int N;\n")
      temp3 <- paste0(temp3, "\tarray[N_station, N_class] int y;\n")
    }
    if(local_stay == TRUE){
      temp3 <- paste0(temp3, "\trow_vector[P_stay] x_stay_station[N_station];\n")
    }

    if(local_stay == FALSE){
      if(length(x_meanstay_covs) > 1) temp3 <- paste0(temp3, "\trow_vector[P_stay] x_meanstay_covs;\n");
    }
    if(length(x_meandensity_covs) > 1) temp3 <- paste0(temp3, "\trow_vector[P_density] x_meandensity_covs;\n");


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
  temp1 <- paste0(temp1, "\tvector[P_density] alpha;\n")

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
  if(stay_family == "lognormal" || stay_family == "gamma" | stay_family == "weibull"){
    temp3 <- paste0("\t","real<lower=0> s;\n")
  }
  # REST-RAD
  if(model == "RAD-REST"){
    temp3 <- paste0(temp3, "\tvector[N_class-1] q;\n")
    temp3 <- paste0(temp3, "\tarray[N_station] vector[N_class-1] eps;\n")
    temp3 <- paste0(temp3, "\treal<lower=0.0> sigma;\n")
  }

  # 活動時間
  temp4 <- temp5 <- ''
  if(activity_model == "mixed"){
    temp4 <- paste0(temp4, "\tvector<lower=0.0, upper=2*pi()>[K] Mu;\n")
    temp4 <- paste0(temp4, "\tvector<lower=0.0, upper=50.0>[K] kappa;\n")
    temp4 <- paste0(temp4, "\tsimplex[K] eta;\n")
  }
  # // 要追加項目
  temp5 <- "\tvector<lower=0.0>[N_station] local_density;\n"
  temp5 <- paste0(temp5, "\treal<lower=0.0> theta;\n")


  para_code <- paste0(para_code,temp1,temp2,temp3,temp4,temp5,"}\n")

  ###transformed parameters
  tp_code <-'transformed parameters{\n'

  temp1 <- "\treal lambda[N_stay];\n"
  temp1 <- paste0(temp1, "\treal predict[N_stay];\n")

  if(stay_family == "lognormal" || stay_family == "gamma" | stay_family == "weibull"){
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
  if(stay_family == "gamma"){
    temp3 <- paste0(temp3,"\tscale = 1/s;\n")
  }else if(stay_family == "lognormal"){
    temp3 <- paste0(temp3,"\tscale = s^2;\n")
  }else if(stay_family == "weibull"){
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

  if(stay_family=="exponential"){
    temp5 <- paste0(temp5,"\t\tlambda[n] = 1 / exp(predict[n]);\n") # 左辺のpredictはパラメータ
  }else if(stay_family=="gamma"){
    temp5 <- paste0(temp5,"\t\tlambda[n] = s / exp(predict[n]);\n") # 分母が期待値
  }else if(stay_family=="lognormal"){
    temp5 <- paste0(temp5,"\t\tlambda[n] = log(exp(predict[n])) - (s*s/2);\n") # 分母が期待値
  }else if(stay_family=="weibull"){
    temp5 <- paste0(temp5,"\t\tlambda[n] = 1 / (log(s) - predict[n]);\n")
  }
  temp5 <-paste0(temp5,"\t}\n")

  # mean_stay                                                                     mean_stay
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
  temp6 <- ""
  if(model == "RAD-REST"){
    temp6 <- paste0(temp6, "\tarray[N_station] vector[N_class] w;\n")
    temp6 <- paste0(temp6, "\tvector[N_station] E_y;\n")
    temp6 <- paste0(temp6, "\tarray[N_station] vector[N_class] E_theta;\n")
    temp6 <- paste0(temp6, "\tvector[N_class] q_all;\n")
    temp6 <- paste0(temp6, "\tfor(j in 1:N_station) {\n")
    temp6 <- paste0(temp6, "\t\tvector[N_class] q_all_j;\n")
    temp6 <- paste0(temp6, "\t\tq_all_j = append_row(rep_vector(0, 1), q[1:(N_class-1)] + eps[j, 1:(N_class-1)]);\n")
    temp6 <- paste0(temp6, "\t\tE_theta[j] = exp(q_all_j)/sum(exp(q_all_j));\n")
    temp6 <- paste0(temp6, "\t\tE_y[j] = 0;\n")
    temp6 <- paste0(temp6, "\t\tfor(k in 1:N_class) {\n")
    temp6 <- paste0(temp6, "\t\t\tE_y[j] += (k-1) * E_theta[j, k];\n")
    temp6 <- paste0(temp6, "\t\t}\n")
    temp6 <- paste0(temp6, "\t}\n")
    temp6 <- paste0(temp6, "\tfor(j in 1:N_station) {\n")
    temp6 <- paste0(temp6, "\t\tw[j, 1] = 0.0;\n")
    temp6 <- paste0(temp6, "\t\tw[j, 2:N_class] = q[1:(N_class-1)] + eps[j, 1:(N_class-1)];\n")
    temp6 <- paste0(temp6, "\t}\n")
  }

  if(activity_model == "mixed"){
    temp6 <- paste0(temp6, "\treal<lower=0.0, upper=1.0> actv;\n")
    temp6 <- paste0(temp6, "\tvector<lower=0.0, upper=2*pi()>[627] f_val;\n")
    temp6 <- paste0(temp6, "\tf_val = VMM_prod(Mu, kappa, eta, K);\n")
    temp6 <- paste0(temp6, "\tactv = activity_calc(f_val);\n")
  }

  # // 要追加項目
  temp7 <- "\tvector<lower=0.0>[N_station] mu;\n"
  if(model == "REST") temp7 <- paste0(temp7, "\tmu = exp(log(S/10^6) + log(eff*24*60*60) + log(local_density) - log(mean_stay) + log(actv));\n")
  if(model == "RAD-REST") temp7 <- paste0(temp7, "\tmu = exp(log(S/10^6) + log(eff*24*60*60) + log(local_density) - log(mean_stay) + log(actv) - log(E_y));\n")


  tp_code <- paste0(tp_code,temp1,temp2,temp3,temp4,temp5,temp6,temp7,"}\n")


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
  if(stay_family=="lognormal" | stay_family=="gamma" | stay_family == "weibull"){ # これは自分で足した（Gammaについて）
    temp3 <- paste0(temp3,"\t","s ~ cauchy(0,",cauchy,");\n")
  }

  temp2 <-'\t'
  if(stay_family == "exponential"){
    temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
    temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
    temp2 <- paste0(temp2, "\t\ttarget += exponential_lpdf(stay[n] | lambda[n]);\n")
    temp2 <- paste0(temp2,"\t} else {\n")
    temp2 <- paste0(temp2, "\t\ttarget += exponential_lccdf(stay[n] | lambda[n]);\n")
  }else if(stay_family=="gamma"){
    temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
    temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
    temp2 <- paste0(temp2, "\t\ttarget += gamma_lpdf(stay[n] | s, lambda[n]);\n")
    temp2 <- paste0(temp2,"\t} else {\n")
    temp2 <- paste0(temp2, "\t\ttarget += gamma_lccdf(stay[n] | s, lambda[n]);\n")
  }else if(stay_family=="weibull"){
    temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
    temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
    temp2 <- paste0(temp2, "\t\ttarget += weibull_lpdf(stay[n] |lambda[n], s);\n")
    temp2 <- paste0(temp2,"\t} else {\n")
    temp2 <- paste0(temp2, "\t\ttarget += weibull_lccdf(stay[n] |lambda[n], s);\n")
  }else if(stay_family=="lognormal"){
    temp2 <- paste0(temp2, "for(n in 1:N_stay){\n")
    temp2 <- paste0(temp2, "\tif(censored[n] == 0) {\n")
    temp2 <- paste0(temp2, "\t\ttarget += lognormal_lpdf(stay[n] |lambda[n], s);\n")
    temp2 <- paste0(temp2, "\t} else {\n")
    temp2 <- paste0(temp2, "\t\ttarget += lognormal_lccdf(stay[n] |lambda[n], s);\n")
  }
  temp2 <- paste0(temp2, "\t\t}\n\t}\n")

  # 活動時間
  temp1 <- ""
  if(activity_model == "mixed"){
  temp1 <- "\tvector[K] ps;\n"
  temp1 <- paste(temp1, "\tfor(n in 1:N_activity){\n\t\tfor(k in 1:K){\n")
  temp1 <- paste(temp1, "\t\t\tps[k] = log(eta[k]) + von_mises_lpdf(activity_data[n] | Mu[k], kappa[k]);\n")

  temp1 <- paste(temp1, "\t\t}\n\t\ttarget += log_sum_exp(ps);\n\t}\n")
  temp1 <- paste(temp1, "\ttarget += uniform_lpdf(Mu|0.0,2*pi());\n")
  temp1 <- paste(temp1, "\ttarget += gamma_lpdf(kappa|0.1,0.1);\n")
  }
  # // 要追加項目
  temp4 <- "\ttarget += cauchy_lpdf(theta | 0.0 , 20.0);\n"
  if(model == "REST"){
    temp4 <- paste0(temp4, "\ttarget += poisson_lpmf(Y | mu);\n")
  } else {
    temp4 <- paste0(temp4, "\ttarget += poisson_lpmf(N | mu);\n")
  }

  temp4 <- paste0(temp4, "\tfor(i in 1:N_station){\n")
  temp4 <- paste0(temp4, "\t\ttarget += gamma_lpdf(local_density[i] | theta, theta / exp(x_density[i]*alpha));\n")
  temp4 <- paste0(temp4, "\t}\n")

  if(model == "RAD-REST"){
    temp4 <- paste0(temp4, "\tfor(j in 1:N_station){\n")
    temp4 <- paste0(temp4, "\t\ttarget += normal_lpdf(eps[j, :] | 0.0, sigma);\n")
    temp4 <- paste0(temp4, "\t\ttarget += multinomial_logit_lpmf(y[j, :] | w[j, :]);\n")
    temp4 <- paste0(temp4, "\t}\n")
  }


  model_code <- paste0(model_code,temp3,temp2,temp1,temp4,"}\n")

  ###generated quantities

  gq_code <-'generated quantities{\n'

  # // 要追加項目
  temp0 <- "\treal<lower=0.0> expected_global_density;\n"
  if(length(x_meandensity_covs) > 1) temp0 <- paste0(temp0, "\texpected_global_density = exp(x_meandensity_covs * alpha);\n")
  if(length(x_meandensity_covs) == 1) temp0 <- paste0(temp0, "\texpected_global_density = exp(alpha[1]);\n")

  if(P_density > 1){
    temp0 <- paste0(temp0, "\tvector<lower=0.0>[N_station] expected_local_density;\n")
    temp0 <- paste0(temp0, "\tfor(i in 1:N_station){\n")
    temp0 <- paste0(temp0, "\t\texpected_local_density[i] = exp(x_density[i]*alpha);\n")
    temp0 <- paste0(temp0, "\t}\n")
  }

  temp1 <- ''
  if(activity_model == "mixed"){
  temp1 <- paste0(temp1, "\tvector[K] ps;\n")
  temp1 <- paste0(temp1, "\tvector[N_activity] loglik_time;\n")
  temp1 <- paste0(temp1, "\tfor(n in 1:N_activity){\n")
  temp1 <- paste0(temp1, "\t\tfor(k in 1:K){\n")
  temp1 <- paste0(temp1, "\t\t\tps[k] = log(eta[k]) + von_mises_lpdf(activity_data[n] | Mu[k], kappa[k]);\n")
  temp1 <- paste0(temp1, "\t\t}\n")
  temp1 <- paste0(temp1, "\t\tloglik_time[n] = log_sum_exp(ps);\n")
  temp1 <- paste0(temp1, "\t}\n")
  }

  if(model == "RAD-REST"){
    temp1 <- paste0(temp1, "\tvector[N_station] loglik_ey;\n")
    temp1 <- paste0(temp1, "\tfor(j in 1:N_station) {\n")
    temp1 <- paste0(temp1, "\t\t\tloglik_ey[j] = multinomial_logit_lpmf(y[j, :] | w[j, :]);\n")
    temp1 <- paste0(temp1, "\t}\n")
  }

  temp2 <- ''
  temp2 <- paste0(temp2, "\n\treal log_lik_stay[N_stay];\n")
  temp2 <- paste0(temp2, '\tfor(n in 1:N_stay){\n')

  if(stay_family == "exponential"){
    temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = exponential_lpdf(stay[n] | lambda[n]);\n")
    temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = exponential_lccdf(stay[n] | lambda[n]);\n")
  }else if(stay_family=="gamma"){
    temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = gamma_lpdf(stay[n] | s, lambda[n]);\n")
    temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = gamma_lccdf(stay[n] | s, lambda[n]);\n")
  }else if(stay_family=="lognormal"){
    temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = lognormal_lpdf(stay[n] | lambda[n], s);\n")
    temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = lognormal_lccdf(stay[n] | lambda[n], s);\n")
  }else if(stay_family=="weibull"){
    temp2 <- paste0(temp2,"\tif(censored[n] == 0) {\n\t\tlog_lik_stay[n] = weibull_lpdf(stay[n] | lambda[n], s);\n")
    temp2 <- paste0(temp2,"\t} else {\n\t\tlog_lik_stay[n] = weibull_lccdf(stay[n] | lambda[n], s);\n")
  }
  temp2 <- paste0(temp2,"\t\t}\n\t}\n")

  temp3 <- ""
  temp3 <- paste0(temp3, "\treal log_lik_station[N_station];\n")
  temp3 <- paste0(temp3, "\tfor(i in 1:N_station){\n")
  if(model == "REST"){
    temp3 <- paste0(temp3, "\t\tlog_lik_station[i] = poisson_lpmf(Y[i] | mu[i]);\n")
  } else {
    temp3 <- paste0(temp3, "\t\tlog_lik_station[i] = poisson_lpmf(N[i] | mu[i]);\n")
  }
  temp3 <- paste0(temp3, "\t}\n")

  gq_code <- paste0(gq_code, temp0, temp1, temp2, temp3, "}\n")

  codestan <- paste(fun_code, data_code,td_code,para_code,tp_code,model_code,gq_code,"\n")

  ###################
  ###running rstan###
  ###################

  rstan_options(auto_write=TRUE)

  options(warn = -1)
  modelname <- paste0(formula_stay_paste, " density", formula_density_paste," [",stay_family,"]")
  cat(paste(modelname, "\nNow compiling!\n"))
  stanmodel <- rstan::stan_model(model_name=modelname,model_code=codestan)
  cat("\nMCMC sampling start.\n")


  fitstan <- rstan::sampling(stanmodel, data=datastan,
                             iter=iter, warmup = warmup,chains=chains,thin=thin,cores=cores)

  ############################
  ###calculating parameters###
  ############################

  ###calculating WAIC
  loglik_stay <- rstan::extract(fitstan,"log_lik_stay")$log_lik_stay
  loglik_station <- rstan::extract(fitstan,"log_lik_station")$log_lik_station
  loglik <- cbind(loglik_stay, loglik_station)
  if(activity_model == "mixed"){
    loglik_activity <- rstan::extract(fitstan,"loglik_activity")$loglik_activity
    loglik <- cbind(loglik_stay, loglik_activity)
  }

  lppd <- sum(log(colMeans(exp(loglik))))
  p_waic <- sum(apply(loglik,2,var))
  waic <- -lppd/N_stay + p_waic/N_stay
  waic2 <- waic * (2*N_stay)

  attr(fitstan,"WAIC") <-c("WAIC" = waic2, "lppd" = lppd, "p_waic" = p_waic)

  attr(fitstan,"formula_stay") <- c("formula" = formula_stay)
  attr(fitstan,"formula_density") <- c("formula" = formula_density)
  if(P_density > 1) res_local_dens <- round(summary(fitstan)$summary[grep("expected_local_density", rownames(summary(fitstan)$summary)), ], 2)
  res_global_dens <- round(summary(fitstan)$summary[grep("expected_global_density", rownames(summary(fitstan)$summary)), ], 2)
  mean_stay <- round(summary(fitstan)$summary[grep("mean_stay", rownames(summary(fitstan)$summary)), ], 2)

  if(P_density > 1) attr(fitstan,"expected_local_density") <- res_local_dens
  attr(fitstan,"expected_global_density") <- res_global_dens
  attr(fitstan,"mean_stay") <- mean_stay

  attr(fitstan,"model") <- stanmodel

  outputwaic <- paste0("\nlppd = ",round(lppd,digits = 4),"\npWAIC = ",round(p_waic,digits=4),
                       "\nWAIC = ", round(waic2,digits = 4))
  cat(outputwaic)

  return(fitstan)
}
Effort <- NULL
