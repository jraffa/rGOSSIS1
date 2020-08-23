#' Apply simple preprocessing of a data.frame in preparation for GOSSIS imputation
#'
#' This function applies some simple preprocessing to a data.frame that is to undergo imputation.
#' Most operations inolve binning categories or handling missing diagnosis data.
#' A full description of transformations can be found in the code comments and below.
#'
#' @param dat a data.frame with the minimal GOSSIS-1 dataset, and correct column names.
#' @return a data.frame that has been processed and ready to begin imputation.
#' @author Jesse D. Raffa
#' @details
#' The function joins the ap_dx data.frame to dat, mapping the diagnoses into group, dx_class, and dx_sub.
#' It also handles missing or mismatched diagnoses.
#' The function also maps the 3 GCS score components into a single 4 level factor variable.
#' @seealso ap_dx
#' @export
#' @importFrom stringr str_split
#' @import dplyr
#' @examples
#' require(GOSSIS)
#' # new_dat <- preprocess_data(dat)
#'
preprocess_data <- function(dat) {
  # Parse APACHE Dx List dividing up the codes by the . e.g., 101.01 = ['101', '01'] for [class, subdx]
  ap_dx$sub_Dx <-
    sapply(str_split(as.character(ap_dx$code), "\\."), length)
  ap_dx$dx_class <-
    sapply(str_split(as.character(ap_dx$code), "\\."), "[[", 1)
  ap_dx$dx_sub <-
    sapply(str_split(as.character(ap_dx$code), "\\."), function(x) {
      return(ifelse(is.na(x[2]), "99", x[2]))
    })

  dat$sub_Dx <-
    sapply(str_split(as.character(dat$apache_3j_diagnosis), "\\."), length)
    dat$dx_class <-
    sapply(str_split(as.character(dat$apache_3j_diagnosis), "\\."), "[[", 1)
  dat$dx_sub <-
    sapply(str_split(as.character(dat$apache_3j_diagnosis), "\\."), function(x) {
      return(ifelse(is.na(x[2]), "99", x[2]))
    })
  # This is for use in a previous version of the model
  dat$d1_creatinine_range <-
    dat$d1_creatinine_max - dat$d1_creatinine_min
  # GCS is put into 4 bins by looking at all 3 score components
  dat$dcs_all_high <-
    as.numeric(dat$gcs_eyes_apache == 4 &
                 dat$gcs_verbal_apache == 5 &
                 dat$gcs_motor_apache == 6)
  dat$dcs_some_high <-
    as.numeric(
      (
        dat$gcs_eyes_apache == 4 |
          dat$gcs_verbal_apache == 5 |
          dat$gcs_motor_apache == 6
      ) & dat$dcs_all_high == 0
    )
  dat$dcs_all_low <-
    as.numeric(dat$gcs_eyes_apache == 1 &
                 dat$gcs_verbal_apache == 1 &
                 dat$gcs_motor_apache == 1)
  dat$dcs_intermed <-
    as.numeric(!with(dat, dcs_all_high | dcs_some_high | dcs_all_low))
  dat$vent <- as.numeric(dat$ventilated_apache | dat$intubated_apache)
  #dat %>% mutate(dx_class)
  dat <-
    dat %>% mutate(
      dcs_group = case_when(
        dcs_all_high == 1 ~ "ALL_HIGH",
        dcs_some_high == 1 ~ "SOME_HIGH",
        dcs_all_low == 1 ~ "ALL_LOW",
        dcs_intermed == 1 ~ "INTERM_DCS"
      ),
      dcs_group = as.factor(dcs_group)
    ) %>%
    select(-dcs_all_high,-dcs_some_high,-dcs_all_low,-dcs_intermed)
  # Remove those without an outcome
  dat <- dat %>% filter(!is.na(hospital_death))
  # We call bodysystem group simply 'group'
  # eicu has blank groups which we put in the "Other medical disorder category"
  dat <-
    dat %>% rename(group = apache_3j_bodysystem) %>% mutate(
      group = as.factor(
        ifelse(
          data_source == "eicu" &
            group == "",
          "Other medical disorders",
          as.character(group)
        )
      ),
      dx_class = ifelse(data_source == "eicu" &
                          dx_class == 0, 1002, dx_class)
    )
  # Make the dx groups factors.
  dat <-
    dat %>% mutate(
      group = as.factor(group),
      data_source = as.factor(data_source),
      icu_type = as.factor(icu_type),
      icu_admit_source = as.factor(icu_admit_source)
    )
  # Drop old levels like ""
  dat <- dat %>%  mutate( group = droplevels(group))
  # Map missing or unknown dx classes to "99" which is "Unknown" or unclassified
  dat <-
    dat %>% mutate(dx_class = as.factor(ifelse(is.na(dx_class), "99", dx_class)), dx_sub =
                        as.factor(dx_sub))
  return(dat)
}

check_data <- function(x) {
  return(TRUE)
}

#' Apply imputation procedure for the Global Open Source Severity of Illness Scale.
#'
#' This function imputes missing data based on models developed for the Global Open Source Severity of Illness Scale.
#' It offers several variants of imputation, including median imputation, with the suggested implementation, algorithm 3
#'
#' @param x a data.frame with the minimal GOSSIS-1 dataset, and correct column names, that has been preprocessed with preprocess_data
#' @param algo the imputation algorithm selected to use.  defaults to option 3.
#' Option 1 is a model where all predictors (except the variable to be imputed and the dx_sub variable) are used in prediction.
#' This algorithm is applied to all missing data.
#' Option 2 is a model which uses all predictors except the variable to be imputed, the variables related to the variable being imputated (e.g., the other _min/_max and _apache variables of the same variable).
#' The algorithm is applied to all missing data.
#' Option 3 is a hybrid of algorithm 1 and 2, applying algorithm 1 to those observations which have a similar variable (e.g., someone with _apache, but missing _max), and algorithm 2 for those missing all of _min, _max and _apache variables of the same type.
#' Option 4 imputes the median or mode of each variable using the training dataset.
#' @param files optional argument to specify the files which to use for imputation.  By default it uses the ones used in the paper.
#' @return a data.frame that has been completed
#' @author Jesse D. Raffa
#' @details
#' The function applies a series of prediction models used for imputation of missing data contained in a GOSSIS dataset.
#' There are currently 4 variants of the algorithm implemented.  All four return a completed (missing data filled in) dataset.
#' Please use this function carefully.  It is not recommended for non-GOSSIS applications.
#' @export
#' @import dplyr
#' @import caret
#' @import stringr
#' @importFrom stats na.pass
#' @examples
#' require(GOSSIS)
#' # new_dat <- preprocess_data(dat)
#' # imp_dat <- impute_data(dat,algo=3)
#'

impute_data <- function(x,algo=3,files=NULL,median_Overide=FALSE,mediandata=NULL) {
  all_vars <- x %>% select(
    starts_with("d1"),
    ends_with("apache"),
    contains("dcs"),-contains("invasive"),-contains("prob"),-contains("death"),-contains("score"),-contains("gcs"),
    group,-d1_creatinine_min,
    icu_admit_source,
    icu_type,
    elective_surgery,
    age,
    dx_class,
    dx_sub,
    vent,
    aids,
    cirrhosis,
    diabetes_mellitus,
    lymphoma,
    solid_tumor_with_metastasis,-contains("urine"),
    hepatic_failure
  ) %>% names()
  if(is.null(files)) {
	imputation_models_ <- c(imputation_models,imputation_models_fac)
  } else if(sum(file.exists(files))==2) {
        load(files[1])
	load(files[2])
	imputation_models_ <- c(o1,o2)
	if(class(imputation_models_)!="list") {
		stop("imputation files not loaded properly")
	}
  }
  vars <- sapply(imputation_models_,"[[","variable")
  xx <- x
  if(algo==1) {  # uses algorithm type I from the imputation models object.
    for(v in vars) {
      idx <- which(v==vars)[1]
      if(sum(is.na(xx[[v]]))>0) {
      xx[[v]][is.na(xx[[v]])] <-
        predict.train(imputation_models_[[idx]]$typeI,newdata=x %>% filter(is.na(x[[v]])), na.action = na.pass)
      }
    }
  }
  if(algo==2) { # uses algorithm type II from the imputation models object.
    for(v in vars) {
      idx <- which(v==vars)[1]
      if(is.null(imputation_models_[[idx]]$typeII) ) {
        if(sum(is.na(xx[[v]]))>0) {
        xx[[v]][is.na(xx[[v]])] <-
          predict.train(imputation_models_[[idx]]$typeI,newdata=x %>% filter(is.na(x[[v]])), na.action = na.pass)
        }
      } else {
        if(sum(is.na(xx[[v]]))>0) {
        xx[[v]][is.na(xx[[v]])] <-
          predict.train(imputation_models_[[idx]]$typeII,newdata=x %>% filter(is.na(x[[v]])), na.action = na.pass)
        }
      }
    }
  }
  if(algo==3) {  # uses a hybrid of ty  pe I and type II depending on what related variables are present.
    for(v in vars) {
      idx <- which(v==vars)[1]
      if(grepl("d1_",v)) {
        root_var <- str_replace(v, "d1_", "")
        root_var <- str_replace(root_var, "_max|_min|_range", "")
        ar_2 <- grep(root_var,all_vars, value = TRUE)
        #print(ar_2)
        if (length(ar_2) == 2) {
          #print(root_var)
          missing_ <-
            which(is.na(x %>% select((v))))
          missing_2 <-(is.na(x %>% select((v))) &
                         x %>% mutate_at(ar_2, is.na) %>% select_at(ar_2) %>% rowSums ==2)
          #print(length(missing_))
          #print(sum(missing_2))
          if(sum(is.na(x[[v]][missing_]))>0) {
            #print("I")
            #print(sum(is.na(xx[[v]][missing_])))
          xx[[v]][missing_] <-
            predict.train(imputation_models_[[idx]]$typeI,newdata=x[missing_,], na.action = na.pass)
          }
          if(sum(is.na(x[[v]][missing_2]))>0) {
            #print(sum(is.na(xx[[v]][missing_2])))
            #print("II")
          xx[[v]][missing_2] <-
            predict.train(imputation_models_[[idx]]$typeII,newdata=x[missing_2,], na.action = na.pass)
          }
        } else if(length(ar_2) == 3) {
          #print(root_var);
          missing_ <-
            which(is.na(x %>% select((v))))
          missing_2 <-
            (
              is.na(x  %>% select((v))) &
                x %>% mutate_at(ar_2, is.na) %>% select_at(ar_2) %>% rowSums ==
                3
            )
          #print(length(missing_))
          #print(sum(missing_2))
          if(sum(is.na(x[[v]][missing_]))>0) {
            #print(sum(is.na(x[[v]][missing_])))
            #print("III")
          xx[[v]][missing_] <-
            predict.train(imputation_models_[[idx]]$typeI,newdata=x[missing_,], na.action = na.pass)
          }
          if(sum(is.na(x[[v]][missing_2]))>0) {
            #print(sum(is.na(x[[v]][missing_2])))
            #print("IV")
          xx[[v]][missing_2] <-
            predict.train(imputation_models_[[idx]]$typeII,newdata=x[missing_2,], na.action = na.pass)
          }
        }

      } else {
        if(sum(is.na(x[[v]]))>0) {
          #print("V")
        xx[[v]][is.na(x[[v]])] <-
          predict.train(imputation_models_[[idx]]$typeI,newdata=x %>% filter(is.na(x[[v]])), na.action = na.pass)
        }
      }
      }
  }
  if(algo==4) {  # Uses the median/mode values computed in the training set.
    if(median_Overide &is.null(mediandata)) {
      get_mode_dat <- function(x) {
        return(names(which.max(table(x))) )
      }
      medians_numeric <- x %>% summarise_at(vars(starts_with("d1_"),age,),funs(median),na.rm=T)
      modes_cat <- x %>% select_at(vars(fac_vars)) %>% summarise_all(funs(get_mode_dat))
    } else if(!is.null(mediandata)) {
      medians_numeric <- mediandata[[1]]
      modes_cat <- mediandata[[2]]
    }
    for(v in vars) {
      if(v %in% names(medians_numeric)) {
        xx[[v]][is.na(xx[[v]])] <- medians_numeric[[v]]
      } else if(v %in% names(modes_cat)) {
        xx[[v]][is.na(xx[[v]])] <- modes_cat[[v]]
      }
    }
  }
  return(xx)
}

load_imputation_model <- function(x) {

}

make_predictions <- function(x,type="hosp_mortality") {
  preds <- NULL
  return(preds)
}


evaluate_prediction <- function(y,pred,strata=NULL) {
  auc <- NULL
  SMR <- NULL
  return(list(auc=auc,SMR=SMR))
}


sums_and_diffs <- function(x) {
  x <- x %>% mutate(d1_diasbp_avg = (d1_diasbp_max + d1_diasbp_min)/2,
                    d1_diasbp_diff = (d1_diasbp_max - d1_diasbp_min),
                    d1_heartrate_avg = (d1_heartrate_max + d1_heartrate_min)/2,
                    d1_heartrate_diff = (d1_heartrate_max - d1_heartrate_min),
                    d1_mbp_avg = (d1_mbp_max + d1_mbp_min)/2,
                    d1_mbp_diff = (d1_mbp_max - d1_mbp_min),
                    d1_resprate_avg = (d1_resprate_max + d1_resprate_min)/2,
                    d1_resprate_diff = (d1_resprate_max - d1_resprate_min),
                    d1_spo2_avg = (d1_spo2_max + d1_spo2_min)/2,
                    d1_spo2_diff = (d1_spo2_max - d1_spo2_min),
                    d1_sysbp_avg = (d1_sysbp_max + d1_sysbp_min)/2,
                    d1_sysbp_diff = (d1_sysbp_max - d1_sysbp_min),
                    d1_temp_avg = (d1_temp_max + d1_temp_min)/2,
                    d1_temp_diff = (d1_temp_max - d1_temp_min),
                    d1_albumin_avg = (d1_albumin_max + d1_albumin_min)/2,
                    d1_albumin_diff = (d1_albumin_max - d1_albumin_min),
                    d1_bilirubin_avg = (d1_bilirubin_max + d1_bilirubin_min)/2,
                    d1_bilirubin_diff = (d1_bilirubin_max - d1_bilirubin_min),
                    d1_bun_avg = (d1_bun_max + d1_bun_min)/2,
                    d1_bun_diff = (d1_bun_max - d1_bun_min),
                    d1_creatinine_avg = (d1_creatinine_max + d1_creatinine_max - d1_creatinine_range)/2,
                    d1_creatinine_diff = (d1_creatinine_range),
                    d1_glucose_avg = (d1_glucose_max + d1_glucose_min)/2,
                    d1_glucose_diff = (d1_glucose_max - d1_glucose_min),
                    d1_hco3_avg = (d1_hco3_max + d1_hco3_min)/2,
                    d1_hco3_diff = (d1_hco3_max - d1_hco3_min),
                    d1_hemaglobin_avg = (d1_hemaglobin_max + d1_hemaglobin_min)/2,
                    d1_hemaglobin_diff = (d1_hemaglobin_max - d1_hemaglobin_min),
                    d1_hematocrit_avg = (d1_hematocrit_max + d1_hematocrit_min)/2,
                    d1_hematocrit_diff = (d1_hematocrit_max - d1_hematocrit_min),
                    d1_inr_avg = (d1_inr_max + d1_inr_min)/2,
                    d1_inr_diff = (d1_inr_max - d1_inr_min),
                    d1_lactate_avg = (d1_lactate_max + d1_lactate_min)/2,
                    d1_lactate_diff = (d1_lactate_max - d1_lactate_min),
                    d1_platelets_avg = (d1_platelets_max + d1_platelets_min)/2,
                    d1_platelets_diff = (d1_platelets_max - d1_platelets_min),
                    d1_potassium_avg = (d1_potassium_max + d1_potassium_min)/2,
                    d1_potassium_diff = (d1_potassium_max - d1_potassium_min),
                    d1_sodium_avg = (d1_sodium_max + d1_sodium_min)/2,
                    d1_sodium_diff = (d1_sodium_max - d1_sodium_min),
                    d1_calcium_avg = (d1_calcium_max + d1_calcium_min)/2,
                    d1_calcium_diff = (d1_calcium_max - d1_calcium_min),
                    d1_wbc_avg = (d1_wbc_max + d1_wbc_min)/2,
                    d1_wbc_diff = (d1_wbc_max - d1_wbc_min),
                    d1_arterial_pco2_avg = (d1_arterial_pco2_max + d1_arterial_pco2_min)/2,
                    d1_arterial_pco2_diff = (d1_arterial_pco2_max - d1_arterial_pco2_min),
                    d1_arterial_ph_avg = (d1_arterial_ph_max + d1_arterial_ph_min)/2,
                    d1_arterial_ph_diff = (d1_arterial_ph_max - d1_arterial_ph_min),
                    d1_arterial_po2_avg = (d1_arterial_po2_max + d1_arterial_po2_min)/2,
                    d1_arterial_po2_diff = (d1_arterial_po2_max - d1_arterial_po2_min),
                    d1_pao2fio2ratio_avg = (d1_pao2fio2ratio_max + d1_pao2fio2ratio_min)/2,
                    d1_pao2fio2ratio_diff = (d1_pao2fio2ratio_max - d1_pao2fio2ratio_min))
  return(x)
}

#' Prepare data frame for model fitting in the Global Open Source Severity of Illness Scale.
#'
#' This function prepares the data for fitting the prediction models
#'
#' @param x a data.frame with the minimal GOSSIS-1 dataset, and correct column names, that has been preprocessed with preprocess_data
#' @return a data.frame that has been prepared
#' @author Jesse D. Raffa
#' @details
#' The function applies a series of prediction models used for imputation of missing data contained in a GOSSIS dataset.
#' There are currently 4 variants of the algorithm implemented.  All four return a completed (missing data filled in) dataset.
#' Please use this function carefully.  It is not recommended for non-GOSSIS applications.
#' @export
#' @import dplyr
#' @import stringr
#' @examples
#' require(GOSSIS)
#' # new_dat <- preprocess_data(dat)
#' # imp_dat <- impute_data(dat,algo=3)
#' # dat <- prepare_fit(imp_dat)
#'


prepare_fit <- function(x) {
  x <- sums_and_diffs(x)
  x <- x %>% mutate(dx_class = as.factor(paste0(dx_class,group)),dx_sub = (paste0(as.character(dx_sub),dx_class,group)))  %>%
    mutate(elective_surgery=as.factor(as.numeric(elective_surgery>0)),
           vent = as.factor(as.numeric(vent>0)),
           lymphoma = as.factor(as.numeric(lymphoma>0)),
           hepatic_failure = as.factor(as.numeric(hepatic_failure>0)),
           solid_tumor_with_metastasis = as.factor(as.numeric(solid_tumor_with_metastasis>0)),
           diabetes_mellitus = as.factor(as.numeric(diabetes_mellitus>0)),
           cirrhosis = as.factor(as.numeric(cirrhosis>0)),
           aids = as.factor(as.numeric(aids>0)),
           arf_apache=as.factor(as.numeric(arf_apache>0)))

  levs_in <- complete_levels[complete_levels %in% unique(as.character(x$dx_sub))]
  levs_out <- complete_levels[!(complete_levels %in% unique(as.character(x$dx_sub)))]
  x$dx_sub <- factor(x$dx_sub,levels=levs_in)
  if(length(levs_out)>0) {
    levels(x$dx_sub)[(length(levs_in)+1):length(c(levs_in,levs_out))] <- levs_out
  }
  #drop_facs <- complete_levels[!(complete_levels %in% levels(mod$model$dx_sub))]
  x$y <- x$hospital_death
  return(x)
}


#' Generate mgcv model forumula to fit the Global Open Source Severity of Illness Scale.
#'
#' This function generates a the model formula to be fit with gam/bam
#'
#' @param basis "\"cr\"" or "\"cr\"" string to be used in fitting
#' @param kay the degrees of freedom for the fit
#' @return a string that can be used as a formula.
#' @author Jesse D. Raffa
#' @details
#' Helper function to generate a model formula.
#' Used mainly in tuning the model during development.
#' @export
#' @examples
#' require(GOSSIS)
#' # new_dat <- preprocess_data(dat)
#' # imp_dat <- impute_data(dat,algo=3)
#' # dat <- prepare_fit(imp_dat)
#' # here
#' # fit  <- gam(formula(formm),data=dp_io2,family="binomial",na.action=na.exclude)
#'

generate_formula <- function(basis,kay) {
  formm <-  paste0("y ~ te(d1_diasbp_avg,d1_diasbp_diff, bs = ", basis, ", k = ", kay, ") + te(d1_heartrate_avg,d1_heartrate_diff, bs = ", basis, ", k = ", kay, ") +
te(d1_mbp_diff, d1_mbp_avg,bs = ", basis, ", k = ", kay, ") + te(d1_resprate_diff,d1_resprate_avg, by=vent,bs = ", basis, ", k = ", kay, ") +
    te(d1_spo2_avg,d1_spo2_diff,bs = ", basis, ", k = ", kay, ") + te(d1_sysbp_avg,d1_sysbp_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_temp_avg,d1_temp_diff,bs = ", basis, ", k = ", kay, ") + te(d1_albumin_avg,d1_albumin_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_bilirubin_avg,d1_bilirubin_diff,bs = ", basis, ", k = ", kay, ") + te(d1_bun_avg,d1_bun_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_creatinine_avg,d1_creatinine_diff,bs = ", basis, ", k = ", kay, ",by=arf_apache) +
    te(d1_glucose_diff,d1_glucose_avg,bs = ", basis, ", k = ", kay, ") + te(d1_hco3_avg,d1_hco3_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_hemaglobin_avg,d1_hemaglobin_diff,bs = ", basis, ", k = ", kay, ") + te(d1_hematocrit_avg,d1_hematocrit_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_inr_avg,d1_inr_diff,bs = ", basis, ", k = ", kay, ") + te(d1_lactate_avg,d1_lactate_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_platelets_avg,d1_platelets_diff,bs = ", basis, ", k = ", kay, ") + te(d1_potassium_avg,d1_potassium_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_sodium_avg,d1_sodium_diff,bs = ", basis, ", k = ", kay, ")  + te(d1_wbc_avg,d1_wbc_diff,bs = ", basis, ", k = ", kay, ") +
    te(d1_arterial_pco2_avg,d1_arterial_pco2_diff,bs = ", basis, ", k = ", kay, ",by=vent)  +
    te(d1_arterial_ph_avg,d1_arterial_ph_diff,bs = ", basis, ", k = ", kay, ",by=vent) +
    te(d1_arterial_po2_avg,d1_arterial_po2_diff,bs = ", basis, ", k = ", kay, ",by=vent) +
    te(d1_pao2fio2ratio_avg,d1_pao2fio2ratio_diff,bs = ", basis, ", k = ", kay, ",by=vent) +
    te(d1_calcium_avg,d1_calcium_diff,bs = ", basis, ", k = ", kay, ") +
    elective_surgery + s(age,bs = ", basis, ", k = ", kay, ")  + dcs_group +
  icu_admit_source + aids + cirrhosis +  diabetes_mellitus + lymphoma +  solid_tumor_with_metastasis +   s(group,bs='re') + s(dx_class,bs='re') + s(dx_sub,bs='re') ")
  return(formm)
}

#' @export
#' @import mgcv
#'
gpredict <- function(mod,newdata,na.action="pass0") {
  mod_dx <- gsub("dx_sub","",grep("dx_sub",names(mod$cmX),value = TRUE))
  out_lev <- levels(newdata$dx_sub)[!(levels(newdata$dx_sub) %in% mod_dx)]
  out <- rep(NA,nrow(newdata))
  out[!(newdata$dx_sub %in% out_lev)] <- plogis(predict(mod,newdata=newdata %>% filter(!(dx_sub %in%out_lev)) ))
  if(sum(newdata$dx_sub %in%out_lev)) {
    message("detected new dx_sub levels; imputing dx_class mean")
  try( { tmp <- predict(mod,newdata=newdata %>% filter((dx_sub %in%out_lev))%>%
                          mutate(dx_sub=levels(newdata$dx_sub)[1]),type="lpmatrix",allow.new.levels=TRUE)
  tmp[,grepl("dx_sub",colnames(tmp))] <- 0
  out[(newdata$dx_sub %in% out_lev)]  <- plogis(tmp%*%coef(mod))


  })
}
  return(out)
}
