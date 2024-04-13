#############################################################
############# Custom Functions of Maki ######################
#############################################################

######## eulerLotka() #####
#' eulerLotka
#'
#' @description
#' based on a bisection method (\code{uniroot} -> searches for 0)
#'
#' The Euler-Lotka equation:
#' 1 = SUM[ l(x) * m(x) * exp(-r * x) ]
#'
#' l = survival rate until age x;
#'
#' m = amount of offspring at age x;
#'
#' x = age of reproduction;
#'
#' r = population growth rate;
#'
#' using the Euler-Lotka equation:
#' when calling this function, use cumulative survival rate for clutch 1-3 (l),
#' amount of offspring at clutch 1-3 (m) and the age of the reproducing individual at clutch 1-3 (x)
#' if all 3 values of those are 0, the function will set one to 0.0001 (otherwise it won't work)
#' if one value is NA, this function will return NA
#'
#' seed is needed for reproducible results (otherwise random numbers will never be the same!)
#' you can take every number you like to
#'
#' @param data Data.frame containing your data.
#' @param col_l1 Column name of the survival rate until age 1.
#' @param col_l2 Column name of the survival rate until age 2.
#' @param col_l3 Column name of the survival rate until age 3.
#' @param col_m1 Column name of the amount of offspring at age 1.
#' @param col_m2 Column name of the amount of offspring at age 2.
#' @param col_m3 Column name of the amount of offspring at age 3.
#' @param col_x1 Column name of the age of reproduction at age 1.
#' @param col_x2 Column name of the age of reproduction at age 2.
#' @param col_x3 Column name of the age of reproduction at age 3.
#'
#' @return Returns a vector that has the same order like your data.frame. You can just add this vector to your data.frame as a new column.
#'
#' @export
#' @importFrom stats uniroot
#'
eulerLotka <- function(data, col_l1, col_l2, col_l3, col_m1, col_m2, col_m3, col_x1, col_x2, col_x3){
  #iterate over dataframe and create a vector that will be returned
  pop_gr <- vector(mode="numeric", length=0)

  for (i in 1:length(data[,col_l1])){
    l1 <- data[,col_l1][i]
    l2 <- data[,col_l2][i]
    l3 <- data[,col_l3][i]
    m1 <- data[,col_m1][i]
    m2 <- data[,col_m2][i]
    m3 <- data[,col_m3][i]
    x1 <- data[,col_x1][i]
    x2 <- data[,col_x2][i]
    x3 <- data[,col_x3][i]

    #if one of those values is NA, return NA and don't calculate anything
    if(is.na(l1) || is.na(l2) || is.na(l3) || is.na(m1) || is.na(m2) || is.na(m3) || is.na(x1) || is.na(x2) || is.na(x3)){
      r <- NA
    }
    else{
      # if all 3 values of a parameter are 0, the function will set one to 0.0001 (otherwise it won't work)
      if(l1 == 0 && l2 == 0 && l3 == 0){
        l1 <- 0.0001
      }
      if(m1 == 0 && m2 == 0 && m3 == 0){
        m1 <- 0.0001
      }
      if(x1 == 0 && x2 == 0 && x3 == 0){
        x1 <- 0.0001
      }

      #this is the Lotka-Euler-Eq -z and z becomes 1 as uniroot finds a 0-value
      y <- function(r, l1, l2, l3, m1, m2, m3, x1, x2, x3, z){((l1*m1*exp(-r*x1)) + (l2*m2*exp(-r*x2)) + (l3*m3*exp(-r*x3))) - z}

      # uniroot finds a 0 value, so offset function, thats why -z in the upper formula
      r <- stats::uniroot(y, l1=l1, l2=l2, l3=l3, m1=m1, m2=m2, m3=m3, x1=x1, x2=x2, x3=x3, z = 1, interval = c(-10, 10), extendInt = "yes")$root #writing only the result of r into variable
    }

    pop_gr[i] <- r
  }
  print("Calculated population growth rate for all lines.")
  # return r into table
  return(pop_gr)

} #end of eulerLotka function



######## DIndex.OneStep() ################
#' @description
#' Calculates the defence index in one step using dTransform, DIndex and dContribution.
#' TODO: not yet implemented. Don't actually know whether this makes sense.
#'
#' @param data the data.frame to be used.
#' @param treatments_column_name The name of the column that contains the treatments as factors.
#' @param control_name The name of the factor for your control group. This is important for standardization.
#' @param column_name_desired_value The name of the value (column name) you like to transform.
#' @param usemc 1 = use mean control to standardize the values. In this case, values will be 
#' standardized for the experiment. Therefore, values are completely independent of the experimental
#' design and other factors. This increases the comparability with other studies. However,
#' when using this for defence index calculation, the weighing (dmax) of the traits will be 1 in each case.
#' Therefore, no real weighing takes place. Defaults to 0 (standardize by the span of 
#' treatment and control (max-min)).
#'
#' @return Returns a list containing the contribution of each trait. For each
#'         treatment a separate data frame is created and stored within the list.
#'         Furthermore, an automated, qualitative interpretation is provided.
#'
#' @export
#' @importFrom 
#'
DIndex.OneStep <- function(data, treatments_column_name, control_name, 
                           column_name_desired_value, usemc = 0, 
                           traits.adapt, traits.mal = ""){
  
  
}


######## simD() ################
#' @description
#' This function simulates random data according to settings and applies the Defence Index framework including the contribution 
#' of each trait to DI. The aim is to get an estimate how precise the algorithm works with different data properties. 
#' A data set containing different settings can be set up and iterated to run this simulation in iterations.
#'
#' @param baseData not necessary any more.
#' @param nTraits The number of traits that should be created.
#' @param nReplicates Number of replicates for each treatment.
#' @param nAdaptive Number of traits that are determined to be adaptive. Those traits will be simulated to exhibit an 
#' increase in the treatment. In the end, those traits should be recognized as adaptive from the algorithm to have a correct count.
#' The first x traits will be set as adaptive.
#' @param nMaladaptive Number of traits that are determined to be maladaptive. Those traits will be simulated to exhibit a 
#' decrease in the treatment. In the end, those traits should be recognized as adaptive from the algorithm to have a correct count.
#' The last x traits will be set as maladaptive.
#' @param nNonNormal Number of non-gaussian distributed traits. The first x traits will be non-gaussian distributed.
#' @param diffBetweenTreats This parameters defines how strong the traits can differ in the treatment from control when adaptive or
#' maladaptive. However, in the algorithm we included a random factor, which actually make this term more variable to a certain extend.
#' This is thought to create more realistic data properties.
#' @param th.change In the estimation of the contribution to DI we apply two thresholds. This parameter is one of them and determines
#' at which change in contribution the respective trait is estimated to be (mal)adaptive. These thresholds were used to optimize the 
#' algorithm. In both cases a threshold of 0.05 achieved the best results.
#' @param th.sd In the estimation of the contribution to DI we apply two thresholds. This parameter is one of them and determines
#' at which sd of the treatments' contribution an interaction is assumed. These thresholds were used to optimize the 
#' algorithm. In both cases a threshold of 0.05 achieved the best results.
#' 
#' @return Returns a list containing the overall distance change of each point to the control mean as a measure for the performance
#' of the transformation method. Additionally, it contains how many estimations were correct and incorrect.
#'
#' @export
#' @importFrom 
#'

simD <- function(baseData, nTraits, nReplicates, nAdaptive, nMaladaptive, nNonNormal = 0, diffBetweenTreats, th.change, th.sd){
  # object to save all values of interest
  obj <- list()
  
  ##### simulate data set ######
  # we will use the same base for each trait and randomly change its magnitude. The SD is set to 10% of the mean + a random number (-> SD will be randomly between 5 and 20 % of the mean)
  
  # create simulated data for each treatment
  data.sim.control <- data.frame(treatment = c(1:nReplicates))
  data.sim.control$treatment <- "control"
  # create control values for each trait
  skew=runif(nTraits, min = 0.5, max = 1) 
  kurt=runif(nTraits, min = 0.5, max = 2)
  
  for(i in 1:nTraits){
    mean <- 1*runif(1, min = 0.01, max = 10000)
    if(nNonNormal >= i){ # if there should be non-normal data, the first traits will be non-normal
      # generate non-normal data using the third-order polynomial method described by Vale & Maurelli (1983) and Fleishman (1978)
      # https://www.rdocumentation.org/packages/SimDesign/versions/2.6/topics/rValeMaurelli
      # skew between 0.5 and 1 
      # kurt between 0.5 and 2 (otherwise it gets too close to normal distribution)
      data.sim.control[,paste("trait", i, sep="_")] <- rValeMaurelli(nReplicates, mean=mean, sigma=mean*(0.1+runif(1, min = -0.05, max = 0.1)), skew=skew[i], kurt=kurt[i])
      
    }
    else{ # create normal distributed data
      
      data.sim.control[,paste("trait", i, sep="_")] <- rnorm(nReplicates, mean = mean, sd = mean*(0.1+runif(1, min = -0.05, max = 0.1)))
      
    }
  }
  # create treatment values for each trait
  # these values should rely on the same values as the control but maybe with a slight difference
  data.sim.treat <- data.frame(treatment = c(1:nReplicates))
  for(i in 1:nTraits){
    # get the factor for the differences between traits (input: 0.2 for 20% -> factor = 0.2+1 = 1.2)
    if(nAdaptive >= i){
      diff <- diffBetweenTreats + 1
    }
    else if(nMaladaptive >= (nTraits - i + 1)){
      diff <- 1 - diffBetweenTreats
    }
    else{
      diff <- 1
    }
    data.sim.treat$treatment <- "treatment"
    
    if(nNonNormal >= i){ # if there should be non-normal data, the first traits will be non-normal
      
      # skew between 0.5 and 1 
      # kurt between 0.5 and 2 (otherwise it gets too close to normal distribution)
      data.sim.treat[,paste("trait", i, sep="_")] <- rValeMaurelli(nReplicates, mean=mean(data.sim.control[, i+1])*diff, sigma=(sd(data.sim.control[, paste("trait", i, sep="_")])*diff)*(0.1+runif(1, min = -0.1, max = 0.1)), skew=skew[i], kurt=kurt[i])
      
    }
    else{ # create normal distributed data
      
      data.sim.treat[,paste("trait", i, sep="_")] <- rnorm(nReplicates, mean = mean(data.sim.control[, i+1])*diff, sd = (sd(data.sim.control[, paste("trait", i, sep="_")])*diff)*(0.1+runif(1, min = -0.1, max = 0.1)))
      
    }
  }
  
  # merge treatments into one data set
  data.sim <- rbind(data.sim.control, data.sim.treat)
  
  ##### transform values ######
  #transform values to d and z (as benchmark)
  for(i in 1:nTraits){
    # d transform values 
    data.sim[,paste("trait", i, "d", sep="_")] <- InARes::Trex(data.sim, "treatment", "control", paste("trait", i, sep="_"))
    # z transform values as reference
    data.sim[,paste("trait", i, "z", sep="_")] <- z_transform(data.sim, "treatment", paste("trait", i, sep="_"))
  }
  
  ##### benchmark transformation ######
  # To benchmark the transformation, and how it affects the data in comparison to the z transformation
  # we calculate the standardized Euclidean distance from each point to the mean of the control
  # prepare data frame and such
  dist_change <- data.frame(trait = c(1:nTraits),
                            dtransformed_mean = 0,
                            dtransformed_sd = 0,
                            ztransformed_mean = 0,
                            ztransformed_sd = 0
  )
  params <- c()
  params_d <- c()
  params_z <- c()
  for(i in 1:nTraits){
    dist_change$trait[i] <- paste("trait", i, sep="_")
    
    params = c(params, paste("trait", i, sep="_"))
    params_d = c(params_d, paste("trait", i, "d", sep="_"))
    params_z = c(params_z, paste("trait", i, "z", sep="_"))
  }
  
  # calculate the change of the standardized euclidean distance of each value to the mean of control from the values before and after transformation
  # container to gather the dist_change of each trait
  overallDistChange_d <- c()
  overallDistChange_z <- c()
  # iterate overt the traits and calculate the distance of each value to the mean of control for each:
  # untransformed, d transformend and z transformed values
  # it is: (value(x) - mean(control))^2 / sd(values(x))^2
  
  for(i in 1:length(params)){
    dist <- ((data.sim[,params[i]] - mean(data.sim[which(data.sim$treatment == "control"),][,params[i]]))^2 ) / sd(data.sim[,params[i]])^2
    dist_d <- ((data.sim[,params_d[i]] - mean(data.sim[which(data.sim$treatment == "control"),][,params_d[i]]))^2 ) / sd(data.sim[,params_d[i]])^2
    dist_z <- ((data.sim[,params_z[i]] - mean(data.sim[which(data.sim$treatment == "control"),][,params_z[i]]))^2 ) / sd(data.sim[,params_z[i]])^2
    
    # then calculate the mean and sd of the change from the untransformed to the transformed data
    # to have an idea of whether it strongly affects the relation of the points, we compare these values with a common z-transformation
    dist_change$dtransformed_mean[i] <- mean(dist - dist_d)
    dist_change$dtransformed_sd[i] <- sd(dist - dist_d)
    dist_change$ztransformed_mean[i] <- mean(dist - dist_z)
    dist_change$ztransformed_sd[i] <- sd(dist - dist_z)
    
  }
  
  # take the mean of all treatments
  # the unit is standard deviation (a value of 1 means the relative position to the control mean changes by one standard deviation in mean)
  overallDistChange_d <- mean(dist_change$dtransformed_mean)
  overallDistChange_z <- mean(dist_change$ztransformed_mean)
  #cat("Mean change when dtransformed: ", mean(dist_change$dtransformed_mean), " +/- ", sd(dist_change$dtransformed_mean))
  #cat("Mean change when ztransformed: ", mean(dist_change$ztransformed_mean), " +/- ", sd(dist_change$ztransformed_mean))
  
  obj$overallDistChange_d <- overallDistChange_d
  obj$overallDistChange_z <- overallDistChange_z
  
  
  
  ###### calculate the defence index D ####
  # to have unbiased estimates, the adaptive and maladaptive traits are included in traits.adapt and as.defined = FALSE (default). Therefore, it is not already provided, which traits are adaptive or maladaptive. 
  traits.adapt = params_d
  #traits.mal = ""
  data.sim$D <- InARes::InARes(data.sim, traits.adapt = traits.adapt)$D
  
  
  ###### evaluate the traits' contribution to D ########
  # ensure that treatment is set as factor (important for next step)
  data.sim$treatment <- as.factor(data.sim$treatment)
  # evaluate contributions
  contrData <- InARes::contARes(data.sim, params_d, "D", "treatment", "control", th.change = th.change, th.sd = th.sd)
  
  
  ####### check correctness ##########
  # values for correct or incorrect evaluation
  # only necessary if one would put all traits together without knowing which trait is maladaptive
  correct <- 0
  incorrect <- 0
  
  # iterate over traits and check in contrData how often the evaluation was correct/incorrect
  # if adaptive traits are desired, the first traits will be adaptive
  # if maladaptive traits are desired, they will be the last traits
  # therefore we can check correctness via the index i
  for(i in 1:nTraits){
    
    switch(contrData[[ paste("treatment", "summary", params_d[i], sep = ".")]],
           
           "An increase in this trait seem to be adaptive. However, there might be interactions with other traits." = {
             if(nAdaptive >= i){
               correct <- correct +1
             }
             else{
               incorrect <- incorrect +1
             }
           },
           "An increase in this trait seem to be adaptive." = {
             if(nAdaptive >= i){
               correct <- correct +1
             }
             else{
               incorrect <- incorrect +1
             }
           },
           "An increase in this trait seem to be maladaptative. However, there might be interactions with other traits." = {
             if(nMaladaptive >= (nTraits - i + 1)){
               correct <- correct +1
             }
             else{
               incorrect <- incorrect +1
             }
           },
           "An increase in this trait seem to be maladaptative." = {
             if(nMaladaptive >= (nTraits - i + 1)){
               correct <- correct +1
             }
             else{
               incorrect <- incorrect +1
             }
           },
           "This trait is indefinite in its change in contribution. There might be an interaction with other traits." = {
             if(nMaladaptive < (nTraits - i + 1) && nAdaptive < i){
               correct <- correct +1
             }
             else{
               incorrect <- incorrect +1
             }
           },
           "There seem no (mal-)adaptive function in this trait." = {
             if(nMaladaptive < (nTraits - i + 1) && nAdaptive < i){
               correct <- correct +1
             }
             else{
               incorrect <- incorrect +1
             }
           }
    )
  } # end check correctness
  obj$correct <- correct
  obj$incorrect <- incorrect
  
  
  
  return(obj)
}



######## z_transform() #################
#' z_transform
#' @description
#' A function for calculating z-transformation (standardization) of certain input (a column of a df):
#' z is a transformed value (standardized) that is calculated as:
#' z = X(individual)-mean(x(treatment)) / sd(x(treatment))
#'
#' @param data A data.frame to be used.
#' @param treatments_column_name The name of the column that contains the treatments as factors.
#' @param value A value (column name) you like to transform.
#'
#' @return Returns a data.frame containing the transformed data.
#'
#' @export
#' @importFrom stats sd
#'
z_transform <- function(data, treatments_column_name, value){
  
  #z-transformation (standardisation):
  # z = X-mean(x) / stabw(x)
  
  treatments <- unique(data[,treatments_column_name])
  
  # z-transform data
  z_data <- data.frame()
  for (i in 1:length(treatments)){
    z_temp <- subset(data, data[,treatments_column_name] == treatments[i])
    
    for (j in 1:nrow(z_temp)){
      z_temp$zval[j] <- (z_temp[,value][j] - mean(z_temp[,value])) / sd(z_temp[,value])
    }
    z_data <- rbind(z_data, z_temp)
  }
  
  #print("z-transformend data (standardized) added as additional column to your data set.")
  return(z_data$zval)
}


######## as_numeric_chain() ###########
#' as_numeric_chain
#' @description
#' A function for set as.numeric for several colums at once.
#'
#' @param data A data.frame containing the data with the desired columns.
#' @param col_name_list A list ('c()') of column names to be selected for conversion to as.character.
#'
#' @return Returns the data.frame with converted columns.
#' @export
#'
as_numeric_chain <- function(data, col_name_list){
  # iterate over all columns of the list
  for (i in 1:length(col_name_list)){
    data[,col_name_list[i]] <- as.numeric(as.character(data[,col_name_list[i]]))
  }
  cat("Changed columns to numeric.")
  return(data)
}



######## as_factor_chain() ###########
#' as_factor_chain
#' @description
#' A function for set as.factor for several colums at once.
#'
#' @param data A data.frame containing the data with the desired columns.
#' @param col_name_list A list ('c()') of column names to be selected for conversion to as.character.
#'
#' @return Returns the data.frame with converted columns.
#' @export
#'
as_factor_chain <- function(data, col_name_list){
  # iterate over all columns of the list
  for (i in 1:length(col_name_list)){
    data[,col_name_list[i]] <- as.factor(data[,col_name_list[i]])
  }
  cat("Changed columns to factor")
  return(data)
}



######## as_character_chain() ###########
#' as_character_chain
#' @description
#' A function for set as.character for several colums at once.
#'
#' @param data A data.frame containing the data with the desired columns.
#' @param col_name_list A list ('c()') of column names to be selected for conversion to as.character.
#'
#' @return Returns the data.frame with converted columns.
#' @export
#'
as_character_chain <- function(data, col_name_list){
  # iterate over all columns of the list
  for (i in 1:length(col_name_list)){
    data[,col_name_list[i]] <- as.character(data[,col_name_list[i]])
  }
  cat("Changed columns to character")
  return(data)
}



######## prep_data_dodged_bargraph() ###########
#' prep_data_dodged_bargraph
#' @description
#' This function will melt your data in a way that allows you to set up a bar graph
#' with dodged bars where the groups are single treatmens and the bars contain data of your selected values
#' NAs will be ignored.
#'
#' @param data A data.frame containing the following parameters
#' @param treatment_column_name The name of the column where your treatments are saved as factor
#' @param value_list A list ('c()') of the values (columns) that you like to include in your graph
#' @param treat_names_list A list ('c()') of treatment names that you would like to rename for your graph
#' @param treat_rename_list A list ('c()') of new treatment names for the treat_names_list. It should have the same order as treat_names_list!
#'
#' @return Returns a data.frame that contains the melted values of your input
#'
#' @export
#' @importFrom reshape2 melt
#' @importFrom stats aggregate
#' @importFrom stats sd
#'
prep_data_dodged_bargraph <- function(data, treatment_column_name, value_list, treat_names_list = c(), treat_rename_list = c()){
  # to plot everything into one bargraph we need reshape 2
  # the data has to be summarized (mean and sd of a each value for each treatment)
  # additionally, the treatments can be renamed
  # NAs will be ignored/removed by aggregate
  #requireNamespace(reshape2, quietly = TRUE)

  data_sum <- data.frame(treatment = levels(data[, treatment_column_name]))

    # aggregate data (mean and sd), and set into new columns
  for(i in 1:length(value_list)){
    data_sum[, paste(value_list[i], "mean", sep="_")] <- stats::aggregate(data[c(value_list)], by = list(data[, treatment_column_name]), FUN = mean, na.rm = TRUE)[, value_list[i]]
    data_sum[, paste(value_list[i], "sd", sep="_")] <- stats::aggregate(data[c(value_list)], by = list(data[, treatment_column_name]), FUN = stats::sd, na.rm = TRUE)[, value_list[i]]
  }

  # define which columns should be melted
  melt_list_mean <- vector(mode="character", length=0)
  melt_list_mean[1] <- "treatment"
  melt_list_sd <- vector(mode="character", length=0)
  melt_list_sd[1] <- "treatment"
  for(i in 1:length(value_list)){
    melt_list_mean[i+1] <- paste(value_list[i], "mean", sep="_")
    melt_list_sd[i+1] <- paste(value_list[i], "sd", sep="_")
  }

  # then data have to be melted correctly
  dfm <- data.frame()
  dfm_sd <- data.frame()
  dfm <- reshape2::melt(data_sum[, melt_list_mean], id.vars = 1)
  dfm_sd <- reshape2::melt(data_sum[, melt_list_sd], id.vars = 1)
  dfm$sd <- dfm_sd$value

  # rename treatments
  for(i in 1:length(treat_names_list)){
    levels(dfm$treatment)[levels(dfm$treatment)==treat_names_list[i]] <- treat_rename_list[i]
  }

  return(dfm)

  # # Now we can plot the data into one graph (however, I didn't manage to reorder the groups, yet)
  # ggplot(data=dfm, aes(x=treatment, y = value, fill = variable))+
  #   geom_bar(stat = "identity", position = "dodge")+
  #   geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width = .5, position = position_dodge(.9))+
  #   theme(axis.text.x = element_text(size=15, angle=45, vjust=0.4), axis.text.y = element_text(size = 15), axis.title = element_text(size = 20), plot.title = element_text(hjust = 0.5), legend.position="right", panel.background = element_rect(fill = "white"), axis.line = element_line(size = 0.5))+
  #   #coord_cartesian(ylim = c(0, 9000))+ # change limits of coordinate system
  #   xlab("Treatment") + # label x-axis
  #   scale_x_discrete(c("control", "kairomone extract", "predator water", "direct induction"))+
  #   ylab("size of measurements in mm")+ # label y-axis
  #   scale_fill_grey(# change legend
  #     #values=c("#999999", "#E69F00", "#56B4E9", ""), # colors to fill bars
  #     name="Measurement", # header
  #     breaks=c("body_length_primpar_?m_mean", "body_width_primpar_?m_mean", "spine_length_primpar_?m_mean", "rel_spine_mean"), # original label
  #     labels=c("body length", "body width", "spine length", "relative spine length") # custom label
  #   )+ # change colors into B/W
  #   ggtitle("Defences of D. magna against Triops")
  # #ggsave("Defences of D. magna against Triops.jpg", height = 25, width = 45, units = "cm", dpi = 300, device = "jpg", path = "C:/Users/marvi_000/Desktop/Masterarbeit/Data R analysis/Plots")
  #

}

######## subset_multiple_factors(data, col_name, selected_factors_list) ###########
#' subset_multiple_factors
#' @description
#' Function to subset a dataset with several selected factors (of one column!). Not used levels will be dropped!
#'
#' @param data A data.frame containing the data of the following parameters.
#' @param col_name The column name of the column containing the factors you want to subset for.
#' @param selected_factors_list A list ('c()') of factors from the selected column you want to subset for.
#'
#' @return Returns a data.frame containing the subsetted data
#' @export
#'
subset_multiple_factors <- function(data, col_name, selected_factors_list){

  sub_data <- data.frame()

  for(i in 1:length(selected_factors_list)){
    sub_data <- rbind(sub_data, droplevels(subset(data, data[, col_name] == selected_factors_list[i])))
  }

  return(sub_data)
}

######## MakiCV() ###########
#' created by Marvin Kiene, use this function at own risk!
#' @description
#' For cross validation we split the dataset into k parts. 
#' Then iterating over these parts, each time using k-1 parts for 
#' model creation and one part for testing Brier score. 
#' Brier score will be our criterion in case of binomial models.
#' RMSE is the criterion in case of other distributions.
#' The data are splitted randomly. Therefore, it is recommended not to use this function
#' with small datasets with a lot of levels. However, setting k = nrow(data) and rept = 1
#' will perform a leave-one-out cross validation (LOOCV), which might be the better choice
#' for small datasets. But consider, that LOOCV is computaionally intense.
#' 
#' @param data a data set to operate with
#' @param mod_func function to be used for modelling. Currently possible model functions are: gamm4::gamm4, lme4::glmer, lme4::lmer, lm, glm
#' @param fam family to be used for the model (has to be provided as used in model function)
#' @param k determine k-fold cross validation (in how many parts the data set should be splitted)
#' @param rept repetitions for the k-fold cross validation
#' @param params parameterisation of fixed terms (including smooths) and random effects in case of lme4 models (except gamm4)
#' @param paramRdm random terms, for nlme models (have to be provided as for model function necessary)
#' @param it number of max iterations (for maxfun)
#' @param response response variable for prediction
#' @param setseed whether to set a seed or not. If T, outcome will be reproducible. Else sampling is random each time. Defaults to F.
#' 
#' @return returns an object with summary and score values of each model run and cross check
#' @export
#'
MakiCV <- function(data, mod_func, fam="gaussian", it=1e6, k=5, rept=3, params, paramRdm = "", response, setseed=F){
  
  if(setseed == T){set.seed(42)}
  
  paramFixed <- params
  # for cross validation we split the dataset into 5 parts. Then iterating over these parts, each time using 4 parts for model creation and one part for testing Brier score. Brier score will be our criterion in case of binomial models.
  # The goal is a k-fold cross validation with rept repeats.
  
  # check which model function should be used and then use the respective cv-algorythm
  # currently possible: gamm4, glmer

  # check whether selected function can be performed
  if(mod_func == "gamm4" || 
     mod_func == "glmer" || mod_func == "lmer" || 
     mod_func == "lm" || mod_func == "glm"){
  # data frame to store scores after model creation and prediction
    bs <- data.frame(x = NA)
    # get the full data set line numbers in a vector
    datafull <- sample(nrow(data), nrow(data))
    
    # vectors for summary values (means and SDs)
    statmeans <- c()
    statsds <- c()
    
    #### Data splitting #####
    for(h in 1:rept){
      
      # randomly split dataset into k parts
      # list to save the splitted dataset
      obj <- list()
      samples <- matrix(nrow= floor(nrow(data)/k), ncol= k)
      obj$samples <- samples
      
      # vector of line numbers that are used in other parts of the splitted data set
      notuseln <- c()
      
      # iterate over k and save line numbers of splitted data set into obj
      for(i in 1:k){
        
        #check which lines are already used and which are still to be used
        for(j in 1:i){
          if(i == 1){} # do nothing in the first iteration (no data selected)
          else{ # get line numbers that are used so far
            notuseln <- c(notuseln, obj$samples[,j])
          }
        }
        # get line numbers of one part of the splitted data set
        # check which lines are already used
        obj$samples[,i] <- sample(setdiff(datafull, 
                                          notuseln), 
                                  floor(nrow(data)/k)
        )
      }# end for i
      
      # add a column for current repetition to bs
      bsrep <- c() 
      
      #### iteration over k models
      # iterate over k, create a model using k-1 subsamples and predict with the left subsample to calculate Brier score
      for(i in 1:k){
        
        # provide info on progress
        cat("Progress of cross validation:", (((i-1)+k*(h-1))/(k*rept))*100, "% \n")
        
        # get data of k-1 of the subsamples for model creation/training by leaving out the kth subsample
        moddata <- data[-c(obj$samples[,i]),]
          
        ##### gamm4 ########
        if(mod_func == "gamm4"){
          # run model
          mod.cv <- gamm4::gamm4(paramFixed,
                                 random = paramRdm,
                                 data = moddata, family = fam, REML=TRUE,
                                 control=lme4::glmerControl(optCtrl=list(maxfun = it)))
          
          # make prediction with this model, but with the leftover data set
          testdata <- data[obj$samples[,i],]
          fit.cv = predict(mod.cv$gam, newdata = testdata, type = 'response', se.fit = FALSE)
          
          if(fam == "binomial"){ # if binomial, calculate Brier score
            # calculate Brier score Brier score for this prediction 
            # Bier score = mean((pred.prob-event_success_bin)^2)
            # following https://stackoverflow.com/questions/25149023/how-to-find-the-brier-score-of-a-logistic-regression-model-in-r
            bsrep <- c(bsrep, mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2))
          }
          else{ # else calculate RMSE
            # RSME = sqrt(mean((fitted-observed)^2))
            bsrep <- c(bsrep, sqrt(mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2)))
          }  
          
        } # end if(mod_func == "gamm4")
        ##### glmer or lmer ######
        if(mod_func == "glmer" || mod_func == "lmer"){
          # run model
          if(fam == "gaussian" || mod_func == "lmer"){
            mod.cv <- lme4::lmer(paramFixed,
                                 data = moddata,
                                 control=lme4::lmerControl(optCtrl=list(maxfun = it)))
          }
          else{
            mod.cv <- lme4::glmer(paramFixed,
                                  data = moddata, family = fam,
                                  control=lme4::glmerControl(optCtrl=list(maxfun = it)))
          }
          # make prediction with this model, but with the leftover data set
          testdata <- data[obj$samples[,i],]
          #cat("predict")
          fit.cv = predict(mod.cv, newdata = testdata, allow.new.levels = TRUE, type = 'response', se.fit = FALSE)
          
          #### calculate and gather scores for k-folds #####
          if(fam == "binomial"){ # if binomial, calculate Brier score
            # calculate Brier score Brier score for this prediction 
            # Bier score = mean((pred.prob-event_success_bin)^2)
            # following https://stackoverflow.com/questions/25149023/how-to-find-the-brier-score-of-a-logistic-regression-model-in-r
            bsrep <- c(bsrep, mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2))
          }
          else{ # else calculate RMSE
            # RSME = sqrt(mean((fitted-observed)^2))
            bsrep <- c(bsrep, sqrt(mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2)))
          }  
        } # end if(mod_func == "glmer" || mod_func == "lmer")  
        ##### lm or glm ######
        if(mod_func == "glm" || mod_func == "lm"){
          # run model
          if(fam == "gaussian" || mod_func == "lm"){
            mod.cv <- lm(paramFixed,
                         data = moddata,
                         control = list(maxit = it, epsilon=1))
          }
          else{
            mod.cv <- glm(paramFixed,
                          data = moddata, family = fam,
                          control = list(maxit = it, epsilon=1))
          }
          # make prediction with this model, but with the leftover data set
          testdata <- data[obj$samples[,i],]
          cat("predict")
          fit.cv = predict(mod.cv, newdata = testdata, allow.new.levels = TRUE, type = 'response', se.fit = FALSE)
          
          #### calculate and gather scores for k-folds #####
          if(fam == "binomial"){ # if binomial, calculate Brier score
            # calculate Brier score Brier score for this prediction 
            # Bier score = mean((pred.prob-event_success_bin)^2)
            # following https://stackoverflow.com/questions/25149023/how-to-find-the-brier-score-of-a-logistic-regression-model-in-r
            bsrep <- c(bsrep, mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2))
          }
          else{ # else calculate RMSE
            # RSME = sqrt(mean((fitted-observed)^2))
            bsrep <- c(bsrep, sqrt(mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2)))
          }  
        } # end if(mod_func == "glmer" || mod_func == "lmer")  
      }# end i
      
      ### gather values of rept ####  
      # bind values together in data frame and rename added column
      bs <- cbind(bs, bsrep)
      colnames(bs)[which(colnames(bs) == "bsrep")] <- c(paste("rep", h, sep=""))
        
    }# end h
    
    bs$x <- NULL
    
    ##### create an object summarizing the performed CV and values #####
    summarycv <- list()
    summarycv$cvInfo <- paste(k,"-fold cross validation with ", rept, " repeats", sep="")
    summarycv$modInfo <- paste("Using model ", mod_func, " with family = ", fam, sep="")
    if(fam == "binomial"){
      summarycv$scoreInfo <- paste("For ", fam, " models, Brier score is used as criterion.", sep="")
    }
    else{
      summarycv$scoreInfo <- paste("For ", fam, " models, RMSE is used as criterion.", sep="")
    }
    statsum <- data.frame(repetition = colnames(bs), mean = NA, sd = NA)
    for(i in 1:rept){
      statsum$mean[i] <- mean(bs[, i])
      statsum$sd[i] <- sd(bs[, i])
    }
    summarycv$statSummary <- statsum
    summarycv$cvResults <- bs
    fullmean <- mean(t(summarycv$cvResults))
    fullsd <- sd(t(summarycv$cvResults))
    if(fam == "binomial"){
      summarycv$summaryInfo <- paste("Over all runs, the model had a mean acuracy of ", 
                                     round((1-(fullmean*2))*100, digits = 2),
                                     " % (Brier score of ",
                                     round(fullmean, digits = 4),
                                     ") and a standard deviation of ",
                                     round(fullsd*2*100, digits = 2),
                                     " % (SD of Brier score of ",
                                     round(fullsd, digits = 4),
                                     ")",
                                     sep="")
    }
    else{
      summarycv$summaryInfo <- paste("Over all runs, the model had a mean RMSE of ",
                                     round(fullmean, digits = 4),
                                     " and a standard deviation of RMSE of ",
                                     round(fullsd, digits = 4),
                                     ")",
                                     sep="")
    }
      
    cat(crayon::blue("100 % done!"))
    return(summarycv)
    
  } # if(mod_func == ...)
  else{
    
    cat(crayon::red("ERROR: The selected model function is not available or misspelled."))
    
  }
}

######## MakiCV.nlme() ###########
#' created by Marvin Kiene, use this function at own risk!
#' @description
#' For cross validation we split the dataset into k parts. 
#' Then iterating over these parts, each time using k-1 parts for 
#' model creation and one part for testing Brier score. 
#' Brier score will be our criterion in case of binomial models.
#' RMSE is the criterion in case of other distributions.
#' The data are splitted randomly. Therefore, it is recommended not to use this function
#' with small datasets with a lot of levels. However, setting k = nrow(data) and rept = 1
#' will perform a leave-one-out cross validation (LOOCV), which might be the better choice
#' for small datasets. But consider, that LOOCV is computaionally intense.
#' 
#' @param data a data set to operate with
#' @param mod_func function to be used for modelling. Currently possible model functions are: mgcv::gamm
#' @param fam family to be used for the model (has to be provided as used in model function)
#' @param k determine k-fold cross validation (in how many parts the data set should be splitted)
#' @param rept repetitions for the k-fold cross validation
#' @param params.fixed parameterisation of fixed terms (including smooths) and random effects in case of lme4 models (except gamm4)
#' @param params.rdm random terms, for nlme models (have to be provided as for model function necessary)
#' @param it number of max iterations (for maxfun)
#' @param response response variable for prediction
#' @param setseed whether to set a seed or not. If T, outcome will be reproducible. Else sampling is random each time. Defaults to F.
#' 
#' @return returns an object with summary and score values of each model run and cross check
#' @export
#'
MakiCV.nlme <- function(data, 
                        mod_func, 
                        fam="gaussian", 
                        it=1e6, 
                        k=5, 
                        rept=3, 
                        params.fixed, 
                        params.rdm = "", 
                        response, 
                        correlation = FALSE, 
                        weights = FALSE, 
                        setseed=F){
  
  if(setseed == T){set.seed(42)}
  
  # for cross validation we split the dataset into 5 parts. Then iterating over these parts, each time using 4 parts for model creation and one part for testing Brier score. Brier score will be our criterion in case of binomial models.
  # The goal is a k-fold cross validation with rept repeats.
  
  # check which model function should be used and then use the respective cv-algorythm
  # currently possible: gamm4, glmer
  
  # check whether selected function can be performed
  if(#mod_func == "lme" || 
     #mod_func == "gls" || 
     #mod_func == "gam" || 
     mod_func == "gamm"
     ){
    # data frame to store scores after model creation and prediction
    bs <- data.frame(x = NA)
    # get the full data set line numbers in a vector
    datafull <- sample(nrow(data), nrow(data))
    
    # vectors for summary values (means and SDs)
    statmeans <- c()
    statsds <- c()
    
    #### Data splitting #####
    for(h in 1:rept){
      
      # randomly split dataset into k parts
      # list to save the splitted dataset
      obj <- list()
      samples <- matrix(nrow= floor(nrow(data)/k), ncol= k)
      obj$samples <- samples
      
      # vector of line numbers that are used in other parts of the splitted data set
      notuseln <- c()
      
      # iterate over k and save line numbers of splitted data set into obj
      for(i in 1:k){
        
        #check which lines are already used and which are still to be used
        for(j in 1:i){
          if(i == 1){} # do nothing in the first iteration (no data selected)
          else{ # get line numbers that are used so far
            notuseln <- c(notuseln, obj$samples[,j])
          }
        }
        # get line numbers of one part of the splitted data set
        # check which lines are already used
        obj$samples[,i] <- sample(setdiff(datafull, 
                                          notuseln), 
                                  floor(nrow(data)/k)
        )
      }# end for i
      
      # add a column for current repetition to bs
      bsrep <- c() 
      
      #### iteration over k models
      # iterate over k, create a model using k-1 subsamples and predict with the left subsample to calculate Brier score
      for(i in 1:k){
        
        # provide info on progress
        cat("Progress of cross validation:", (((i-1)+k*(h-1))/(k*rept))*100, "% \n")
        
        # get data of k-1 of the subsamples for model creation/training by leaving out the kth subsample
        moddata <- data[-c(obj$samples[,i]),]
        
        ##### gamm ########
        if(mod_func == "gamm"){
          # run/create model
          # I have to check whether correlations or weights are included. Otherwise the models don't work with empty parameters
          if(correlation == FALSE && weights == FALSE){
            mod.cv <- mgcv::gamm(params.fixed,
                                   random = params.rdm,
                                   data = moddata, 
                                   family = fam, 
                                   method = "REML",
                                  control=mgcv::gam.control(maxit = it))
          }
          else if(correlation != FALSE && weights == FALSE){
            mod.cv <- mgcv::gamm(params.fixed,
                                 random = params.rdm,
                                 data = moddata, 
                                 family = fam, 
                                 method = "REML",
                                 correlation = correlation,
                                 control=mgcv::gam.control(maxit = it))
          }
          else if(correlation == FALSE && weights != FALSE){
            mod.cv <- mgcv::gamm(params.fixed,
                                 random = params.rdm,
                                 data = moddata, 
                                 family = fam, 
                                 method = "REML",
                                 weights = weights,
                                 control=mgcv::gam.control(maxit = it))
          }
          else{ # correlation and weights != FALSE
            mod.cv <- mgcv::gamm(params.fixed,
                                 random = params.rdm,
                                 data = moddata, 
                                 family = fam, 
                                 method = "REML",
                                 weights = weights,
                                 correlation = correlation,
                                 control=mgcv::gam.control(maxit = it))
          }
          
          
          
          # make prediction with this model, but with the leftover data set
          testdata <- data[obj$samples[,i],]
          fit.cv = predict(mod.cv$gam, newdata = testdata, type = 'response', se.fit = FALSE)
          
          if(fam == "binomial"){ # if binomial, calculate Brier score
            # calculate Brier score Brier score for this prediction 
            # Bier score = mean((pred.prob-event_success_bin)^2)
            # following https://stackoverflow.com/questions/25149023/how-to-find-the-brier-score-of-a-logistic-regression-model-in-r
            bsrep <- c(bsrep, mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2))
          }
          else{ # else calculate RMSE
            # RSME = sqrt(mean((fitted-observed)^2))
            bsrep <- c(bsrep, sqrt(mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2)))
          }  
          
        } # end if(mod_func == "gamm")
        ##### glmer or lmer ######
        if(mod_func == "glmer" || mod_func == "lmer"){
          # run model
          if(fam == "gaussian" || mod_func == "lmer"){
            mod.cv <- lme4::lmer(paramFixed,
                                 data = moddata,
                                 control=lme4::lmerControl(optCtrl=list(maxfun = it)))
          }
          else{
            mod.cv <- lme4::glmer(paramFixed,
                                  data = moddata, family = fam,
                                  control=lme4::glmerControl(optCtrl=list(maxfun = it)))
          }
          # make prediction with this model, but with the leftover data set
          testdata <- data[obj$samples[,i],]
          #cat("predict")
          fit.cv = predict(mod.cv, newdata = testdata, allow.new.levels = TRUE, type = 'response', se.fit = FALSE)
          
          #### calculate and gather scores for k-folds #####
          if(fam == "binomial"){ # if binomial, calculate Brier score
            # calculate Brier score Brier score for this prediction 
            # Bier score = mean((pred.prob-event_success_bin)^2)
            # following https://stackoverflow.com/questions/25149023/how-to-find-the-brier-score-of-a-logistic-regression-model-in-r
            bsrep <- c(bsrep, mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2))
          }
          else{ # else calculate RMSE
            # RSME = sqrt(mean((fitted-observed)^2))
            bsrep <- c(bsrep, sqrt(mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2)))
          }  
        } # end if(mod_func == "glmer" || mod_func == "lmer")  
        ##### lm or glm ######
        if(mod_func == "glm" || mod_func == "lm"){
          # run model
          if(fam == "gaussian" || mod_func == "lm"){
            mod.cv <- lm(paramFixed,
                         data = moddata,
                         control = list(maxit = it, epsilon=1))
          }
          else{
            mod.cv <- glm(paramFixed,
                          data = moddata, family = fam,
                          control = list(maxit = it, epsilon=1))
          }
          # make prediction with this model, but with the leftover data set
          testdata <- data[obj$samples[,i],]
          cat("predict")
          fit.cv = predict(mod.cv, newdata = testdata, allow.new.levels = TRUE, type = 'response', se.fit = FALSE)
          
          #### calculate and gather scores for k-folds #####
          if(fam == "binomial"){ # if binomial, calculate Brier score
            # calculate Brier score Brier score for this prediction 
            # Bier score = mean((pred.prob-event_success_bin)^2)
            # following https://stackoverflow.com/questions/25149023/how-to-find-the-brier-score-of-a-logistic-regression-model-in-r
            bsrep <- c(bsrep, mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2))
          }
          else{ # else calculate RMSE
            # RSME = sqrt(mean((fitted-observed)^2))
            bsrep <- c(bsrep, sqrt(mean((fit.cv - unlist(data[obj$samples[,i],][, response]))^2)))
          }  
        } # end if(mod_func == "glmer" || mod_func == "lmer")  
      }# end i
      
      ### gather values of rept ####  
      # bind values together in data frame and rename added column
      bs <- cbind(bs, bsrep)
      colnames(bs)[which(colnames(bs) == "bsrep")] <- c(paste("rep", h, sep=""))
      
    }# end h
    
    bs$x <- NULL
    
    ##### create an object summarizing the performed CV and values #####
    summarycv <- list()
    summarycv$cvInfo <- paste(k,"-fold cross validation with ", rept, " repeats", sep="")
    summarycv$modInfo <- paste("Using model ", mod_func, " with family = ", fam, sep="")
    if(fam == "binomial"){
      summarycv$scoreInfo <- paste("For ", fam, " models, Brier score is used as criterion.", sep="")
    }
    else{
      summarycv$scoreInfo <- paste("For ", fam, " models, RMSE is used as criterion.", sep="")
    }
    statsum <- data.frame(repetition = colnames(bs), mean = NA, sd = NA)
    for(i in 1:rept){
      statsum$mean[i] <- mean(bs[, i])
      statsum$sd[i] <- sd(bs[, i])
    }
    summarycv$statSummary <- statsum
    summarycv$cvResults <- bs
    fullmean <- mean(t(summarycv$cvResults))
    fullsd <- sd(t(summarycv$cvResults))
    if(fam == "binomial"){
      summarycv$summaryInfo <- paste("Over all runs, the model had a mean acuracy of ", 
                                     round((1-(fullmean*2))*100, digits = 2),
                                     " % (Brier score of ",
                                     round(fullmean, digits = 4),
                                     ") and a standard deviation of ",
                                     round(fullsd*2*100, digits = 2),
                                     " % (SD of Brier score of ",
                                     round(fullsd, digits = 4),
                                     ")",
                                     sep="")
    }
    else{
      summarycv$summaryInfo <- paste("Over all runs, the model had a mean RMSE of ",
                                     round(fullmean, digits = 4),
                                     " and a standard deviation of RMSE of ",
                                     round(fullsd, digits = 4),
                                     ")",
                                     sep="")
    }
    
    cat(crayon::blue("100 % done!"))
    return(summarycv)
    
  } # if(mod_func == ...)
  else{
    
    cat(crayon::red("ERROR: The selected model function is not available or misspelled."))
    
  }
}



######## get_thresholds() ###########
#' created by Marvin Kiene, use this function at own risk!
#' @description
#' Calculate marginal effects (the derivative of the curve) to assess the change per 
#' concentration (the slope of the curve). 
#' So far, the automated version just works, if there is one significant decrease without any break!
#' Problems occur, if the uncertainty in the slope is very wiggly. You can check with margins::cplot().
#' 
#' @param data a data set to operate with. It can be a new predictor sequence for the model, not necessarily the
#' original data.
#' @param mod a model object to be used. Should be a binomial GLM.
#' @param eop 'effect of predictor': name the predictor that should be evaluated in ''.
#' @param acc the accuracy of the estimate. The higher the more accurate it will be estimated (takes more time).
#' @param lev the level where Ks should be calculated (default 0.5, but if another value is required, you can change).
#'
#' @return returns an object list with a one-row data frame incorporating the threshold values and their uncertainty 
#' and the data frame with the marginal effects and calculated values (mainly for plotting purpose).
#' @export
#'
get_thresholds <- function(mod, data, eop, acc = 1, lev = 0.5){
  # ms <- margins(mod.x, data = newdf) # could be used, but shows no confidence interval
  # this shows also confidence intervals
  R.devices::suppressGraphics({ # this wrapper function suppresses the graphics, thus increasing the performance
    cp <- margins::cplot(mod, eop, what = "effect", data = data, n = 1000*acc)
  })
  
  # check where there is a significant change. That is, where the CI does not touches 0
  cp$sig <- NA
  for(i in 1:nrow(cp)){
    if(cp$upper[i] < 0 || cp$lower[i] > 0){
      cp$sig[i] <- cp$yvals[i]
    }
  }
  
  # df to store the values later on
  thresholds <- data.frame(param = eop,
                           upper = NA,
                           # upper.sd = NA,
                           # upper.L = NA,
                           # upper.R = NA,
                           lower = NA,
                           # lower.sd = NA,
                           # lower.L = NA,
                           # lower.R = NA,
                           Ks = NA # bei lev = 0.5 der Response (hier 50 % Mortality)
                           # Ks.sd = NA,
                           # Ks.L = NA,
                           # Ks.R = NA
  )
  
  # row of apex of curve, necessary not to get two x values at a certain y value
  apex <- which(abs(cp$yvals) == max(abs(cp$yvals)), arr.ind=TRUE) # abs() as I don't know wether it is mortality (negative slope) or survival (positive slope)
  
  # get thresholds and uncertainty at the point where significant change starts/ends
  # Thus, I need to get the corresponding x values for the area of the CI 
  ## get the point where a significant change starts:
  for(i in 2: (nrow(cp)-1)){
    if(!is.na(cp$sig[i]) && is.na(cp$sig[i-1])){ # check the upper threshold point
      # upTh.y <- cp$yvals[i]
      # upTh.y.row <- i-1
      # upTh.y.confLy <- cp$upper[i] # for the upper threshold, the upper CI is left of the point and the lower CI right
      # upTh.y.confRy <- cp$lower[i]
      upTh.x <- cp$xvals[i] # this is the threshold value
    }
    if(is.na(cp$sig[i]) && !is.na(cp$sig[i-1])){ # check the lower threshold point
      # lwrTh.y <- cp$yvals[i]
      # lwrTh.y.row <- i
      # lwrTh.y.confLy <- cp$lower[i] # for the lower threshold, the lower CI is left of the point and the upper CI right
      # lwrTh.y.confRy <- cp$upper[i]
      lwrTh.x <- cp$xvals[i] # this is the threshold value
    }
  }
  
  # # uncertainty of upper threshold
  # lbr <- which.min(abs(cp$yvals[1:upTh.y.row] - upTh.y.confLy))
  # rbr <- upTh.y.row + which.min(abs(cp$yvals[upTh.y.row:apex] - upTh.y.confRy))
  # 
  # # get sd
  # upTh.sd <- sd(cp$xvals[c(lbr:rbr)])
  # upTh.xL <- cp$xvals[lbr]
  # upTh.xR <- cp$xvals[rbr]
  # 
  # # uncertainty of lower threshold
  # lbr <- apex + which.min(abs(cp$yvals[c(apex:lwrTh.y.row)] - lwrTh.y.confLy))
  # rbr <- lwrTh.y.row + which.min(abs(cp$yvals[lwrTh.y.row:nrow(cp)]) - lwrTh.y.confRy)
  # 
  # # get sd
  # lwrTh.sd <- sd(cp$xvals[c(lbr:rbr)])
  # lwrTh.xL <- cp$xvals[lbr]
  # lwrTh.xR <- cp$xvals[rbr]
  
  # gather values into df
  thresholds$upper[1] <- upTh.x
  # thresholds$upper.sd[1] <- upTh.sd
  # thresholds$upper.L[1] <- upTh.xL
  # thresholds$upper.R[1] <- upTh.xR
  
  thresholds$lower[1] <- lwrTh.x
  # thresholds$lower.sd[1] <- lwrTh.sd
  # thresholds$lower.L[1] <- lwrTh.xL
  # thresholds$lower.R[1] <- lwrTh.xR
  
  # now get Ks and uncertainty
  R.devices::suppressGraphics({ # this wrapper function suppresses the graphics, thus increasing the performance
    cp.p <- margins::cplot(mod, eop, what = "prediction", data = data, n = 1000*acc)
  })
  lineKs <- which.min(abs(lev-cp.p$yvals))
  Ks <- cp.p$xvals[lineKs]
  # Ks.yR <- cp$lower[lineKs]
  # Ks.yL <- cp$upper[lineKs]
  
  # # uncertainty of Ks
  # lbr <- which.min(abs(cp$yvals[1:lineKs] - Ks.yL))
  # rbr <- lineKs + which.min(abs(cp$yvals[lineKs:nrow(cp)] - Ks.yR))
  
  # # get sd
  # Ks.sd <- sd(cp$xvals[c(lbr:rbr)])
  # Ks.L <- cp$xvals[lbr]
  # Ks.R <- cp$xvals[rbr]
  
  thresholds$Ks[1] <- Ks
  # thresholds$Ks.sd[1] <- Ks.sd
  # thresholds$Ks.L[1] <- Ks.L
  # thresholds$Ks.R[1] <- Ks.R
  
  
  obj <- list(thresholds, cp)
  
  return(obj)
}



######### get_thresholds_gams() ##############
get_thresholds_gams <- function(data, Cholc, EPAc, mod, eval, poi){
newD <- expand.grid(Cholc = Cholc,
                      EPAc = EPAc)
  want <- seq(1, nrow(data), length.out = 200)
  #mod <- mod.y
  pdat <- newD
  
  p2 <- predict(mod, newdata = pdat, type = "terms", se.fit = TRUE)
  pdat <- transform(pdat, p2 = p2$fit[,eval], se2 = p2$se.fit[,eval]) # 1 = EPAc, 2 = Cholc, 3 = (EPAc,Cholc)
  
  # Aus dem Code vom Kurs -> wird nicht gebraucht, da wir ja keine Treatments haben
  # for(i in 1:length(pdat$Treatment)){
  #   if(pdat$Treatment[i] == "Triops"){
  #     pdat$p2[i] <- p2$fit[,3][i]
  #     pdat$se2[i] <- p2$se.fit[,3][i]
  #   }
  # }
  
  df.res <- df.residual(mod)
  crit.t <- qt(0.025, df.res, lower.tail = FALSE)
  pdat <- transform(pdat,
                    upper = p2 + (crit.t * se2),
                    lower = p2 - (crit.t * se2))
  # create Xp matrix with the difference of the point and the shifted point
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  eps <- 1e-7
  newD[,poi] <- newD[,poi] + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix") 
  Xp <- (X1 - X0) / eps
  
  # Aus dem Kurs: brauchst du hier nicht, da wir nur auf einen bestimmten Spline aus sind (jeweils)
  # Xp.r <- NROW(Xp)
  # Xp.c <- NCOL(Xp)
  
  ## which cols of xp relate to splines of interest?
  #c1 <- grepl('(Cholc)', colnames(Xp))
  #c2 <- grepl('(Cholc,EPAc)', colnames(Xp))
  ## which rows of xp relate to sites of interest?
  #r1 <- with(pdat, Treatment == '(Cholc)')
  #r2 <- with(pdat, Treatment == '(Cholc,EPAc)')
  
  ## get Xp matrix for control and Triops separately
  # Xp.control <- Xp[r1, ]
  
  # estimate derivatives 
  t.labs <- attr(mod$terms, "term.labels")
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  
  m4.d <- lD # here are the derivatives
  Term <- poi
  m4.dci <- confint(m4.d, term = Term)
  
  # check whether confidence interval touches 0
  m4.dsig <- signifD(pdat$p2, d = m4.d[[Term]]$deriv,
                     m4.dci[[Term]]$upper, m4.dci[[Term]]$lower)
  
  # sammeln der Daten der Ableitung
  test <- data.frame(poi = pdat[,poi],
                     deriv = m4.d[[poi]]$deriv,
                     upper = m4.dci[[poi]]$upper,
                     lower = m4.dci[[poi]]$lower,
                     sigincr = m4.d[[poi]]$deriv,
                     sigdecr = m4.d[[poi]]$deriv
  )
  test$sigincr[! c(!is.na(unlist(m4.dsig$incr)))] <- NA
  test$sigdecr[! c(!is.na(unlist(m4.dsig$decr)))] <- NA
  
  gathered.thresholds.sgr <- data.frame(param = poi,
                                        upper = NA,
                                        lower = NA,
                                        Ks.median = NA, # Ks in der Mitte des signifikanten Bereichs
                                        Ks.apex = NA # das ist der Ks am Wendepunkt der Kurve (bevorzugt, aber vllt nicht immer berechenbar)
  )
  
  gathered.thresholds.sgr$upper <- test[min( which( !is.na(test$sigincr) ) ),]$poi
  gathered.thresholds.sgr$lower <- test[max( which( !is.na(test$sigincr) ) )+1,]$poi
  gathered.thresholds.sgr$Ks.median <- test[round( median( which( !is.na(test$sigincr) ) ) ),]$poi
  #gathered.thresholds.sgr$Ks.apex <- test[which( test$deriv == max(test$deriv[ which( !is.na(test$sigincr) ) ]), arr.ind = TRUE ),]$poi
  
  obj <- list()
  obj$gathered.thresholds.sgr <- gathered.thresholds.sgr
  obj$data <- test
  
  return(obj)
}


