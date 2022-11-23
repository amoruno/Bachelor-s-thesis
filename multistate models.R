###############################
# Author: Aitor Moruno Cuenca #
# Date: 28JUN2022             #
###############################

# REQUIRE PACKAGES AND READ DATA

library(msm)
library(dplyr)
library(mstate)
library(ggplot2)
library(gplots)
library(gt)
library(lubridate)
library(gtsummary)
library(xtable)
webshot::install_phantomjs()

data <- read.csv("LifeTimeSurvivalTFG.csv", sep= ";", header = TRUE)

##########################################
#####         DATA CLEANING          #####
##########################################

# Generate states

# Function to create states
# State = 1, tuberculosis-free
# State = 2, infected of tuberculosis
# State = 3, dead.
###############################################################################

# If Death = NA or Positive = NA, we cannot determine state. State = 99 
# If Death = 1 and/or Positive = 1/0/NA, then State = 3.
# If Death = 0 and Positive = 1, then State = 2.
# If Death = 0 and Positive = 0/NA, then State = 1.

classify <- function(death, positive, euth){
  if(euth == 1){
    return(3)
  }else{
    if((!is.na(death)) && (!is.na(positive))){
      if(death == 1){
        return(4)
      }else{
        if(positive == 1){
          return(2)
        }
        return(1)
      }
    }
  }
  return(99)
}

data$event <- mapply(classify, data$dth, data$Val.Pos, data$euth)

# GENERATE VARIABLE Time.Block
# Remember time blocks:
# Block 1: September-November 2014
# Block 2: December 2014-February 2015
# Block 3: July-October 2015
# Block 4: January-March 2016
# Block 5: July-September 2016

#Convert to datetime
data$Sample.Date <- gsub("/", "", data$Sample.Date)
data$Sample.Date <- dmy(data$Sample.Date)

timeblock <- function(x){
  date <- "Out of time windows"
  if(year(x) == 2014){
      if(9 <= month(x) && month(x) <= 11){
        date <- "Block 1"
      }
      if(month(x) == 12){
        date <- "Block 2"
      }
  }
  if(year(x) == 2015){
      if(1 <= month(x) && month(x) <= 2){
        date <- "Block 2"
      }
      if(7 <= month(x) && month(x) <= 10){
        date <- "Block 3"
      }
  }
  if(year(x) == 2016){
    if(1 <= month(x) && month(x) <= 3){
        date <- "Block 4"
    }
    if(7 <= month(x) && month(x) <= 9){
        date <- "Block 5"
    }
  }
  return(date)
}

data$Time.block <- sapply(data$Sample.Date, timeblock)
raw_data <- data
# Check which subjects are out of window due to euthanasia
table(data %>% filter(Time.block == "Out of time windows") %>% select(euth))

# Impute missing block times via minimum distance with respect to time windows

min_block_distance <- function(x){
  reference <- data.frame(initial = as.Date(c("2014-09-30", "2014-12-31", 
                                              "2015-07-31", "2016-01-31", 
                                              "2016-07-31")), 
                          final = as.Date(c("2014-10-30", "2015-02-28", 
                                            "2015-10-31", "2016-03-31", 
                                            "2016-09-30")))
    min <- abs(as.Date(x)-reference[1,1])
    argfila_min <- 1
    for(i in 1:4){
      for(j in 1:2){
        dist <- abs(as.Date(x)-reference[i,j])
        if(dist < min){
          min <- dist
          argfila_min <- i
        }
      }
    }
  return(argfila_min)
}

rows <- as.numeric(rownames(data[data$Time.block == "Out of time windows",]))
for(e in rows){
  block <- min_block_distance(data$Sample.Date[e])
  data$Time.block[e] <- paste("Block", as.character(block), sep = " ")
}


# MSM needs at least two observations per id to run, hence we will delete
# all subjects with only one observation scheme, as it does not 
# offer any information

# Remove subjects with only one non-missing observed event

data <- data[data$event != 99,]
table_id <- table(data$id)

for(i in 1:length(table_id)){
  if(table_id[i] == 1){
    id_to_eliminate <- data[data$id == names(table_id[i]),]
    data <- data[!(data$id %in% c(id_to_eliminate)),]
  }
}

# Check for consistency as fixed time covariate, for each variable

check_fixed_time <- function(df, v){
  id <- unique(df$id)
  not_consistent <- c()
  for(e in id){
    aux <- df[df$id == e, v]
    t <- table(aux)
    if(nrow(t) >= 2){
      not_consistent <- c(not_consistent, e)
    }
  }
  return(not_consistent)
}

check_fixed_time(data, "sex")
check_fixed_time(data, "tx_grp")

# DESCRIPTIVE STATS AFTER AND BEFORE DATA CLEANING

summaries <- function(x, block){
  new_data <- x %>% subset(Time.block == block)
  new_data <- new_data[!duplicated(new_data, by = c(id, event)),]
  d<- new_data %>% select(sex, tx_grp, event) %>% 
                    tbl_summary(statistic = 
                                  c(sex, tx_grp, event) ~ "{n} ({p}%)")
  return(d)
}

d1_clean <- summaries(data, "Block 1")
d2_clean <- summaries(data, "Block 2")
d3_clean <- summaries(data, "Block 3")
d4_clean <- summaries(data, "Block 4")
d5_clean <- summaries(data, "Block 5")

merge_clean <-
  tbl_merge(
    tbls = list(d1_clean, d2_clean,d3_clean,d4_clean,d5_clean),
    tab_spanner = c("**Block 1**", "**Block 2**", "**Block 3**", 
                    "**Block 4**", "**Block 5**")
  ) %>% modify_footnote(everything() ~ "n/N (%)") %>% 
    modify_header(label = "**Variable**") %>%
    as_gt() %>%  tab_header(
    title = md("**Descriptive statistics**"),
    subtitle = md("*Based on LifeSurvivalDataset after data cleaning process*")
   )

gt::gtsave(merge_clean, file = file.path(getwd(), "clean_descriptive.png"))

d1_raw <- summaries(raw_data, "Block 1")
d2_raw <- summaries(raw_data, "Block 2")
d3_raw <- summaries(raw_data, "Block 3")
d4_raw <- summaries(raw_data, "Block 4")
d5_raw <- summaries(raw_data, "Block 5")
d6_raw <- summaries(raw_data, "Out of time windows")

merge_raw <-
  tbl_merge(
    tbls = list(d1_raw, d2_raw,d3_raw,d4_raw,d5_raw,d6_raw),
    tab_spanner = c("**Block 1**", "**Block 2**", "**Block 3**", "**Block 4**", 
                    "**Block 5**", "**Out of time windows**")
  ) %>% modify_footnote(everything() ~ "n/N (%)") %>% 
      modify_header(label = "**Variable**") %>% 
      as_gt() %>% tab_header(
      title = md("**Descriptive statistics**"),
      subtitle = md("*Based on raw dataset LifeSurvivalDataset*") 
    ) 

gt::gtsave(merge_raw, file = file.path(getwd(), "raw_descriptive.png"))

# Number of subjects at each time block

summaries_ind <- function(x){
  new_data <- x[!duplicated(x[c("id", "Time.block")]),]
  d<- new_data %>% select(Time.block) %>% 
    tbl_summary(statistic = 
                  c(Time.block) ~ "{n}")
  return(d)
}

subjects_clean <- summaries_ind(data)
subjects_raw <- summaries_ind(raw_data)

merge_clean_2 <-
  tbl_merge(
    tbls = list(subjects_clean, subjects_raw),
    tab_spanner = c("**Clean**", "**Raw**")
  ) %>% modify_footnote(everything() ~ "Number of subjects at time block") %>% 
  modify_header(label = "**Variable**", stat_0_1 = "Subjects", stat_0_2 = "Subjects") %>% 
  as_gt() %>% tab_header(
    title = md("**Number of subjects at each time block**"),
    subtitle = md("*Based on dataset LifeSurvivalDataset*") 
  ) 

gt::gtsave(merge_clean_2, file = file.path(getwd(), "subjects.png"))

# Order by id and start time

data <- data[with(data, order(id, start)),]
raw_data <- raw_data[with(raw_data, order(id, start)),]

balloon_plot <- function(x, path, name){
  states.descriptive <- statetable.msm(event, id, data=x)
  dt <- as.table(as.matrix(states.descriptive))
  dev.copy(png,paste0(path, name))
        balloonplot(t(dt), main = "",xlab ="   To  ", ylab="From
              ", colsrt = 1,
              label = TRUE, show.margins = TRUE, dotcolor = "grey", 
              scale.method = "volume", scale.range = "relative",
              colmar = 1, rowmar = 1.5)
  dev.off()
  return(1)
}

balloon_plot(data, getwd(), "/balloon_clean.jpg")

# Descriptive plots for treatment group
balloon_plot(data[data$tx_grp == 1,], getwd(), "/balloon_t1.png")
balloon_plot(data[data$tx_grp == 2,], getwd(), "/balloon_t2.png")
balloon_plot(data[data$tx_grp == 3,], getwd(), "/balloon_t3.png")


# FIT MODEL 1


# Which interval times correspond to each time block?

matrix <- matrix(data = NA, ncol = 2, nrow = 5)

for(i in 1:5){
  summary_start <- summary(data[data$Time.block == 
                                  paste("Block", as.character(i), sep = " "),
                                  "start"])
  summary_stop <- summary(data[data$Time.block == 
                                 paste("Block", as.character(i), sep = " "),
                               "stop"])
  matrix[i,1] <- summary_start[1]
  matrix[i,2] <- summary_stop[2]
}
matrix[1,1] <- 0

# Initialize transition intensity matrix
Q <- rbind(c(0, 1/3, 1/3, 1/3),
           c(1/3, 0, 1/3, 1/3),
           c(0, 0, 0, 0),
           c(0, 0, 0, 0)
)
Q.crude <- crudeinits.msm(event ~ start, id, data = data, qmatrix = Q)

# Fit model

# covariates as factors with reference level

data$sex <- relevel(factor(data$sex), ref = "F")
data$tx_grp <- relevel(factor(data$tx_grp), ref = 3)

data <- data[with(data, order(id, start)),]
model.msm <- msm(event ~ start, subject = id, exacttimes = TRUE,data = data, 
                 qmatrix = Q.crude, covariates = ~ tx_grp + sex, 
                 control = list(trace = 1, REPORT=1))

# MLE
model.msm

# Q matrix
q.matrix <- qmatrix.msm(model.msm)

# stratified q matrix by treatment

q.matrix.treatment <- list(qmatrix.msm(model.msm, covariates = list(tx_grp = 1)),
                           qmatrix.msm(model.msm, covariates = list(tx_grp = 2)),
                           qmatrix.msm(model.msm, covariates = list(tx_grp = 3))
                           )

# P matrix by time block
pmatrix.msm(model.msm, t = c(matrix[1,1], matrix[1,2]))
pmatrix.msm(model.msm, t = c(matrix[2,1], matrix[2,2]))
pmatrix.msm(model.msm, t = c(matrix[3,1], matrix[3,2]))
pmatrix.msm(model.msm, t = c(matrix[4,1], matrix[4,2]))
pmatrix.msm(model.msm, t = c(matrix[5,1], matrix[5,2]))

plot(model.msm)

# mean sojourn times

l1 <- cbind(data.frame(State = c("State 1", "State 2"),
                       data.frame(Treatment = "High-susceptibility"), 
                       sojourn.msm(model.msm, covariates = list(tx_grp = 1))))
l2 <- cbind(State = c("State 1", "State 2"),
            data.frame(Treatment = "High-contact"), 
            sojourn.msm(model.msm, covariates = list(tx_grp = 2)))
l3 <- cbind(State = c("State 1", "State 2"),
            data.frame(Treatment = "Control"), 
            sojourn.msm(model.msm, covariates = list(tx_grp = 3)))
df <- rbind(l1,l2,l3, make.row.names = FALSE)
df <- df[order(df$State),]
xtable(df)

# total length of stay

l1 <- totlos.msm(model.msm, covariates = list(tx_grp = 1))
l2 <- totlos.msm(model.msm, covariates = list(tx_grp = 2))
l3 <- totlos.msm(model.msm, covariates = list(tx_grp = 3))

df <- rbind(l1,l2,l3)
rownames(df) <- c("High-susceptibility", "High-contact", "Control")
xtable(df)

# hazard ratios

hazard.msm(model.msm, cl = 0.95)

# survival plots

# euthanasia
plot(model.msm, from = c(1,2),to = 3, legend.pos = c(0,0.5), 
     covariates = list(tx_grp = 1), range = c(0,600), 
     main = "High-susceptibility (group 1)")
plot(model.msm, from = c(1,2), to = 3, legend.pos = c(0,0.5), 
     covariates = list(tx_grp = 2), range = c(0,600), 
     main = "High-contact (group 2)")
plot(model.msm, from = c(1,2), to = 3, legend.pos = c(0,0.5), 
     covariates = list(tx_grp = 3), range = c(0,600), 
     main = "Control (group 3)")

#death
plot(model.msm, from = c(1,2),to = 4, legend.pos = c(0,0.5), 
     covariates = list(tx_grp = 1), range = c(0,600), 
     main = "High-susceptibility (group 1)")
plot(model.msm, from = c(1,2), to = 4, legend.pos = c(0,0.5), 
     covariates = list(tx_grp = 2), range = c(0,600), 
     main = "High-contact (group 2)")
plot(model.msm, from = c(1,2), to = 4, legend.pos = c(0,0.5), 
     covariates = list(tx_grp = 3), range = c(0,600), 
     main = "Control (group 3)")

# model assessment

plot.prevalence.msm(model.msm, mintime= 0, maxtime= 1000)


# FIT MODEL 2

data$event <- ifelse(data$event == 4, 3, data$event)
Q <- rbind(c(0, 1/2, 1/2),
           c(1/2, 0, 1/2),
           c(0, 0, 0))
Q.crude <- crudeinits.msm(event ~ start, id, data = data, qmatrix = Q)

model.msm <- msm(event ~ start, subject = id,  data = data, qmatrix = Q.crude, covariates = ~ tx_grp + sex, 
                 control = list(trace = 1, REPORT=1), obstype = 2)

# MLE
model.msm

# Q matrix
q.matrix <- qmatrix.msm(model.msm)

# stratified q matrix by treatment

q.matrix.treatment <- list(qmatrix.msm(model.msm, covariates = list(tx_grp = 1)),
                           qmatrix.msm(model.msm, covariates = list(tx_grp = 2)),
                           qmatrix.msm(model.msm, covariates = list(tx_grp = 3))
)

# P matrix by time block
pmatrix.msm(model.msm, t = c(matrix[1,1], matrix[1,2]))
pmatrix.msm(model.msm, t = c(matrix[2,1], matrix[2,2]))
pmatrix.msm(model.msm, t = c(matrix[3,1], matrix[3,2]))
pmatrix.msm(model.msm, t = c(matrix[4,1], matrix[4,2]))
pmatrix.msm(model.msm, t = c(matrix[5,1], matrix[5,2]))

plot(model.msm)

# mean sojourn times

l1 <- cbind(data.frame(State = c("State 1", "State 2"),
                       data.frame(Treatment = "High-susceptibility"), 
                       sojourn.msm(model.msm, covariates = list(tx_grp = 1))))
l2 <- cbind(State = c("State 1", "State 2"),
            data.frame(Treatment = "High-contact"), 
            sojourn.msm(model.msm, covariates = list(tx_grp = 2)))
l3 <- cbind(State = c("State 1", "State 2"),
            data.frame(Treatment = "Control"), 
            sojourn.msm(model.msm, covariates = list(tx_grp = 3)))
df <- rbind(l1,l2,l3, make.row.names = FALSE)
df <- df[order(df$State),]
xtable(df)

# total length of stay

l1 <- totlos.msm(model.msm, covariates = list(tx_grp = 1))
l2 <- totlos.msm(model.msm, covariates = list(tx_grp = 2))
l3 <- totlos.msm(model.msm, covariates = list(tx_grp = 3))

df <- rbind(l1,l2,l3)
rownames(df) <- c("High-susceptibility", "High-contact", "Control")
xtable(df)

# hazard ratios

hazard.msm(model.msm, cl = 0.95)

# survival plots

# to absorbing state
plot(model.msm, from = c(1,2),to = 3, legend.pos = c(0,0.42), 
     covariates = list(tx_grp = 1), range = c(0,700), 
     main = "High-susceptibility (group 1)")
plot(model.msm, from = c(1,2), to = 3, legend.pos = c(0,0.42), 
     covariates = list(tx_grp = 2), range = c(0,700), 
     main = "High-contact (group 2)"
plot(model.msm, from = c(1,2), to = 3, legend.pos = c(0,0.42), 
     covariates = list(tx_grp = 3), range = c(0,700), 
     main = "Control (group 3)")

# model assessment

plot.prevalence.msm(model.msm, mintime= 0, maxtime= 1000)



