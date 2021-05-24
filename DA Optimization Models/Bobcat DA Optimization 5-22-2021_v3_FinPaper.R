## Day Ahead Bidding Optimization under CVAR risk metric & ANN Classification Boost
## Date: 04-01-2021

rm(list=ls()) 
cat("\014")   # Clears the console (same as ctrl + L)

library(lubridate); library(roll)
library(xlsx); library(readxl); 
library(glpkAPI); library(Matrix)
# source("V:/Risk Management/Risk Management Book/R Functions/PeakType.R")
# source("V:/Risk Management/Risk Management Book/R Functions/DataAttributeGrid.R")
# source("V:/Risk Management/Risk Management Book/R Functions/XYZPivotQuantiles.R")
source("//SDHQFILE01.enxco.com/PowM$/Risk Management/Risk Management Book/R Functions/GetYesData.R")
source("V:/Risk Management/Risk Management Book/R Functions/ImportHedgeTerms.R")
# source("V:/Risk Management/Risk Management Book/R Functions/Deseasonalize.R")

## User inputs
startdate <- '2017-01-1'
enddate <- '2020-12-31'
AggLevel <- 'HOUR'
item_IDs <- 'ERCOT_SCED_GEN_60D_PLT:10002368788,RTLMP:10002236456,DALMP:10002236456,RTLMP:10000697078,DALMP:10000697078'
#Opex_Monthly_Cost <- -140000 ## enter as negative value (weekly)
#VAR_LIMIT <- -575000/30   ## Weekly Loss Limit (positive for loss)
VAR_LIMIT_MONTHLY <- -575000  # monthly loss limit
# VAR_QUANTILE <- .90
#VAR_LIMIT <- 5000   ## Weekly Loss Limit (positive for loss)
ZERO_HEDGE <- TRUE

# group.attributes <- c("MONTH", "HOURENDING")
ISO <- "ERCOT"
directory <- 'C:/Users/matthew.hunter/Dropbox/UW CFRM/601 Independent Study/Bobcat Bluff Day Ahead'
export.filename <- 'BCAT_ERCOTN_2017-Apr2021.csv'
ANN_prediction_filename <- 'ANN_PREDICTION_17-19.csv'
# attributes <-  c("WEEKNUM", "YEAR")
attributes <-  c("MARKETDAY", "YEAR")
Include_ANN_Prediction <- TRUE
DA_Prediction_Type <- "RF_Prediction_Node"
SD_Adder_Bids_NodalSale <- c(-3, -1.5, 0, 1.5, 3)  ## sigma adders
SD_Adder_Bids_HubBuy <- -SD_Adder_Bids_NodalSale   ## sigma adders
SD_Adder_Bids_PTP <- SD_Adder_Bids_HubBuy     ## sigma adders
Num_DA_Products <- 3 # DA node, PTP and Hub products

optim_startdate <- '2019-01-01'
optim_enddate <- '2019-12-31'
optim_startdate_test <- '2020-01-01'
optim_enddate_test <- '2020-12-31'

setwd(directory)

## USER-DEFINED FUNCTIONS ***********************
perc.rank <- function(x, xo)  length(x[x <= xo])/length(x)

# compute confusion matrix of PnL for actual and predictions for DA-RT spread
classification_pnl <- function(actual, pred, return.matrix = TRUE){
  tmp <- cbind(actual,pred)
  act_binary <- as.numeric(actual>0)
  pred_binary <- as.numeric(pred>0)
  
  results <- matrix(NA,2,2)
  colnames(results) <- c("Actual_TRUE", "Actual_FALSE")
  rownames(results) <- c("Pred_TRUE", "Pred_FALSE")
  
  results[1,1] <- sum(actual[(act_binary==TRUE) & (pred_binary==TRUE)])
  results[1,2] <- sum(actual[(act_binary==FALSE) & (pred_binary==TRUE)])
  results[2,1] <- sum(actual[(act_binary==TRUE) & (pred_binary==FALSE)])
  results[2,2] <- sum(actual[(act_binary==FALSE) & (pred_binary==FALSE)])
  
  total_accuracy <- abs(results[1,1]) + abs(results[2,2]) - abs(results[1,2]) - abs(results[2,1])
  
  if (return.matrix){return(results)} 
  else {return(total_accuracy)}
}

# Compute a rolling metric on vector x, but only conditioned on unique values of y
roll_cond_fun <- function(x, y, width, FUN=mean, new.col.name="Z.MEAN"){
  
  N <- length(x)
  unique.y <- unique(y)
  df <- as.data.frame(cbind(x,y))
  df$index <- seq(1,N)
  df$result <- NA
  
  for (n in (width+1):N){
    idx <- (df$index<n) & (df$index>(n-width)) & (df$y==df$y[n])
    df$result[n] <- FUN(df$x[idx])
  }
  
  return (df$result)
}

linearOpt <- function(Constraint_Matrix_Data, clower, cupper, LHSvec, RHSvec, Objvec, Minimize=TRUE, write_all_output=TRUE){
  
  ## initialize model
  lp <- initProbGLPK()
  
  ## Convert Constraint matrix to sparse formate for optimization
  matSparse <- as(Constraint_Matrix_Data, "ngCMatrix")
  
  # Build (row, col) pairs from non-zero entries
  matSparse_Summary <- as.data.frame(summary(matSparse))
  
  ## model dimensions
  nrows <- nrow(Constraint_Matrix_Data) 
  ncols <- ncol(Constraint_Matrix_Data)
  
  ## use sparse format for model data -- start with matrix of zeros and define non-zero values
  nnonzero <- nrow(matSparse_Summary)
  rindexnonzero <- matSparse_Summary$i
  cindexnonzero <- matSparse_Summary$j
  matSparse_Summary$valuesnonzero <- NA
  
  for (n in 1:nrow(matSparse_Summary)){
    matSparse_Summary$valuesnonzero[n] <- Constraint_Matrix_Data[matSparse_Summary$i[n], matSparse_Summary$j[n]]
  }
  valuesnonzero <- matSparse_Summary$valuesnonzero
  
  # maximize objective GLP_Max (minimize with GLP_MIN)
  if (Minimize) {
    setObjDirGLPK(lp, GLP_MIN)
  } else {
    setObjDirGLPK(lp, GLP_MAX)
  }
  
  # tell model how many rows and columns
  addRowsGLPK(lp, nrows)
  addColsGLPK(lp, ncols)
  
  # add column limits
  setColsBndsGLPK(lp,c(1:ncols), clower, cupper)
  setRowsBndsGLPK(lp,c(1:nrows),LHSvec,RHSvec)
  setObjCoefsGLPK(lp,c(1:ncols),Objvec)
  
  ## customization of variable types X are continuous, y's are binary
  setColsKindGLPK(lp,1:ncols,rep(1,ncols))
  
  # load constraint matrix coefficients
  loadMatrixGLPK(lp,nnonzero, rindexnonzero, cindexnonzero,valuesnonzero)
  
  # solve LP problem using Simplex Method
  solveSimplexGLPK(lp) ## ignores integer constraints
  #solveInterior(lp)
  #solveMIPGLPK(lp)   ##to solve integer model, must solve with Simplex first as starting points before integer constraints are applied cfrm 507: 11/10 lecture 0:06
  
  # get results of solution
  # solve status 5 = optimal solution found
  # getSolStatGLPK(lp)
  # status_codeGLPK(getSolStatGLPK(lp))
  
  # objective function value
  # getObjValGLPK(lp)
  
  # value of variables in optimal solution
  Solution.Variables <- as.data.frame(getColsPrimGLPK(lp))
  
  # status of each variable in optimal solution 1 = basic variable
  # getColsStatGLPK(lp)
  
  # get dual values for each row/constraint 
  # getRowsDualGLPK(lp)
  
  # get dual values for each variables // reduced costs
  # getColsDualGLPK(lp)
  
  # mipObjValGLPK(lp) 
  # printMIPGLPK(lp,fname="mipReport.txt")
  
  if (write_all_output){
    write.csv(Solution.Variables, "Solution_Variables.csv", row.names = TRUE)
    write.csv(getColsStatGLPK(lp), "Basic_Variables.csv", row.names = TRUE)
    write.csv(getRowsDualGLPK(lp), "Dual_Values_Constraint_Shadows.csv", row.names = TRUE)
    write.csv(getColsDualGLPK(lp), "Variable_Reduced_Costs.csv", row.names = TRUE)
  }
  
  # return (as.numeric(Solution.Variables))
  return (Solution.Variables)
}

make_constraints <- function(Optimization_Inputs_Mtx, VAR_QUANTILE, VAR_LIMIT, Bids_Segment_Num=5, export_constraint_data = FALSE, Include_ANN_Prediction = TRUE){
  
  N <- nrow(Optimization_Inputs_Mtx)
  Num_DA_Products <- 3  # Node, PTP and Hub DA products
  Num_DA_Combinations <- Num_DA_Products * Bids_Segment_Num * 2  ## x2 for ANN inclusion
  
  # Objective Function: DA position weights (0-1),  Project Rev (static), VAR position (1), E[k] slack values
  Objvec <- c(rep(0,Num_DA_Combinations), 0, VAR_LIMIT, rep(1/(N*(1-VAR_QUANTILE)),N))
  
  # Sum of each DA product segments must be >0 and <= 1
  for (j in 1:Num_DA_Products){
    tmp1 <- rep(0,length(Objvec))
    tmp1[((j-1)*Bids_Segment_Num*2+1):((j-1)*Bids_Segment_Num*2+Bids_Segment_Num*2)] <- 1 
    
    if (j==1) {
      constraint_3 <- tmp1
    }
    else {
      constraint_3 <- rbind(constraint_3, tmp1)
    }
  }
  
  # N constraints such that DA PnL + VAR LIMIT + E[k] >= 0 
  constraint_4 <- cbind(Optimization_Inputs_Mtx, diag(1,N,N))
  Constraint_Matrix_Data <- rbind(constraint_3, constraint_4)
  
  # setwd(directory)
  if (export_constraint_data){
    write.csv(Constraint_Matrix_Data,"ConstraintData.csv",row.names = FALSE)
  }
  
  # Limits on constraints, upper and lower bounds
  RHSvec <- c(rep(1,nrow(constraint_3)), rep(Inf,nrow(constraint_4))) #constraints upper bound 
  LHSvec <- c(rep(0,nrow(constraint_3)), rep(0, nrow(constraint_4))) #constraints lower bound
  # RHSvec <- c(rep(Bids_Segment_Num*2,nrow(constraint_3)), rep(Inf,nrow(constraint_4))) #if running naive and predictive combined optimization
  
  # column upper and lower bounds
  # decision variable non-negativity
  cupper <- c(rep(1,Num_DA_Combinations),1,1,rep(Inf,N))   # upper bounds on parameters
  clower <- c(rep(0,Num_DA_Combinations),1,1,rep(0,N))   # lower bounds on parameters
  
  # FORCE Optimization to choose between Naive and RF/ANN Strategy, for each DA product
  if (Include_ANN_Prediction==TRUE){
    for (j in 1:Num_DA_Products) {
      cupper[((j-1)*2*Bids_Segment_Num+1):((j-1)*2*Bids_Segment_Num+Bids_Segment_Num)] <- 0  ## FORCE NAIVE POSITION MAX = 0
    }
  } else {
    for (j in 1:Num_DA_Products) {
      cupper[(j*Bids_Segment_Num*2-Bids_Segment_Num+1):(j*Bids_Segment_Num*2)] <- 0  ## FORCE ANN POSITION MAX = 0
    }
  }
  
  return(list(Constraint_Matrix_Data, clower, cupper, LHSvec, RHSvec, Objvec))
}

emp_cvar_rng <-function(x, iters=10000, probs=c(.05,.5,.95), sample_per_iter = 365, limit){
  
  emp_cvar_sim <- matrix(NA,nrow=iters,ncol=1)
  for (i in 1:iters){
    s <- sample(x, size=sample_per_iter, replace = TRUE)
    cvar_s <- mean(limit - s[s<=limit])
    emp_cvar_sim[i,1] <- cvar_s 
  }
  return (quantile(emp_cvar_sim,probs,na.rm=TRUE))
}

emp_mean_rng <-function(x, iters=10000, probs=c(.05,.5,.95), sample_per_iter = 365){
  emp_sim <- matrix(NA,nrow=iters,ncol=1)
  for (i in 1:iters){
    s <- sample(x, size=sample_per_iter, replace = TRUE)
    emp_sim[i,1] <- mean(s)
  }
  return (quantile(emp_sim,probs,na.rm=TRUE))
}


compute_cvar <- function(data, x, limit){
  r <- data %*% x
  m <- mean(limit - r[r<=limit])
  return(m)
}

## Compute statistics on Attributes of a Timeseries 
AttributeFunction <- function(data, Price.Name, FUN, group.names = c("MONTH_NUM", "PEAKTYPE"),  return.attributes=TRUE){
  
  P <- length(Price.Name)
  
  data.lean <- data[,c(group.names)]
  data$GridLookup <- apply(data.lean,1,Paste) 
  
  group.seq <- MakeGroupSeq(data=data, group.names=group.names)
  Grid.Summary <- CreateGrid(group.seq=group.seq, group.names = group.names)
  
  # Build summary table for results
  # c <- 1
  Price.Data.Summary <- matrix(NA, nrow=nrow(Grid.Summary), ncol=P)
  
  for (n in 1:nrow(Grid.Summary)){
    data.tmp <- subset(data, data$GridLookup==Grid.Summary$GridLookup[n])
    
    for (p in 1:P){
      # error checks
      if (is.null(data[[Price.Name[p]]])) {return (NA)}  
      
      Price.Data.Summary[n,p] <- FUN(data.tmp[[Price.Name[p]]])
    }
  }
  
  Price.Data.Summary <- as.data.frame(cbind(Grid.Summary, Price.Data.Summary))
  
  colnames(Price.Data.Summary) <-  c(colnames(Grid.Summary), Price.Name)
  
  if (return.attributes){
    return (Price.Data.Summary) 
  } else {
    return (Price.Data.Summary[,match(Price.Name,colnames(Price.Data.Summary))]) 
  }
  
}

MakeGrid <- function(...){
  x <- list(...) # THIS WILL BE A LIST STORING EVERYTHING:
  expand.grid(...)      
}

MakeLookup <- function(...){
  x <- list(...) # THIS WILL BE A LIST STORING EVERYTHING:
  paste(..., sep="_")       
}

MakeList <- function(l){
  c_ <- ncol(l)
  tmp_ <- list()
  for (i in 1:c_){
    tmp_[[i]] <- l[[i]]
  }
  return (tmp_)
}

MakeGroupSeq <- function(data, group.names){
  group.num <- length(group.names)
  group.seq = list()
  
  # Create unique lists of attributes, and a grid of all combinations
  for (i in 1:group.num){
    tmp_ <- unique(data[[group.names[i]]])
    if (is.numeric(tmp_)){tmp_ <- tmp_[order(tmp_)]}
    group.seq[[i]] <- tmp_
  }
  return (group.seq)
}

Paste <- function(x) paste(na.omit(x),collapse="_")

CreateGrid <- function(group.seq, group.names){
  result_ <- MakeGrid(group.seq)
  colnames(result_) <- group.names
  result_$GridLookup <- apply(result_,1,Paste)
  return (result_)
}  

## LOAD DATA
## Import ANN prediction results
ann_pred_data <- read.csv(ANN_prediction_filename)

## Retrieve Data (Use HOUR or DAY AggLevel ONLY)
df <- getYesEnergyData(startdate=startdate, enddate=enddate, items=item_IDs, Interval='HOUR', complete_case=FALSE)
df$ISO_TIMEBLOCK <- return.peaktype(sim.datetimes.HE = df$DATETIME.POSIX , ISO=ISO)
df$WEEKNUM <- week(df$DATETIME.POSIX)
df$Reference <- seq(0,nrow(df)-1,1)

#rename price columns for generic use
colnames(df)[match("BCATWD_WD_1..RTLMP.",colnames(df))] <- "Price.Node.RT"
colnames(df)[match("HB_NORTH..RTLMP.",colnames(df))] <- "Price.Hub.RT"
colnames(df)[match("BCATWD_WD_1..DALMP.",colnames(df))] <- "Price.Node.DA"
colnames(df)[match("HB_NORTH..DALMP.",colnames(df))] <- "Price.Hub.DA"

## Add hedge terms to df
df <- importHedgeTerms(df, project.name = "Bobcat", time.convention = "HE1.24", append.hedge.term = TRUE, hedge.startdate.override="1/1/2017")

# Zero out hedge volume, if applicable (for DA nodal optimization only)
if (ZERO_HEDGE){df$Bobcat_Hedge_Volume <- 0}

## Import P50 Forecasted Generation from 3Tier API 
## Assumes 3Tier data is Hour Starting 0-23 local CPT, adjusts to become HE 1-24
## This can be automated in the future via Python integration
P50_Table <- read.csv("//SDHQFILE01.enxco.com/PowM$/Risk Management/Risk Management Book/Market Data/Vaisala Forecasts/Bobcat/historical_DA_forecasts_Apr2021Bobcat Bluff08-00-00.csv")
P50_Table$DATETIME.POSIX <- as.POSIXct(P50_Table$ForecastFlowDatetime, format="%Y-%m-%d %H:%M:%S")
P50_Table$DATETIME.POSIX.HE <- P50_Table$DATETIME.POSIX + hours(1)
P50_Table <- P50_Table[ , which(names(P50_Table) %in% c("DATETIME.POSIX.HE", "Power"))]
colnames(P50_Table)[match("DATETIME.POSIX.HE",colnames(P50_Table))] <- "DATETIME.POSIX"

## Import & merge results from ANN Prediction
df <- merge(df, ann_pred_data, by="Reference", all.x=TRUE)

# remove NAs
df$BCAT_DART[is.na(df$BCAT_DART)] <- 0
df$CatReg_Prediction_Node[is.na(df$CatReg_Prediction_Node)] <- 1
df$RF_Prediction_Node[is.na(df$RF_Prediction_Node)] <- 1
df$ANN_Prediction_Node[is.na(df$ANN_Prediction_Node)] <- 1

# Merge Forecasted P50 Generation Data into df 
df <- merge(df, P50_Table, by="DATETIME.POSIX")

# ********* START HERE: MANUAL LOAD DF data, instead of API pull and execution steps above ****************************
setwd(directory)
df <- read.csv("DF_Final.csv")
df$DATETIME.POSIX <- as.POSIXct(df$DATETIME.POSIX, format="%m/%d/%Y %H:%M")

# Compute Nodal Revenue and Hedge PnL
df$Rev_Node <- df$Price.Node.RT * df$Bobcat.Bluff.Wind.Project.LLC..ERCOT_SCED_GEN_60D_PLT.
df$PNL_Hedge <- df$Bobcat_Hedge_Volume * (df$Bobcat_Hedge_Price - df$Price.Hub.RT)

# Append volumes for DA Nodal Sale, PTP and Hedge Buyback
df$Volume_Max_NodalSale <- pmax(df$Power-df$Bobcat_Hedge_Volume,0)
df$Volume_Max_HubBuy <- pmax(df$Bobcat_Hedge_Volume-df$Power,0)  
df$Volume_Max_PTP <- pmin(df$Power,df$Bobcat_Hedge_Volume) 

df$Price.PTP.DA <- df$Price.Hub.DA-df$Price.Node.DA
df$Price.PTP.RT <- df$Price.Hub.RT-df$Price.Node.RT

df$PRICE_DART_NODE_ALL <- (df$Price.Node.DA - df$Price.Node.RT) 
df$DART_NODE_PNL_ALL <- df$PRICE_DART_NODE  * df$Volume_Max_NodalSale


## TESTING ON NODE ANN PREDICTION CONFUSION MATRIX
# df2 <- df[complete.cases(df),]
# optim_period <- (df2$DATETIME.POSIX >= optim_startdate) & (df2$DATETIME.POSIX <= optim_enddate)
# 
# classification_pnl(actual = df2$BCAT_DART[optim_period], pred = df2$ANN_Prediction_Node[optim_period], return.matrix = FALSE)
# classification_pnl(actual = df2$PRICE_DART_NODE[optim_period], pred = df2$ANN_Prediction_Node[optim_period], return.matrix = FALSE)
# 
# classification_pnl(actual = df2$BCAT_DART[optim_period], pred = df2$RF_Prediction_Node[optim_period], return.matrix = FALSE)
# classification_pnl(actual = df2$PRICE_DART_NODE[optim_period], pred = df2$RF_Prediction_Node[optim_period], return.matrix = FALSE)
# 
# classification_pnl(actual = df2$BCAT_DART[optim_period], pred = df2$CatReg_Prediction_Node[optim_period], return.matrix = FALSE)
# classification_pnl(actual = df2$PRICE_DART_NODE[optim_period], pred = df2$CatReg_Prediction_Node[optim_period], return.matrix = FALSE)


## Forecasted DA Price and st.dev -- Improve with timeseries ARMAX quantile
window_days <- 14
window_length <- 24*window_days
df$Price_NODE_FCST_RT_MEAN <- roll_cond_fun(x=df$Price.Node.DA, y = df$HOURENDING, width=window_length, FUN=mean) 
df$Price_NODE_FCST_RT_SD <- roll_cond_fun(x=df$Price.Node.DA, y = df$HOURENDING, width=window_length, FUN=sd) / sqrt(window_days) 

df$Price_PTP_FCST_RT_MEAN <- roll_cond_fun(x=df$Price.PTP.DA, y = df$HOURENDING, width=window_length, FUN=mean) 
df$Price_PTP_FCST_RT_SD <- roll_cond_fun(x=df$Price.PTP.DA, y = df$HOURENDING, width=window_length, FUN=sd) / sqrt(window_days) 

df$Price_HUB_FCST_RT_MEAN <- roll_cond_fun(x=df$Price.Hub.DA, y = df$HOURENDING, width=window_length, FUN=mean) 
df$Price_HUB_FCST_RT_SD <- roll_cond_fun(x=df$Price.Hub.DA, y = df$HOURENDING, width=window_length, FUN=sd) / sqrt(window_days) 

## remove NA observations
df <- df[complete.cases(df),]

## Create offer prices
# Sigma adders by DA transaction type

Bids_Segment_Num <- length(SD_Adder_Bids_PTP)
Bids_Segment_Perc <- rep(1,Bids_Segment_Num) ## allow each bid segment to reach 100% max, and 100% total across all

Price_Bids_Node <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)
Price_Bids_Hub <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)
Price_Bids_PTP <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)

Price_Bids_Node <- df$Price_NODE_FCST_RT_SD %*% t(SD_Adder_Bids_NodalSale) + df$Price_NODE_FCST_RT_MEAN
Price_Bids_Node <- cbind(Price_Bids_Node, Price_Bids_Node * df$DA_Classifier_Node)

Price_Bids_Hub <- df$Price_HUB_FCST_RT_SD %*% t(SD_Adder_Bids_HubBuy) + df$Price_HUB_FCST_RT_MEAN
Price_Bids_Hub <- cbind(Price_Bids_Hub, Price_Bids_Hub )  ## placeholder for ANN Classifier on Hub buybacks

Price_Bids_PTP <- df$Price_PTP_FCST_RT_SD %*% t(SD_Adder_Bids_PTP) + df$Price_PTP_FCST_RT_MEAN
Price_Bids_PTP <- cbind(Price_Bids_PTP, Price_Bids_PTP ) ## placeholder for PTP Classifier on Hub buybacks
 
## Rename columns
colnames(Price_Bids_Node) <- c(paste("Price_Bids_Node_Seg",seq(1,Bids_Segment_Num),sep=""),paste("Price_Bids_ANN_Node_Seg",seq(1,Bids_Segment_Num),sep=""))
colnames(Price_Bids_Hub) <- c(paste("Price_Bids_Hub_Seg",seq(1,Bids_Segment_Num),sep=""),paste("Price_Bids_ANN_Hub_Seg",seq(1,Bids_Segment_Num),sep=""))
colnames(Price_Bids_PTP) <- c(paste("Price_Bids_PTP_Seg",seq(1,Bids_Segment_Num),sep=""),paste("Price_Bids_ANN_PTP_Seg",seq(1,Bids_Segment_Num),sep=""))


add_DA_PNL <- function(df, Bids_Segment_Num, Bids_Segment_Perc){
  
  ## Compute awarded PnL
  PNL_AwardedBids_Node <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)
  PNL_AwardedBids_Hub <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)
  PNL_AwardedBids_PTP <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)
  
  ## DA Node PnL for Naive Strategy
  PNL_AwardedBids_Node[,1:Bids_Segment_Num] <- (Price_Bids_Node[,1:Bids_Segment_Num] <= df$Price.Node.DA) * (df$Price.Node.DA - df$Price.Node.RT) * (df$Volume_Max_NodalSale %*% t(Bids_Segment_Perc))
  
  ## DA PnL for ANN Adjust Strategy
  PNL_AwardedBids_Node[,(Bids_Segment_Num+1):(Bids_Segment_Num*2)] <- PNL_AwardedBids_Node[,1:Bids_Segment_Num] * df$DA_Classifier_Node * Include_ANN_Prediction
  
  ## DA Hub PnL for Naive Strategy
  PNL_AwardedBids_Hub[,1:Bids_Segment_Num] <- (Price_Bids_Hub[,1:Bids_Segment_Num] >= df$Price.Hub.DA) * (df$Price.Hub.RT - df$Price.Hub.DA) * (df$Volume_Max_HubBuy %*% t(Bids_Segment_Perc))
  
  ## DA Hub PnL for ANN Adjust Strategy
  PNL_AwardedBids_Hub[,(Bids_Segment_Num+1):(Bids_Segment_Num*2)] <- PNL_AwardedBids_Hub[,1:Bids_Segment_Num] * Include_ANN_Prediction
  
  ## DA PTP PnL for Naive Strategy
  PNL_AwardedBids_PTP[,1:Bids_Segment_Num] <- (Price_Bids_PTP[,1:Bids_Segment_Num] >= (df$Price.Hub.DA-df$Price.Node.DA)) * ((df$Price.Hub.RT-df$Price.Node.RT) - (df$Price.Hub.DA-df$Price.Node.DA)) * (df$Volume_Max_PTP %*% t(Bids_Segment_Perc))
  
  ## DA PTP PnL for ANN Adjust Strategy
  PNL_AwardedBids_PTP[,(Bids_Segment_Num+1):(Bids_Segment_Num*2)] <- PNL_AwardedBids_PTP[,1:Bids_Segment_Num] * Include_ANN_Prediction
  
  colnames_PNL_Node <- c(paste("PNL_Node_Seg",seq(1,Bids_Segment_Num),sep=""), paste("PNL_ANN_Node_Seg",seq(1,Bids_Segment_Num),sep=""))
  colnames_PNL_Hub <- c(paste("PNL_Hub_Seg",seq(1,Bids_Segment_Num),sep=""), paste("PNL_ANN_Hub_Seg",seq(1,Bids_Segment_Num),sep=""))
  colnames_PNL_PTP <- c(paste("PNL_PTP_Seg",seq(1,Bids_Segment_Num),sep=""), paste("PNL_ANN_PTP_Seg",seq(1,Bids_Segment_Num),sep=""))
  
  colnames(PNL_AwardedBids_Node) <- colnames_PNL_Node
  colnames(PNL_AwardedBids_Hub) <- colnames_PNL_Hub
  colnames(PNL_AwardedBids_PTP) <- colnames_PNL_PTP
  
  result <- list()
  
  df2 <- cbind(df, Price_Bids_Node, PNL_AwardedBids_Node, Price_Bids_Hub, PNL_AwardedBids_Hub, Price_Bids_PTP, PNL_AwardedBids_PTP)
  
  # summarize total Pnl by DA product (Node, Hub, PTP, All)
  df2$PNL_DA_Node <- rowSums(df2[,match(colnames_PNL_Node,colnames(df2))])
  df2$PNL_DA_Hub <- rowSums(df2[,match(colnames_PNL_Hub,colnames(df2))])
  df2$PNL_DA_PTP <- rowSums(df2[,match(colnames_PNL_PTP,colnames(df2))])
  df2$PNL_DA_Total <- df2$PNL_DA_Node + df2$PNL_DA_Hub + df2$PNL_DA_PTP
  df2$Project_Rev <- df2$Rev_Node + df2$PNL_Hedge
  
  result[[1]] <- df2
  result[[2]] <- c(colnames_PNL_Node, colnames_PNL_Hub, colnames_PNL_PTP, "Project_Rev")
  
  return(result)
  
}
  


## *************************************************************************************************

## Add look over attributes HOURENDING, MARKETDAY, WEEKNUM, MONTH_NUM

make_optimization_inputs_mtx <- function(optim_startdate, optim_enddate, df, VAR_LIMIT, attributes, mtx_col_names){
  
  optim_period <- (df$DATETIME.POSIX > optim_startdate) & (df$DATETIME.POSIX <= optim_enddate)
  df_grp <- AttributeFunction(data=df[optim_period,], 
                              Price.Name=mtx_col_names, 
                              FUN=sum, group.names = attributes,  return.attributes=TRUE)
  
  df_grp$VAR_LIMIT <- VAR_LIMIT 
  df_grp <- df_grp[!df_grp$Project_Rev==0,]  ## remove future months without nodal revenue
  
  result <- list()
  result[[1]]  <- as.matrix(df_grp[,!(colnames(df_grp) %in% c(attributes, "GridLookup"))]) # Optimization_Inputs_Mtx
  result[[2]]  <- df_grp # df_grp
  return (result)
}

# prepare results matrix for cvar analysis of all strategies and periods
attribute_list <- c("DATETIME.POSIX", "MARKETDAY", "WEEKNUM", "MONTH_NUM")
attribute_list_periods_per_month <- c(8760/12, 30, 4, 1)
attribute_list_periods_per_year <- c(8760, 365, 52, 12)

c.names <- c("Data Set", "ML Prediction Type", "Period", "CVAR RT Test Stat", "CVAR ML", "CVAR Naive", "CVAR RT", "MEAN NAIVE TESTSTAT", "MEAN ML", "MEAN Naive", "MEAN RT")
cvar_result <- data.frame(matrix(NA, nrow=2*length(DA_Prediction_Type_list)*length(attribute_list), ncol=length(c.names)))
colnames(cvar_result) <- c.names
test_stat_percentile <- 0.05

## Add loop over each DA prediction type ***************************************
DA_Prediction_Type_list <- c("RF_Prediction_Node", "ANN_Prediction_Node")
r1 = 0
for (i in 1:length(DA_Prediction_Type_list)){
  
  # update ml prediction type
  DA_Prediction_Type <- DA_Prediction_Type_list[i]
  df$DA_Classifier_Node <- 1  ## 1 will not adjust the naive DA strategy
  df$DA_Classifier_Node[df[DA_Prediction_Type]<0] <- 0 ## zero will withdraw DA offers based on ANN predictions
  
  # update DA Pnl with new ml prediction identifier
  tmp.pnl <- add_DA_PNL(df, Bids_Segment_Num, Bids_Segment_Perc)
  df2 <- tmp.pnl[[1]]
  mtx_col_names <- tmp.pnl[[2]]
  
  # Loop over optimization of different periods (daily, weekly, monthly)
  for (j in 1:length(attribute_list)){
    attributes <-  c(attribute_list[j], "YEAR")
    VAR_LIMIT <- VAR_LIMIT_MONTHLY / attribute_list_periods_per_month[j] 
    
    Optimization_Inputs_Mtx_Train <- make_optimization_inputs_mtx(optim_startdate, optim_enddate, 
                                                                  df2, VAR_LIMIT, attributes, mtx_col_names)[[1]]
    Optimization_Inputs_Mtx_Test <- make_optimization_inputs_mtx(optim_startdate_test, optim_enddate_test, 
                                                                 df2, VAR_LIMIT, attributes, mtx_col_names)[[1]]
    
    VAR_QUANTILE <- perc.rank(make_optimization_inputs_mtx(optim_startdate, optim_enddate, df2, VAR_LIMIT, attributes, mtx_col_names)[[2]]$Project_Rev, -VAR_LIMIT)
    
    constraints <- make_constraints(Optimization_Inputs_Mtx_Train, VAR_QUANTILE, VAR_LIMIT, 
                                    Bids_Segment_Num=Bids_Segment_Num, export_constraint_data = FALSE, 
                                    Include_ANN_Prediction = TRUE)
    
    # ALLOWS OPTIMIZATION TO CHOOSE ANY NAIVE OR rf/ann POSITIONS FOR EACH DA PRODUCT
    N <- nrow(Optimization_Inputs_Mtx_Train)
    Num_DA_Combinations <- Num_DA_Products * Bids_Segment_Num * 2  ## x2 for ANN inclusion
    cupper_naive <- c(rep(c(rep(1,Bids_Segment_Num), rep(0,Bids_Segment_Num)),time = Num_DA_Products),1,1,rep(Inf,N))   # upper bounds on parameters
    
    # decision parameters from optimization solutions
    sol_ml <- linearOpt(Constraint_Matrix_Data=constraints[[1]], clower=constraints[[2]], cupper=constraints[[3]], LHSvec=constraints[[4]], RHSvec=constraints[[5]], Objvec=constraints[[6]], Minimize=TRUE, write_all_output=FALSE)
    sol_naive <- linearOpt(Constraint_Matrix_Data=constraints[[1]], clower=constraints[[2]], cupper=cupper_naive, LHSvec=constraints[[4]], RHSvec=constraints[[5]], Objvec=constraints[[6]], Minimize=TRUE, write_all_output=FALSE)
    
    # Compute CVAR on Training Set solutions, don't include VAR column from optimization_input_mtx 
    C <- Num_DA_Combinations + 1
    
    mean_ml_train <- round(mean(Optimization_Inputs_Mtx_Train[,1:C] %*% sol_ml[1:C,]),2)
    mean_none_train <- round(mean( Optimization_Inputs_Mtx_Train[,1:C] %*% c(rep(0,Num_DA_Combinations),1) ),2)
    mean_naive_train <- round(mean( Optimization_Inputs_Mtx_Train[,1:C] %*% sol_naive[1:C,]),2)
    mean_naive_teststat_train <- round(as.numeric(emp_mean_rng(x = Optimization_Inputs_Mtx_Train[,1:C] %*% sol_naive[1:C,], probs=c(1-test_stat_percentile)), sample_per_iter = attribute_list_periods_per_year[j]),2)
    
    cvar_ml_train <- round(compute_cvar(data = Optimization_Inputs_Mtx_Train[,1:C], x = sol_ml[1:C,], limit=-VAR_LIMIT),2)
    cvar_none_train <- round(compute_cvar(data = Optimization_Inputs_Mtx_Train[,1:C], x = c(rep(0,Num_DA_Combinations),1), limit=-VAR_LIMIT),2)
    cvar_naive_train <- round(compute_cvar(data = Optimization_Inputs_Mtx_Train[,1:C], x = sol_naive[1:C,], limit=-VAR_LIMIT),2)
    
    # compute empirical CVAR on no DA strategy
    cvar_naive_teststat_train <- round(as.numeric(emp_cvar_rng(x = Optimization_Inputs_Mtx_Train[,1:C] %*% c(1,rep(0,Num_DA_Combinations-1),1), probs=c(test_stat_percentile), limit=-VAR_LIMIT), sample_per_iter = attribute_list_periods_per_year[j]),2)
    
    # Compute CVAR on Testing Set solutions, don't include VAR column from optimization_input_mtx 
    mean_ml_test <- round(mean(Optimization_Inputs_Mtx_Test[,1:C] %*% sol_ml[1:C,]),2)
    mean_none_test <- round(mean( Optimization_Inputs_Mtx_Test[,1:C] %*% c(rep(0,Num_DA_Combinations),1) ),2)
    mean_naive_test <- round(mean( Optimization_Inputs_Mtx_Test[,1:C] %*% sol_naive[1:C,]),2)
    
    mean_naive_teststat_test <- round(as.numeric(emp_mean_rng(x = Optimization_Inputs_Mtx_Test[,1:C] %*% sol_naive[1:C,], probs=c(1-test_stat_percentile)), sample_per_iter = attribute_list_periods_per_year[j]),2)
    
    cvar_ml_test <- round(compute_cvar(data = Optimization_Inputs_Mtx_Test[,1:C], x = sol_ml[1:C,], limit=-VAR_LIMIT),2)
    cvar_none_test <- round(compute_cvar(data = Optimization_Inputs_Mtx_Test[,1:C], x = c(rep(0,Num_DA_Combinations),1), limit=-VAR_LIMIT),2)
    cvar_naive_test <- round(compute_cvar(data = Optimization_Inputs_Mtx_Test[,1:C], x = sol_naive[1:C,], limit=-VAR_LIMIT),2)
    
    ## compute empirical CVAR on no DA strategy
    cvar_naive_teststat_test <- round(as.numeric(emp_cvar_rng(x = Optimization_Inputs_Mtx_Test[,1:C] %*% c(1,rep(0,Num_DA_Combinations-1),1), probs=c(test_stat_percentile), limit=-VAR_LIMIT), sample_per_iter = attribute_list_periods_per_year[j]),2)
    
    # store results for training & testing period in optimization
    r1 <- r1+1
    r2 <- r1 + length(attribute_list)*length(attribute_list)
    cvar_result[r1,] <- c("Training Data", DA_Prediction_Type, attributes[1], 
                          cvar_naive_teststat_train, cvar_ml_train, cvar_naive_train, cvar_none_train,
                          mean_naive_teststat_train, mean_ml_test, mean_naive_test, mean_none_test)
    cvar_result[r2,] <- c("Testing Data", DA_Prediction_Type, attributes[1], 
                          cvar_naive_teststat_test, cvar_ml_test, cvar_naive_test, cvar_none_test,
                          mean_naive_teststat_test, mean_ml_test, mean_naive_test, mean_none_test)

  }
  
}

cvar_result <- cvar_result[complete.cases(cvar_result),]
write.csv(cvar_result, "cvar_result.csv", row.names = FALSE)



# Return Naive Revenue Series for comparison to Normal
df$DA_Classifier_Node <- 1
tmp.pnl <- add_DA_PNL(df, Bids_Segment_Num, Bids_Segment_Perc)
df2 <- tmp.pnl[[1]]
mtx_col_names <- tmp.pnl[[2]]
rev_naive <- list()

for (j in 1:length(attribute_list)){
  attributes <-  c(attribute_list[j], "YEAR")
  VAR_LIMIT <- VAR_LIMIT_MONTHLY / attribute_list_periods_per_month[j] 
  
  Optimization_Inputs_Mtx_Test <- make_optimization_inputs_mtx(optim_startdate_test, optim_enddate_test, 
                             df2, VAR_LIMIT, attributes, mtx_col_names)[[1]]
  
  rev_naive[[j]] <- Optimization_Inputs_Mtx_Test[,1:C] %*% c(1,rep(0,Num_DA_Combinations-1),1)
}

rev_naive[[1]]


library(stats)

j <- 3
attribute_list[j]
x <- rev_naive[[j]] 
x_mu <- mean(x)
x_sd <- sd(x)
x.norm <- (x-x_mu)/x_sd

par(mfrow=c(1,2))
qqnorm(x.norm, datax = TRUE,
       xlab = "normal quantile", # this is for the theoretical quantiles
       ylab = "normalized naive revenue", # this is for the sample quantiles
       main = "normal qq probability plot", ylim=c(-10,30))
qqline(x.norm, datax=TRUE, col = 2)
d <- density(x.norm, adjust = 1, na.rm = TRUE)
plot(d, type = "n", xlim=c(-20,20), main="KDE for Normalized Naive Revenue 
and N(0,1) density") # more sophisticated plot
polygon(d, col = "wheat")
z = seq(from=-20,to=20,by=.01)
lines(z,dnorm(z), lty=2,lwd=3,col="red")



