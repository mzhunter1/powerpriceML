## Day Ahead Bidding Optimization under CVAR risk metric & ANN Classification Boost
## Date: 04-01-2021

rm(list=ls()) 
cat("\014")   # Clears the console (same as ctrl + L)

library(lubridate); library(roll)
library(xlsx); library(readxl); 
library(glpkAPI); library(Matrix)
source("V:/Risk Management/Risk Management Book/R Functions/PeakType.R")
source("V:/Risk Management/Risk Management Book/R Functions/DataAttributeGrid.R")
source("V:/Risk Management/Risk Management Book/R Functions/XYZPivotQuantiles.R")
source("//SDHQFILE01.enxco.com/PowM$/Risk Management/Risk Management Book/R Functions/GetYesData.R")
source("V:/Risk Management/Risk Management Book/R Functions/ImportHedgeTerms.R")
source("V:/Risk Management/Risk Management Book/R Functions/Deseasonalize.R")

## User inputs
startdate <- '2017-01-1'
enddate <- '2020-12-31'
AggLevel <- 'HOUR'
item_IDs <- 'ERCOT_SCED_GEN_60D_PLT:10002368788,RTLMP:10002236456,DALMP:10002236456,RTLMP:10000697078,DALMP:10000697078'
#Opex_Monthly_Cost <- -140000 ## enter as negative value (weekly)
VAR_LIMIT <- -575000/30   ## Weekly Loss Limit (positive for loss)
VAR_QUANTILE <- .90
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
setwd(directory)

optim_startdate <- '2019-01-01'
optim_enddate <- '2019-12-31'

optim_startdate_test <- '2020-01-01'
optim_enddate_test <- '2020-12-31'

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

emp_cvar_rng <-function(x, iters=10000, probs=c(.05,.5,.95), sample_per_iter = 365, limit){
  
  emp_cvar_sim <- matrix(NA,nrow=iters,ncol=1)
  for (i in 1:iters){
    s <- sample(x, size=sample_per_iter, replace = TRUE)
    cvar_s <- mean(limit - s[s<=limit])
    emp_cvar_sim[i,1] <- cvar_s 
  }
  return (quantile(emp_cvar_sim,probs))
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

# df$ANN_Prediction_Node[is.na(df$ANN_Prediction_Node)] <- 1
df$DA_Classifier_Node <- 1  ## 1 will not adjust the naive DA strategy
df$DA_Classifier_Node[df[DA_Prediction_Type]<0] <- 0 ## zero will withdraw DA offers based on ANN predictions

# Merge Forecasted P50 Generation Data into df 
df <- merge(df, P50_Table, by="DATETIME.POSIX")

# ********* START HERE: MANUAL LOAD DF data, instead of API pull and execution steps above ****************************
setwd(directory)
df <- read.csv("DF_Final.csv")

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
SD_Adder_Bids_NodalSale <- c(-3, -1.5, 0, 1.5, 3)  ## sigma adders
SD_Adder_Bids_HubBuy <- c(3, 1.5, 0, -1.5, -3)   ## sigma adders
SD_Adder_Bids_PTP <- SD_Adder_Bids_HubBuy     ## sigma adders
Bids_Segment_Num <- length(SD_Adder_Bids_PTP)
Bids_Segment_Perc <- rep(1,Bids_Segment_Num) ## allow each bid segment to reach 100% max, and 100% total across all

Price_Bids_Node <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)
Price_Bids_Hub <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)
Price_Bids_PTP <- matrix(NA, nrow=nrow(df), ncol=Bids_Segment_Num*2)

Price_Bids_Node <- df$Price_NODE_FCST_RT_SD %*% t(SD_Adder_Bids_NodalSale) + df$Price_NODE_FCST_RT_MEAN
Price_Bids_Node <- cbind(Price_Bids_Node, Price_Bids_Node * df$ANN_Classifier_Node)

Price_Bids_Hub <- df$Price_HUB_FCST_RT_SD %*% t(SD_Adder_Bids_HubBuy) + df$Price_HUB_FCST_RT_MEAN
Price_Bids_Hub <- cbind(Price_Bids_Hub, Price_Bids_Hub )  ## placeholder for ANN Classifier on Hub buybacks

Price_Bids_PTP <- df$Price_PTP_FCST_RT_SD %*% t(SD_Adder_Bids_PTP) + df$Price_PTP_FCST_RT_MEAN
Price_Bids_PTP <- cbind(Price_Bids_PTP, Price_Bids_PTP ) ## placeholder for PTP Classifier on Hub buybacks
 
## Rename columns
colnames(Price_Bids_Node) <- c(paste("Price_Bids_Node_Seg",seq(1,Bids_Segment_Num),sep=""),paste("Price_Bids_ANN_Node_Seg",seq(1,Bids_Segment_Num),sep=""))
colnames(Price_Bids_Hub) <- c(paste("Price_Bids_Hub_Seg",seq(1,Bids_Segment_Num),sep=""),paste("Price_Bids_ANN_Hub_Seg",seq(1,Bids_Segment_Num),sep=""))
colnames(Price_Bids_PTP) <- c(paste("Price_Bids_PTP_Seg",seq(1,Bids_Segment_Num),sep=""),paste("Price_Bids_ANN_PTP_Seg",seq(1,Bids_Segment_Num),sep=""))


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

df <- cbind(df, Price_Bids_Node, PNL_AwardedBids_Node, Price_Bids_Hub, PNL_AwardedBids_Hub, Price_Bids_PTP, PNL_AwardedBids_PTP)

# summarize total Pnl by DA product (Node, Hub, PTP, All)
df$PNL_DA_Node <- rowSums(df[,match(colnames_PNL_Node,colnames(df))])
df$PNL_DA_Hub <- rowSums(df[,match(colnames_PNL_Hub,colnames(df))])
df$PNL_DA_PTP <- rowSums(df[,match(colnames_PNL_PTP,colnames(df))])
df$PNL_DA_Total <- df$PNL_DA_Node + df$PNL_DA_Hub + df$PNL_DA_PTP

df$Project_Rev <- df$Rev_Node + df$PNL_Hedge

# write.csv(df,"df.csv", row.names = FALSE)

## Summarize results by month-year in new df 
optim_period <- (df$DATETIME.POSIX > optim_startdate) & (df$DATETIME.POSIX <= optim_enddate)
df_grp <- AttributeFunction(data=df[optim_period,], Price.Name=c(colnames_PNL_Node, colnames_PNL_Hub, colnames_PNL_PTP, "Project_Rev"), FUN=sum, group.names = attributes,  return.attributes=TRUE)
quantile(df_grp$Project_Rev,.1)

#df_grp$Opex <- Opex_Monthly_Cost  ## Avg. Monthly Operating Expense
df_grp$VAR_LIMIT <- VAR_LIMIT 
df_grp <- df_grp[!df_grp$Project_Rev==0,]  ## remove future months without nodal revenue


## Prepare summary data for optimization
Optimization_Inputs_Mtx <- as.matrix(df_grp[,!(colnames(df_grp) %in% c(attributes, "GridLookup"))])

## OPTIMIZATION ***********************************************************************************
## CVAR

setwd(directory)
write.csv(df_grp,"df_grp_rf_daily_nodeonly_17-19.csv",row.names = FALSE)
#write.csv(Optimization_Inputs_Mtx,"Optimization_Data.csv",row.names = FALSE)

N <- nrow(Optimization_Inputs_Mtx)
Num_DA_Products <- 3  # Node, PTP and Hub DA products
Num_DA_Combinations <- Num_DA_Products * Bids_Segment_Num * 2  ## x2 for ANN inclusion

# Objective Function: DA position weights (0-1),  Project Rev (static), VAR position (1), E[k] slack values
VAR_QUANTILE <- perc.rank(df_grp$Project_Rev, -VAR_LIMIT)
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
write.csv(Constraint_Matrix_Data,"ConstraintData.csv",row.names = FALSE)

# Limits on constraints, upper and lower bounds
RHSvec <- c(rep(1,nrow(constraint_3)), rep(Inf,nrow(constraint_4))) #constraints upper bound 
LHSvec <- c(rep(0,nrow(constraint_3)), rep(0, nrow(constraint_4))) #constraints lower bound
# RHSvec <- c(rep(Bids_Segment_Num*2,nrow(constraint_3)), rep(Inf,nrow(constraint_4))) #if running naive and predictive combined optimization

# column upper and lower bounds
# decision variable non-negativity
cupper <- c(rep(1,Num_DA_Combinations),1,1,rep(Inf,N))   # upper bounds on parameters
clower <- c(rep(0,Num_DA_Combinations),1,1,rep(0,N))   # lower bounds on parameters

# ALLOWS OPTIMIZATION TO CHOOSE ANY NAIVE OR rf/ann POSITIONS FOR EACH DA PRODUCT
cupper_naive <- c(rep(c(rep(1,Bids_Segment_Num), rep(0,Bids_Segment_Num)),time = Num_DA_Products),1,1,rep(Inf,N))   # upper bounds on parameters

# FORCE Optimization to choose between Naive and RF/ANN Strategy, for each DA product
if (Include_ANN_Prediction==TRUE){
  for (j in 1:Num_DA_Products) {
    cupper[((j-1)*2*Bids_Segment_Num+1):((j-1)*2*Bids_Segment_Num+Bids_Segment_Num)] <- 0  ## FORCE NAIVE POSITION MAX = 0
  }
} else {
  for (j in 1:Num_DA_Products) {
    cupper[(j*Bids_Segment_Num*2-Bids_Segment_Num+1):(j*Bids_Segment_Num*2)] <- 0  ## FORCE NAIVE POSITION MAX = 0
  }
}

sol1 <- linearOpt(Constraint_Matrix_Data, clower, cupper, LHSvec, RHSvec, Objvec, Minimize=TRUE, write_all_output=FALSE)
sol2 <- linearOpt(Constraint_Matrix_Data, clower, cupper_naive, LHSvec, RHSvec, Objvec, Minimize=TRUE, write_all_output=FALSE)

# Compute CVAR on solutions, don't include VAR column from optimization_input_mtx 
C <- Num_DA_Combinations+1

compute_cvar <- function(data, x, limit){
  r <- data %*% x
  m <- mean(limit - r[r<=limit])
  return(m)
}


cvar_rf <- compute_cvar(data = Optimization_Inputs_Mtx[,1:C], x = sol1[1:C,], limit=-VAR_LIMIT)
cvar_none <- compute_cvar(data = Optimization_Inputs_Mtx[,1:C], x = c(rep(0,Num_DA_Combinations),1), limit=-VAR_LIMIT)
cvar_naive <- compute_cvar(data = Optimization_Inputs_Mtx[,1:C], x = c(1,rep(0,Num_DA_Combinations-1),1), limit=-VAR_LIMIT)


## compute empirical CVAR on no DA strategy
cvar_naive_rng <- emp_cvar_rng(x = Optimization_Inputs_Mtx[,1:C] %*% c(1,rep(0,Num_DA_Combinations-1),1), limit=-VAR_LIMIT)




## Testing period optimization results

# Summarize results by month-year in new df 
optim_period_test <- (df$DATETIME.POSIX > optim_startdate_test) & (df$DATETIME.POSIX <= optim_enddate_test)
df_grp_test <- AttributeFunction(data=df[optim_period_test,], Price.Name=c(colnames_PNL_Node, colnames_PNL_Hub, colnames_PNL_PTP, "Project_Rev"), FUN=sum, group.names = attributes,  return.attributes=TRUE)
df_grp_test$VAR_LIMIT <- VAR_LIMIT 
df_grp_test <- df_grp_test[!df_grp_test$Project_Rev==0,]  ## remove future months without nodal revenue

## Prepare summary data for optimization
Optimization_Inputs_Mtx <- as.matrix(df_grp_test[,!(colnames(df_grp_test) %in% c(attributes, "GridLookup", "Opex"))])

cvar_rf_test <- compute_cvar(data = Optimization_Inputs_Mtx[,1:C], x = sol1[1:C,], limit=-VAR_LIMIT)
cvar_none_test <- compute_cvar(data = Optimization_Inputs_Mtx[,1:C], x = c(rep(0,Num_DA_Combinations),1), limit=-VAR_LIMIT)
cvar_naive_test <- compute_cvar(data = Optimization_Inputs_Mtx[,1:C], x = c(1,rep(0,Num_DA_Combinations-1),1), limit=-VAR_LIMIT)
 

## compute empirical CVAR on no DA strategy
cvar_naive_rng_test <- emp_cvar_rng(x = Optimization_Inputs_Mtx[,1:C] %*% c(1,rep(0,Num_DA_Combinations-1),1), limit=-VAR_LIMIT)



setwd(directory)
write.csv(df_grp_test,"df_grp_test_rf_daily_nodeonly.csv",row.names = FALSE)













## OPTIMIZATION ***********************************************************************************
## Mean - Variance optimization of monthly gross margin

library(quadprog)

## Convert all revenue int0 $ Millions
Optimization_Inputs_Mtx <- round(Optimization_Inputs_Mtx / (abs(Opex_Monthly_Cost)),4)

## Decision Variables = [DA_Node 1- 5, DA_Hub 1-5, DA PTP 1- 5, Nodal Rev, Hedge PnL, Opex]
N <- ncol(Optimization_Inputs_Mtx)
# num_fixed_parameters <- sum((colnames(Optimization_Inputs_Mtx) %in% c("Rev_Node", "PNL_Hedge")))
num_fixed_parameters <- 0
num_var_parameters <- N - num_fixed_parameters

# Expected revenue by column
dvec <- apply(Optimization_Inputs_Mtx,MARGIN = 2, FUN=mean)

# Covariance matrix with gamma penalty for risk aversion
gamma <- 1  #scalar for risk tolerance. Increase gamma to penalize variance more
# gamma <- gamma2
cov_mtx <- cov(Optimization_Inputs_Mtx)
C <- chol(cov_mtx)
Dmat <- t(C) %*% diag(sqrt(gamma),N) %*% C

## Constraints: A^t b >= b_0
# constraint coefficients Amat
c1 <- cbind(matrix(0,num_fixed_parameters,num_var_parameters),diag(1,num_fixed_parameters,num_fixed_parameters)) # nodal rev and hedge pnl = 1, fixed parameters
c2 <- cbind(diag(1,num_var_parameters,num_var_parameters), matrix(0,num_var_parameters,num_fixed_parameters)) # non-neg position weight >= 0
c3 <- cbind(diag(-1,num_var_parameters,num_var_parameters), matrix(0,num_var_parameters,num_fixed_parameters)) # position weight <= 1
Amat <- rbind(c1,c2,c3) 

# constraint limit values: b_0
bvec <- as.vector(c(rep(1,num_fixed_parameters),rep(0,num_var_parameters),rep(-1,num_var_parameters)))

## Solve quadratic program
sol <- solve.QP(Dmat*2, dvec, t(Amat), bvec=bvec, meq=2, factorized=FALSE)
results <- round(sol$solution,3)
results

## mean result
res.mean <- dvec %*% results
res.var <- t(results) %*% Dmat %*% results
res.gamma <- res.var / abs(res.mean)  


all_strat <- rep(1,N)
all.mean <- dvec %*% all_strat
all.var <- t(all_strat) %*% Dmat %*% all_strat
all.gamma <- all.var / abs(all.mean) 
all.gamma

NumPorts = 100
portraverse = NULL
portrtol = NULL
portweights = array(0,c(NumPorts,N))
portreturn = NULL
portcov = NULL
objval = NULL
portlagr = array(0,c(NumPorts,N+1))

for(nport in 1:NumPorts)
{
  #raverse =  0.2 + nport * 0.025
  gamma <- nport*50
  
  Dmat <- t(C) %*% diag(sqrt(gamma),N) %*% C
  solution <- solve.QP(Dmat*2, dvec, t(Amat), bvec=bvec, meq=2, factorized=FALSE)
  
  portweights[nport,] <- solution$solution
  portreturn <- c(portreturn,sum(solution$solution * dvec))
  objval <- c(objval,solution$value)
  portcov <- c(portcov,solution$value + sum(solution$solution * dvec))
 
}

portweights <- round(portweights,4)


# draw the efficient frontier
plot(sqrt(portcov), portreturn, type="l",
     xaxt = "n", yaxt = "n", frame.plot = FALSE)
axis(1, pos=0)
axis(2, pos=0)



## Linear Optimization of flattened monthly revenue

N <- nrow(Optimization_Inputs_Mtx)
Optimization_Inputs_Vec <- c(Optimization_Inputs_Mtx)
Num_DA_Products <- 3  # Node, PTP and Hub DA products

# write.csv(Optimization_Inputs_Mtx,"OptConstMtx.csv",row.names = FALSE)
# write.csv(Optimization_Inputs_Vec,"OptConstVec.csv",row.names = FALSE)

# develop constraints for linear optimization
num_fixed_parameters <- sum((colnames(Optimization_Inputs_Mtx) %in% c("Rev_Node", "PNL_Hedge", "Opex")))
tmp1 <- cbind(rep(1,N-1),diag(x=-1, nrow=N-1, ncol=N-1))
constraints_1 <- matrix(0,(N-1)*Num_DA_Products*Bids_Segment_Num, N*Num_DA_Products*Bids_Segment_Num)  # constrain all months of same DA product type and bid segment to be solved as the equivalent participation

for (j in 1:(Bids_Segment_Num*Num_DA_Products)){
  # r1 = (j-1)*N+1
  r2 = j*N
  constraints_1 [((j-1)*(N-1)+1):(j*(N-1)),(r2-N+1):r2] <- tmp1
}

zeros_Nx3 <- matrix(0,nrow(constraints_1),num_fixed_parameters*N)

# All month-year positions must be the same for each product type, and must include 
# tmp.lvl1 <- cbind(tmp2, Z, Z, zeros_Nx3) 
# tmp.lvl2 <- cbind(Z, tmp2, Z, zeros_Nx3) 
# tmp.lvl3 <- cbind(Z, Z, tmp2, zeros_Nx3) 
Constraint_Matrix_Data <- cbind(constraints_1, zeros_Nx3)

setwd(directory)
write.csv(Constraint_Matrix_Data,"ConstraintData.csv",row.names = FALSE)


## Add last three variables to include RT Nodal Rev, Hedge PnL and Opex
# Objvec <- rep(1,N*Num_DA_Products*Bids_Segment_Num + num_fixed_parameters*N)
Objvec <- Optimization_Inputs_Vec
# varnames <- colnames(Data)

RHSvec <- c(rep(0,nrow(constraints_1))) #constraints upper bound (zero, all equivalent)
LHSvec <- c(rep(0,nrow(constraints_1))) #constraints lower bound (zero, all equivalent)

# column upper and lower bounds
# decision variable non-negativity
cupper <- c(rep(1,N*Num_DA_Products*Bids_Segment_Num),rep(1,num_fixed_parameters*N))   # upper bounds on parameters
clower <- c(rep(0,N*Num_DA_Products*Bids_Segment_Num),rep(1,num_fixed_parameters*N))   # lower bounds on parameters


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
setObjDirGLPK(lp, GLP_MAX)

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
getSolStatGLPK(lp)
status_codeGLPK(getSolStatGLPK(lp))

# objective function value
getObjValGLPK(lp)

# value of variables in optimal solution
Solution.Variables <- as.data.frame(getColsPrimGLPK(lp))



sum(Solution.Variables)


rownames(Solution.Variables) <- colnames(Data)
write.csv(Solution.Variables, "Solution_Variables.csv", row.names = TRUE)

# status of each variable in optimal solution 1 = basic variable
getColsStatGLPK(lp)
write.csv(getColsStatGLPK(lp), "Basic_Variables.csv", row.names = TRUE)

# get dual values for each row/constraint 
getRowsDualGLPK(lp)
write.csv(getRowsDualGLPK(lp), "Dual_Values_Constraint_Shadows.csv", row.names = TRUE)

# get dual values for each variables // reduced costs
getColsDualGLPK(lp)
write.csv(getColsDualGLPK(lp), "Variable_Reduced_Costs.csv", row.names = TRUE)


mipObjValGLPK(lp) 
printMIPGLPK(lp,fname="mipReport.txt")