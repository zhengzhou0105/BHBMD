# Settings--------
seed <- 47405
iter <- 1E4
warmup <- 0.5
Nchain <- 3
thin <- 1

source("BHBMD model utilities.R")

# Input Raw Data-----

df_raw_bladder_meta <- read.csv("bladder cancer ~ ADD all 0520.csv")
df_raw_bladder_cc <- filter(df_raw_bladder_meta,
                            Design == "Case-Control")
df_raw_bladder_Shao <- read.csv("Bladder_Cancer_Input.csv")
df_raw_bladder_Allen <- read.csv("Allen 2020 iAs bladder cancer data.csv")
df_raw_example <- readxl::read_xlsx("Validation Hill 0524.xlsx",
                                    sheet = "test1")

#----- Bladder cancer--------
# Before unifying format, calculate SD of RR 
df_raw_bladder_meta <- df_raw_bladder_meta %>% mutate(
  ysdL = log(RR.u/RR.l)/(2*1.96)
) %>% mutate(
  ysdL = ifelse(is.nan(ysdL),0,ysdL)
)
df_bladder_meta <- RRDataUnify(df_raw = df_raw_bladder_meta,
                               Index = "?..Index",Study = "Study",N = "N",
                               dose = "ADD.Total_ug.kg.D",
                               RR = "RR",ysdL = "ysdL",
                               RRlower = "RR.l",RRupper = "RR.u")

# #* add single study normalization
# # standardized dose if as single study
# tempdose_minmax <- vector(mode = "numeric",length = 0)
# tempdose_meansd <- tempdose_minmax
# for(g in 1:10){
#   tempdf <- bladder_meta_all %>% filter(Index == g)
#   tempADD <- tempdf$ADD 
#   tempdose_minmax <- c(tempdose_minmax, dnormalize_minmax(tempADD))
#   tempdose_meansd <- c(tempdose_meansd, dnormalize_meanSD(tempADD))
# }
# 
# bladder_meta_all$dose_minmax_single <- tempdose_minmax
# bladder_meta_all$dose_meansd_single <- tempdose_meansd
# 
# # compare the ratio between standardized dose between meta and single study
# bladder_meta_all %>% mutate(
#   ratio = dose_minmax_single / dose_minmax
# ) %>% select(Index, dose_minmax, dose_minmax_single, ratio)
# 
# # Output the ADD, minmax standardized dose of each study
# # comparison between single study and meta standardization
# ListSingleVSMeta <- vector("list",10)        # create a list with given length to save time
# # add names for elements in the list
# names(ListSingleVSMeta) <- c(paste0("Study #",1:10," ",unique(bladder_meta_all$Study))[1:8],
#                              "Study #9 Sawada (2013)", "Study #10 Bates (1995)")
# for(g in 1:10){
#   # extract and store df by study index
#   tempdf <- bladder_meta_all %>% filter(Index == g)
#   # create a dataframe for each study, storing ADD, two doses and the ratio
#   ListSingleVSMeta[[g]] <- data.frame(
#     ADD = tempdf$ADD,
#     Dose_minmax_single = tempdf$dose_minmax_single,
#     Dose_minmax_meta = tempdf$dose_minmax,
#     Ratio = tempdf$dose_minmax_single / tempdf$dose_minmax
#   )
# }
# # Output each study to individual worksheets
# filename <- paste("Dose_minmax Single vs Meta.xlsx")  # store the filename
# wb <- createWorkbook()      # create empty workbook
# sheets <- lapply(names(ListSingleVSMeta),createSheet,wb = wb)  # create and define worksheets name
# voids <- Map(addDataFrame,ListSingleVSMeta,sheets)  # mapping elements to sheets pairwise 
# # saveWorkbook(wb, filename)   # output the workbook

#*Case control only-----
#*# Before unifying format, calculate SD of RR
df_raw_bladder_cc <- df_raw_bladder_cc %>% mutate(
  ysdL = log(RR.u/RR.l)/(2*1.96)
) %>% mutate(
  ysdL = ifelse(is.nan(ysdL),0,ysdL)
)
df_bladder_meta_cc <- RRDataUnify(df_raw = df_raw_bladder_cc,
                               Index = "?..Index",Study = "Study",N = "N",
                               dose = "ADD.Total_ug.kg.D",
                               RR = "RR",ysdL = "ysdL",
                             RRlower = "RR.l",RRupper = "RR.u"
                             )
# note the study index needs to be formatted
df_bladder_meta_cc$Index <- c(rep(1,6),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3),rep(7,4))


#*Allen 2010 data--------
df_raw_bladder_Allen <- df_raw_bladder_Allen %>% mutate(
  ysdL = log(RR.u / RR.l) / (2*1.96)
) %>% mutate(
  ysdL = ifelse(is.na(ysdL),0,ysdL)
)
df_bladder_Allen <- RRDataUnify(df_raw = df_raw_bladder_Allen,
                                Index = "Index",Study = "Study",N = "N",
                                dose = "ADD", RR = "RR",ysdL = "ysdL",
                                RRlower = "RR.l",RRupper = "RR.u"
                                )
# lower and upper bounds of RR at reference were NAs
df_bladder_Allen <- mutate(df_bladder_Allen,
    RR.l = ifelse(is.na(RR.l),0,RR.l),
    RR.u = ifelse(is.na(RR.u),0,RR.u)
                           )


# *bladder_Shao----------
OrderShao <- c("Baris 2016","Bates 2004","Chen 2010","Karagas 2004",
               "Kurttio 1999","Meliker 2010","Mostafa 2015","Stainmaus 2013")
df_raw_bladder_Shao <- df_raw_bladder_Shao %>% group_by(Study_ID) %>% 
  mutate(Study = OrderShao[Study_ID])
df_bladder_Shao <- RRDataUnify(df_raw = df_raw_bladder_Shao,
                               Index = "Study_ID",N="N",Study = "Study",
                               RRlog = "Mean",ysdL = "SD",dose = "Dose")

# #Example data continuous
# df_example <- RRDataUnify(df_raw = df_raw_example,
#                           Index = NULL,N = "N",Study = NULL,
#                           dose = "Dose",RR = "Response",ysd = "SD")
# 


# *Use only first 3 doses  for all data--------
df_bladder_Allen_nonrag <- df_bladder_Allen %>% group_by(Index) %>% 
  mutate(group = row_number()) %>% filter(group <= 3)

df_bladder_meta_nonrag <- df_bladder_meta %>% group_by(Index) %>% 
  mutate(group = row_number()) %>% filter(group <= 3) 

df_bladder_meta_cc_nonrag <- df_bladder_meta_cc %>% group_by(Index) %>% 
  mutate(group = row_number()) %>% filter(group <= 3)

df_bladder_Shao_nonrag <- df_bladder_Shao %>% group_by(Index) %>% 
  mutate(group = row_number()) %>% filter(group <= 3)

# Output------------
save(
  df_bladder_meta,df_bladder_meta_cc,df_bladder_Allen,df_bladder_Shao,
  df_bladder_Allen_nonrag,df_bladder_meta_nonrag,df_bladder_meta_cc_nonrag,
  df_bladder_Shao_nonrag,
  file = "Hill meta base data.Rdata"
)
