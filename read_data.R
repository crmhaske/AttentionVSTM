## =============================================================================
## Load all subject data as one data frame v.1
## By Christie Haskell, 27-07-2014 (v.1)
## Source:
## =============================================================================

read.data<-function(dir) {
   subjs<-list.files(dir)
   
   if (length(subjs) == 1) {
     df<-read.csv(file.path(dir,subjs))
   }
   else {
    for (i in 1:length(subjs)) {
        if(i==1) {
           df<-read.csv(file.path(dir,subjs[i]))
        }
        else {
          temp<-read.csv(file.path(dir,subjs[i]))
          df<-rbind(df,temp)
          df<-as.data.frame(df)
        }
     }
    return(df)
   }
}