read.data<-function(dir) {
   subjs<-list.files(dir)
   
   for (i in 1:length(subjs)) {
      if(i==1) {
         df<-read.csv(file.path(dir,subjs[i]))
      }
      temp<-read.csv(file.path(dir,subjs[i]))
      df<-rbind(df,temp)
      df<-as.data.frame(df)
   }
   return(df)
}