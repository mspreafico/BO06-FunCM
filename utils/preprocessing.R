library(data.table)

impute_cyc6<-function(data,varname){
  pts<-unique(data[is.na(get(varname)) & cycno==6]$patid)
  for(i in 1:length(pts)){
    value<-unlist(data[patid==pts[i] & cycno==1, ..varname])
    sdev<-sd(unlist(data[cycno==1, ..varname]),na.rm=T)
    patients<-unique(data[cycno==1 & get(varname)>=value-sdev & get(varname)<=value+sdev]$patid)
    med<-median(unlist(data[cycno==6 & patid %in% patients, ..varname]),na.rm=T)
    data[patid==pts[i] & cycno==6, c(varname):= med]
  }
  return(data)
}

impute_ntest<-function(data,varname,n.test){
  pts<-unique(data[is.na(get(varname)) & ntest==n.test]$patid)
  for(i in 1:length(pts)){
    value<-unlist(data[patid==pts[i] & ntest==1, ..varname])
    sdev<-sd(unlist(data[ntest==1, ..varname]),na.rm=T)
    patients<-unique(data[ntest==1 & get(varname)>=value-sdev & get(varname)<=value+sdev]$patid)
    med<-median(unlist(data[ntest==n.test & patid %in% patients, ..varname]),na.rm=T)
    data[patid==pts[i] & ntest==n.test, c(varname):= med]
  }
  return(data)
}
