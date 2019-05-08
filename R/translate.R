
getConverter <- function(googId = "1CNafsRfF8y79RQ1KzaBMBZUahMggsYBe8_SbPt5fQHU"){

  fname <- file.path(tempdir(),stringr::str_c(googId,".csv"))

#check time
if(is.na(lubridate::now()-file.mtime(fname)) | lubridate::now()-file.mtime(fname) > 30){

#download name conversion
convo <- googledrive::as_id(googId) %>%
  googledrive::drive_get() %>%
  googledrive::drive_download(path = fname,overwrite = T) %>%
  dplyr::select(local_path) %>%
  as.character() %>%
  readr::read_csv()

}else{

  convo <- readr::read_csv(fname)

}

return(convo)
}



library(readxl)
library(magrittr)
library(dplyr)

path <- "~/Dropbox/CLIMATE12k excel formatted/Chironomid/"

fname <- "Antonsson_2006_chironomids_Fennoscandia_Gilltjärnen_checkedSE.xlsx"

xl <- read_xlsx(file.path(path,fname))

fa <- stringr::str_extract(fname,pattern = "^[^_]+(?=_)")
py <- stringr::str_extract(fname,pattern = "(?<=_)[^_]+(?=_)")



# read metadata -----------------------------------------------------------

xlm <- xl[,1:2]
names(xlm) <- c("key","value")

xlm <- filter(xlm,!is.na(value))

#get converter
convo <- getConverter()

nts <- list()

#convert and store in TS
for(i in 1:nrow(xlm)){
  ind <- which(xlm$key[i] == convo$climate12kName)
  if(length(ind)>1){
    stop("multiple matches with key")
  }else if(length(ind)==0){
    warning(stringr::str_c("no conversion match for key: ", xlm$key[i],", skipping..."))
  }else{
    #check type
    varType <- convo$type[ind]
    if(varType == "character"){
      varFun <- as.character
    }else if(varType == "numeric"){
      varFun <- as.numeric
    }else if(varType == "boolean"){
      varFun <- as.logical
    }else{
      stop("variable type not recognized")
    }
    nts[[convo$tsName[ind]]] <- varFun(xlm$value[i])
  }
}

#create dsn
sn <- stringr::str_replace_all(string = nts$geo_siteName,pattern = "[^A-Za-z]","")
nts$dataSetName <- stringr::str_c(fa,".",sn,".",py)


# get timeseries ----------------------------------------------------------

#find where to start
rc <- which(xl == "Original Sample ID",arr.ind = TRUE)

#chop out the data
xlt <- xl[rc[1]:nrow(xl),rc[2]:ncol(xl)]

#remove columns that are all NAs
not_all_na <- function(x) {!all(is.na(x[-1]))}

xlp <- select_if(xlt,not_all_na)
names(xlp) <- xlp[1,]
xlp <- xlp[-1,]

#separate into paleo and chron
c2 <- which(names(xlp) == "Original Date ID")

#isolate chron
ct <- xlp[,c2:ncol(xlp)]

#isolate paleo
pt <- xlp[,-c((c2-1):ncol(xlp))]


#make paleo TS
ts <- vector(mode = "list",length = ncol(pt))

for(i in 1:length(ts)){
  #copy in metadata
  ts[[i]] <- nts

  ts[[i]]$paleoData_values <- pt[,1]

  #parse name
  ts[[i]]$paleoData_variableName <- stringr::str_replace_all(string = names(pt)[i],pattern =" ","") %>%
    stringr::str_extract(pattern = "^[^(]+(?=())")

  #try to get units
  ts[[i]]$paleoData_units <- stringr::str_replace_all(string = names(pt)[i],pattern =" ","") %>%
    stringr::str_extract(pattern = "(?<=[(])[^_]+(?=[)])")

  #generate TSid
  ts[[i]]$paleoData_TSid <- lipdR::createTSid()


}



