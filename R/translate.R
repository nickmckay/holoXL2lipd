
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

#
# library(lipdR)
# library(lipdverseR)
# library(readxl)
# library(magrittr)
# library(dplyr)
#
# path <- "~/Dropbox/CLIMATE12k excel formatted/Chironomid/"
# outPath <- "~/Dropbox/CLIMATE12k excel formatted/ChironomidLipd/"
#
# dir.create(outPath)
#
# fname <- list.files(path,pattern = "*.xlsx")
# good = !stringr::str_detect(fname, "[~$]")
# fname <- fname[good]
#
# for(i in fname){
#   L <- climate12k_excel_to_lipd_converter(path = path,fname = i)
#   writeLipd(L,path = outPath)
#}

climate12k_excel_to_lipd_converter <- function(path,fname){

  print(paste("Converting",fname))
  convo <- getConverter()
  xl <- readxl::read_xlsx(file.path(path,fname))

  #clean up special characters
  rosetta <- lipdverseR:::rosettaStone()
  xl <- purrr::map_df(xl,lipdverseR:::replaceSpecialCharacters,rosetta)



  fa <- stringr::str_extract(fname,pattern = "^[^_]+(?=_)")
  py <- stringr::str_extract(fname,pattern = "(?<=_)[^_]+(?=_)")



  # read metadata -----------------------------------------------------------

  xlm <- xl[,1:2]
  names(xlm) <- c("key","value")

  xlm <- dplyr::filter(xlm,!is.na(value))

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
  sn <- lipdverseR:::replaceSpecialCharacters( nts$geo_siteName,rosetta)
  nts$dataSetName <- stringr::str_c(sn,".",fa,".",py) %>%
    stringr::str_replace_all(" ","")


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

  toCheck <- c("Top Depth of Date (cm)",
               "Bottom Depth of Date (cm)",
               "Date Type")
  w <- 1
  while(length(c2)==0){
    c2 <- which(names(xlp) == toCheck[w])
    w = w+1
    if(w>100){
      c2 <- NULL
      break
    }
  }

  if(length(c2)==1){
    hasChron = TRUE
  }else{
    hasChron = FALSE
  }

  if(hasChron){
    #isolate chron
    ct <- xlp[,c2:ncol(xlp)]

    ct <- ct[-which(rowSums(!is.na(ct)) == 0), ]

    #isolate paleo
    pt <- xlp[,-c((c2-1):ncol(xlp))]
  }else{
    pt <- xlp[,-ncol(xlp)]
  }




  #make paleo TS
  ts <- vector(mode = "list",length = ncol(pt))

  for(i in 1:length(ts)){
    #copy in metadata
    ts[[i]] <- nts

    ts[[i]]$paleoData_values <- as.matrix(pt[,i])

    #parse name
    ts[[i]]$paleoData_variableName <- stringr::str_replace_all(string = names(pt)[i],pattern =" ","") %>%
      stringr::str_extract(pattern = "^[^(]+(?=())") %>%
      stringr::str_replace_all(pattern ="[^A-Za-z0-9]","")


    #try to get units
    ts[[i]]$paleoData_units <- stringr::str_replace_all(string = names(pt)[i],pattern =" ","") %>%
      stringr::str_extract(pattern = "(?<=[(])[^_]+(?=[)])")

    #generate TSid
    ts[[i]]$paleoData_TSid <- lipdR::createTSid()


    #look for special metadata
    #Depth
    if(ts[[i]]$paleoData_variableName == "Depth"){
      if(!is.na(xl[8,6])){
        ts[[i]]$paleoData_sampleThickness <- as.numeric(xl[8,6])
      }
      if(!is.na(xl[9,6])){
        ts[[i]]$paleoData_depthReference <- as.character(xl[9,6])
      }
      if(!is.na(xl[10,6])){
        ts[[i]]$paleoData_notes <- as.character(xl[10,6])
      }
    }

    #TempRecon1
    if(ts[[i]]$paleoData_variableName == "TemperatureReconstruction1"){
      ts[[i]]$interpretation1_variable <- "T"
      ts[[i]]$interpretation1_direction <- "positive"
      ts[[i]]$interpretation1_scope <- "climate"
      if(ts[[i]]$timeseriesType == "Uncalibrated"){
        ts[[i]]$paleoData_units <- NA
      }else{
        ts[[i]]$paleoData_units <- "degC"
      }

      if(!is.na(xl[6,12])){
        ts[[i]]$interpretation1_seasonality <- as.character(xl[6,12])
      }
      if(!is.na(xl[7,12])){
        ts[[i]]$calibration_uncertaintyType <- as.character(xl[7,12])
      }
      if(!is.na(xl[8,12])){
        ts[[i]]$calibration_method <- as.character(xl[8,12])
      }
      if(!is.na(xl[9,12])){
        ts[[i]]$paleoData_modernTemperature <- as.character(xl[9,12])
      }
      if(!is.na(xl[10,12])){
        ts[[i]]$paleoData_notes <- as.character(xl[10,12])
      }
    }

    #TempRecon2
    if(ts[[i]]$paleoData_variableName == "TemperatureReconstruction2"){
      ts[[i]]$interpretation1_variable <- "T"
      ts[[i]]$interpretation1_direction <- "positive"
      ts[[i]]$interpretation1_scope <- "climate"



      if(!is.na(xl[6,18])){
        ts[[i]]$interpretation1_seasonality <- as.character(xl[6,18])
      }
      if(!is.na(xl[7,18])){
        ts[[i]]$calibration_uncertaintyType <- as.character(xl[7,18])
      }
      if(!is.na(xl[8,18])){
        ts[[i]]$calibration_method <- as.character(xl[8,18])
      }
      if(!is.na(xl[9,18])){
        ts[[i]]$paleoData_modernTemperature <- as.character(xl[9,18])
      }
      if(!is.na(xl[10,18])){
        ts[[i]]$paleoData_notes <- as.character(xl[10,18])
      }
    }

    #TempRecon3
    if(ts[[i]]$paleoData_variableName == "TemperatureReconstruction3"){
      ts[[i]]$interpretation1_variable <- "T"
      ts[[i]]$interpretation1_direction <- "positive"
      ts[[i]]$interpretation1_scope <- "climate"
      ts[[i]]$paleoData_useInGlobalTemperatureAnalysis <- "?"

      if(!is.na(xl[6,24])){
        ts[[i]]$interpretation1_seasonality <- as.character(xl[6,24])
      }
      if(!is.na(xl[7,24])){
        ts[[i]]$calibration_uncertaintyType <- as.character(xl[7,24])
      }
      if(!is.na(xl[8,24])){
        ts[[i]]$calibration_method <- as.character(xl[8,24])
      }
      if(!is.na(xl[9,24])){
        ts[[i]]$paleoData_modernTemperature <- as.character(xl[9,24])
      }
      if(!is.na(xl[10,24])){
        ts[[i]]$paleoData_notes <- as.character(xl[10,24])
      }
    }
  }

  L <- collapseTs(ts,force = TRUE)

  #assign in some metadata
  L$lipdVersion <-1.3
  L$createdBy <- "holoXL2lipd"

  if(hasChron){
    #make up a chronData
    L$chronData <- L$paleoData


    cts <- vector(mode = "list",length = ncol(ct))


    for(i in 1:length(cts)){
      #copy in metadata
      cts[[i]] <- nts
      cts[[i]]$paleoData_values <- as.matrix(ct[,i])

      #parse name
      cts[[i]]$paleoData_variableName <- stringr::str_replace_all(string = names(ct)[i],pattern =" ","") %>%
        stringr::str_extract(pattern = "^[^(]+(?=())") %>%
        stringr::str_replace_all(pattern ="[^A-Za-z0-9]","")


      #try to get units
      cts[[i]]$paleoData_units <- stringr::str_replace_all(string = names(ct)[i],pattern =" ","") %>%
        stringr::str_extract(pattern = "(?<=[(])[^_]+(?=[)])")

      #generate TSid
      cts[[i]]$paleoData_TSid <- lipdR::createTSid()


      #look for special metadata
      #Depth
      if(cts[[i]]$paleoData_variableName == "DateBP"){
        if(!is.na(xl[8,30])){
          cts[[i]]$paleoData_ageModelSource <- as.character(xl[8,30])
        }
        if(!is.na(xl[9,30])){
          cts[[i]]$paleoData_depthReference <- as.character(xl[9,30])
        }
        if(!is.na(xl[10,30])){
          cts[[i]]$paleoData_notes <- as.character(xl[10,30])
        }
      }

    }
    C <- collapseTs(cts,force = TRUE)

    L$chronData <- C$paleoData
  }

  return(L)
}



