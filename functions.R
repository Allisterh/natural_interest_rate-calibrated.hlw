
shiftQuarter <- function(original.start,shift){

    if (shift > 0) {
        new.start = c(0,0)
        sum = original.start[2] + shift
    
        if (sum <= 4) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] + ceiling(sum/4) - 1
        }

        if (sum %% 4 > 0) {
            new.start[2] = sum %% 4
        }
        else {
            new.start[2] = sum %% 4 + 4
        }
    }

    else {
        new.start = c(0,0)
        diff = original.start[2] - abs(shift)
    
        # Get the year value
        if (diff > 0) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] - (1 + floor(abs(diff)/4))
        }

        if (diff %% 4 > 0) {
            new.start[2] = diff %% 4
        }
        else {
            new.start[2] = diff %% 4 + 4
        }
    }
        
return(new.start)}


shiftMonth <- function(original.start,shift){

    if (shift > 0) {
        new.start = c(0,0)
        sum = original.start[2] + shift
    
        if (sum <= 12) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] + ceiling(sum/12) - 1
        }

        if (sum %% 12 > 0) {
            new.start[2] = sum %% 12
        }
        else {
            new.start[2] = sum %% 12 + 12
        }
    }

    else {
        new.start = c(0,0)
        diff = original.start[2] - abs(shift)
    
        if (diff > 0) {
            new.start[1] = original.start[1]
        }
        else {
            new.start[1] = original.start[1] - (1 + floor(abs(diff)/12))
        }

        if (diff %% 12 > 0) {
            new.start[2] = diff %% 12
        }
        else {
            new.start[2] = diff %% 12 + 12
        }
    }
        
return(new.start)}


getFRED <- function(url, freq = "Quarterly") {
    
    txt.file.name <- paste0("rawData/",substr(url, regexpr('[a-zA-z0-9]*.txt',url),1000))
    if (!file.exists(txt.file.name)){
        # Download the data from FRED
        #download.file(url, destfile = 'FREDtemp.txt', method = "wget")
        system(paste0('wget --no-check-certificate "', url, '"'))
        system(paste('mv',substr(url, regexpr('[a-zA-z0-9]*.txt',url),1000),txt.file.name))
    }
    FREDraw <- readLines(txt.file.name) 

    freq.FRED <- gsub(' ', '',substr(FREDraw[which(regexpr('Frequency', FREDraw)==1)],
                                     (nchar('Frequency')+2),100))    

    datastart = which(gsub(' ', '',FREDraw)=='DATEVALUE') - 2

    data <- read.table(txt.file.name, skip = datastart, header = TRUE)

    first.year  <- as.numeric(format(as.Date(data$DATE[1]),'%Y'))
    first.month <- as.numeric(format(as.Date(data$DATE[1]),'%m'))
    
    if (freq.FRED == 'Quarterly'){
        first.q  <- (first.month-1)/3 + 1
        data.tis <- tis(data$VALUE, start = c(first.year, first.q), tif = 'quarterly')
    } else if (freq.FRED == 'Monthly') {
        data.tis <- tis(data$VALUE, start = c(first.year, first.month), tif = 'monthly')
    }

    if (freq.FRED == 'Monthly' & freq == 'Quarterly') {
        data.tis <- convert(data.tis, tif = 'quarterly', method = 'constant', observed. = 'averaged')
    }

    return(data.tis)
} 


splice <- function(s1, s2, splice.date, freq) {
  
    t <- splice.date
    if (freq == "quarterly" | freq == "Quarterly") {
        t.minus.1 <- shiftQuarter(t,-1)
    }
    else if (freq == "monthly" | freq == "Monthly") {
        t.minus.1 <- shiftMonth(t,-1)
    }
    else { stop("You must enter 'quarterly' or 'monthly' for freq.") }
    ratio <- as.numeric(window(s2,start = t, end = t)/
                        window(s1,start = t, end = t))

    return(mergeSeries(ratio*window(s1,end = t.minus.1),window(s2, start = t)))
}

gradient <- function(f, x, delta = x * 0 + 1.0e-5) {

    g <- x * 0
    for (i in 1:length(x)) {
        x1 <- x
        x1[i] <- x1[i] + delta[i]
        f1 <- f(x1)
        x2 <- x
        x2[i] <- x2[i] - delta[i]
        f2 <- f(x2)
        g[i] <- (f1 - f2) / delta[i] / 2
    }
    return(g)
}
