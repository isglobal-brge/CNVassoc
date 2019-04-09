## replace blanks o multiple of 4 spaces.

#### R

files <- list.files("./R", full=TRUE)

for (ff in files){
  
  print(ff)

  code <- scan(ff, what="character", sep="\n")
  ii <- grep("^ ", code)
  ii
  code2 <- code
  
  for (i in ii){
    linia <- code[i]
    linia.split <- strsplit(linia,"")[[1]]
    counts <- 0
    k <- 1
    while(k<=length(linia.split) & linia.split[k]==" "){
      k <- k+1
    }
    blanks <- (k%/%4)*4+4*((k%%4)>0)
    code2[i] <- sub("^[ ]+",paste(rep(" ",blanks), collapse=""), linia)
  }
  
  
  write(code2, file=ff)

}


#### man

files <- list.files("./man", full=TRUE)


for (ff in files){
  
  print(ff)
  
  code <- scan(ff, what="character", sep="\n")
  ii <- grep("^ ", code)
  ii
  code2 <- code
  
  for (i in ii){
    linia <- code[i]
    linia.split <- strsplit(linia,"")[[1]]
    counts <- 0
    k <- 1
    while(k<=length(linia.split) & linia.split[k]==" "){
      k <- k+1
    }
    blanks <- (k%/%4)*4+4*((k%%4)>0)
    code2[i] <- sub("^[ ]+",paste(rep(" ",blanks), collapse=""), linia)
  }
  
  
  write(code2, file=ff)
  
}


#### vignette

ff <- "./vignettes/CNVassoc_vignette.Rmd"

code <- scan(ff, what="character", sep="\n")
ii <- grep("^ ", code)
ii
code2 <- code
  
for (i in ii){
  linia <- code[i]
  linia.split <- strsplit(linia,"")[[1]]
  counts <- 0
  k <- 1
  while(k<=length(linia.split) & linia.split[k]==" "){
    k <- k+1
  }
  blanks <- (k%/%4)*4+4*((k%%4)>0)
  code2[i] <- sub("^[ ]+",paste(rep(" ",blanks), collapse=""), linia)
}
  
write(code2, file=ff)


