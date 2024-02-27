# 
# Copyright © 2017, Université catholique de Louvain
# All rights reserved.
# 
# Copyright © 2017 Forschungszentrum Jülich GmbH
# All rights reserved.
# 
# Developers: Guillaume Lobet
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted under the GNU General Public License v3 and provided that the following conditions are met:
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#Disclaimer
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# You should have received the GNU GENERAL PUBLIC LICENSE v3 with this file in license.txt but can also be found at http://www.gnu.org/licenses/gpl-3.0.en.html
# NOTE: The GPL.v3 license requires that all derivative work is distributed under the same license. That means that if you use this source code in any other program, you can only distribute that program with the full source code included and licensed under a GPL license.




# FUNCTION TO READ THE PPARAM FILES FROM CROOTBOX AND STORE IT INTO A DATAFRAME
read_rparam <- function(path){

  fileName <- path
  param <- read_file(fileName)
  
  
  param <- strsplit(param, "#")
  dataset_init <- NULL
  for(k in c(2:length(param[[1]]))){
    spl <- strsplit(param[[1]][k], "\n")
    type <- ""
    name <- ""
    for(i in c(1:length(spl[[1]]))){
      temp <- spl[[1]][i]
      pos <- regexpr("//", temp)
      if(pos != -1) temp <- substr(temp, 0, pos-1)
      if(nchar(temp) > 0){
        temp <- strsplit(temp, "\t")
        temp2 <- data.frame("type" = character(0), "name" = character(0), 
                            "param" = character(0), "val1" = numeric(0),
                            #Addition of val4
                            "val2" = numeric(0), "val3" = numeric(0), "val4" = numeric(0), stringsAsFactors = F)
        
        if(temp[[1]][1] == "type"){ type <- temp[[1]][2]
        } else if(temp[[1]][1] == "name"){ name <- temp[[1]][2]
        } else if(grepl("Param", temp[[1]][1])){
        } else if(temp[[1]][1] == "tropism") {
          temp2[[1,3]] <- "n_tropism"
          temp2$val1 <- temp[[1]][3]
          temp2$type <- type
          temp2$name <- name
          dataset_init <- rbind(dataset_init, temp2)
          temp2$param <- "sigma_tropism"
          temp2$val1 <- temp[[1]][4]
          temp2$type <- type
          temp2$name <- name
          dataset_init <- rbind(dataset_init, temp2)  
          temp2$param <- "tropism"
          temp2$val1 <- temp[[1]][2]
          temp2$type <- type
          temp2$name <- name
          dataset_init <- rbind(dataset_init, temp2)  
        } else {
          for(j in c(1:5)){
            temp2[[1,j+2]] <- temp[[1]][j]
            temp2$type <- type
            temp2$name <- name
          }
          dataset_init <- rbind(dataset_init, temp2)
        }
      }
    }
  } 
  
  return(dataset_init)
}


# FUNCTION TO READ THE PPARAM FILES FROM CROOTBOX AND STORE IT INTO A DATAFRAME
read_pparam <- function(path){
  ## READ THE PARAMETER FILE AND STORE THE DATA IN A DATAFRAME
  data <- read_file(path)
  # READ THE PARAMETER FILE AND STORE THE DATA IN A DATAFRAME
  plant_init <- NULL
  spl <- strsplit(data, "\n")
  for(i in c(1:length(spl[[1]]))){
    temp <- spl[[1]][i]
    if(nchar(temp) > 0){
      temp <- strsplit(temp, "\t")
      temp2 <- data.frame( "param" = character(0), "val1" = numeric(0), stringsAsFactors = F)
      for(j in c(1:2)){
        temp2[[1,j]] <- temp[[1]][j]
      }
      plant_init <- rbind(plant_init, temp2)
    }
  }      
  
  colnames(plant_init) <- c("param", "val1")  
  return(plant_init)
}




# FUNCTION TO WRITE THE RPARAM FILES FROM CROOTBOX AND STORE IT INTO A DATAFRAME

write_rparam <- function(dataset, files){
  
  types <- unique(dataset$type)
  text <- NULL
  for(t in types){
    if(is.null(text)){text <- "# Parameter set for type"
    }else{
      text <- paste(text, "# Parameter set for type", sep="\n")
    }
    
    temp <- dataset[dataset$type == t,]
    
    str <- paste("type", temp$type[1], sep="\t")
    text <- paste(text, str, sep="\n")
    
    str <- paste("name", temp$name[1], sep="\t")
    text <- paste(text, str, sep="\n")
    
    for(i in c(1:nrow(temp))){
      if(temp[i, 3] == "n_tropism"){
        str <- paste("tropism", temp[i+2, 4], temp[i, 4], temp[i+1, 4], sep="\t")
        text <- paste(text, str, sep="\n")
      }else if(temp[i, 3] == "sigma_tropism" | temp[i, 3] == "tropism"){
      }else if(temp[i, 3] == "dx"){
        str <- paste(temp[i, 3], temp[i, 4], sep="\t")
        text <- paste(text, str, sep="\n")
      }else{
        str <- paste(temp[i, 3], temp[i, 4], temp[i, 5], temp[i, 6], temp[i, 7], sep="\t")
        text <- paste(text, str, sep="\n")
      }
    }
    
  }
  text <- gsub("\tNA", "", text)
  for(f in files){
    cat(text, file=f)
  }

}

# FUNCTION TO WRITE THE PPARAM FILES FROM CROOTBOX AND STORE IT INTO A DATAFRAME

write_pparam <- function(plant, files){
  text <- NULL
  for(i in c(1:nrow(plant))){
    str <- paste(plant[i, 1], plant[i, 2], sep="\t")
    text <- paste(text, str, sep="\n")
  }
  
  text <- gsub("\tNA", "", text)

  for(f in files){
    cat(text, file=f)
  }
}


# Needs the node table from SmartRoot
# rs = DataFrame with the node data
# order = the order of root for which the sigma needs to be computed
getSigma <- function(rs, order = 1){
  
  require(plyr)
  
  rs <- rs %>% 
    filter(root_order == order) %>% 
    ddply(.(#image, 
      root), summarise, 
          diff_o = abs(diff(theta)), 
          diff_l = diff(distance_from_base)) %>% 
    filter(diff_o > 0) %>%
    ddply(.(image), summarise, theta = mean(diff_o / diff_l))
  
  rs$sigma <- exp(1.484 - 0.00798 * rs$theta)
  
  return(rs)
  
}

getSigma2_0 <- function(rs, name = "Lat"){
  
  # require(plyr)
  
  rs <- rs %>% 
    filter(root_name == name) %>% 
    group_by(image, root) %>% 
    mutate(diff_o = c(0,abs(diff(theta))), 
           diff_l = c(0,diff(distance_from_base))) %>% 
    filter(diff_o > 0) %>%
    ungroup %>% 
    group_by(image) %>% 
    summarise(theta = mean(diff_o / diff_l)) %>% 
    mutate(sigma = exp(1.484 - 0.00798 * theta))
  # 
  # rs3 <- rs %>% 
  #   filter(root_name == name) %>% 
  #   ddply(.(image, root), summarise, diff_o = abs(diff(theta)), 
  #          diff_l = diff(distance_from_base)) %>% 
  #   filter(diff_o > 0) %>%
  #   ddply(.(image), summarise, theta = mean(diff_o / diff_l))
  
  # rs$sigma <- exp(1.484 - 0.00798 * rs$theta)
  
  return(rs)
  
}
  
  

