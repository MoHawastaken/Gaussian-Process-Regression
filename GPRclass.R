#needs package R6

GPR <- R6::R6Class("GPR",
                   private = list(
                     .X = NA,
                     .k = NA,
                     .y = NA,
                     .noise = NA
                   ),
                   public = list(
                     initialize = function(X, k, y, noise){
                       .X <-  X
                       .k <-  k
                       .y <- y
                       .noise <- noise
                     }
                   ),
                   active = list(
                     X = function(value){
                       if(missing(value)){
                         private$.X
                       } else{
                         stop("`$X` is read only", call. = FALSE)
                       }
                     },
                     k = function(value){
                       if(missing(value)){
                         private$.k
                       } else{
                         stop("`$k` is read only", call. = FALSE)
                       }
                     },
                     y = function(value){
                       if(missing(value)){
                         private$.y
                       } else{
                         stop("`$y` is read only", call. = FALSE)
                       }
                     },
                     noise = function(value){
                       if(missing(value)){
                         private$.noise
                       } else{
                         stop("`$noise` is read only", call. = FALSE)
                       }
                     }
                   )
  
)

GPR.linear <- 