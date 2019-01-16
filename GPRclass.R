#needs package R6

GPR <- R6::R6Class("GPR",
                   private = list(
                     .X = NA,
                     .k = NA,
                     .y = NA,
                     .noise = NA
                   ),
                   public = list(
                     initialize = function(X, y, k, noise){
                       private$.X <-  X
                       private$.k <-  k
                       private$.y <- y
                       private$.noise <- noise
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

GPR.constant <- R6::R6Class("GPR.constant",
                          inherit = GPR,
                          public = list(
                            initialize = function(X, y, c, noise){
                              k <- function(x, y) c
                              super$initialize(X, y, k, noise)
                            }
                          )
)                          