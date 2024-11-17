
    options(error=stop)

    logit <- function(prob){ log(prob) - log(1-prob)}
    
    expit <- function(logodds){ 1/(1+exp(-logodds))}
    
    
    getlogop <- function(p0,p1){
        log(p0) + log(p1) - log(1-p0) - log(1-p1)
    }
    
    getlogor <- function(p0,p1){
        -log(p0) + log(p1) + log(1-p0) - log(1-p1)
    }
    
    getlogrr <- function(p0,p1){
        log(p1) - log(p0)
    }
    
    #### Function for checking if two things are equal 
    #### within numerical precision
    same <- function(x,y){isTRUE(all.equal(x,y))}
    
    
    ### Function for finding probabilities from logrr and logop
    ###
    getProbScalarRDiff <- function(atanhrd,logop,weight = 1){
        
        rd = tanh(atanhrd) * weight
        
        if(logop>350){
            if (atanhrd < 0){
                p0 <- 1
                p1 <- p0+rd
            }else{
                p1<-1
                p0 <- p1-rd}   		
        }
        else{ ## not on boundary
            if(same(logop,0)){## logop = 0; solving linear equations 
                p0<-0.5*(1-rd)}
            else{
                p0 <- (-(exp(logop)*(rd-2)-rd)-
                           sqrt((exp(logop)*(rd-2)-rd)^2+4*exp(logop)*(1-rd)*(1-exp(logop))))/(2*(exp(logop)-1))}
            p1 <- p0 + rd}
        return(c(p0,p1))
    }
    
  