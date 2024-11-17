
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
    
    
    # MyAttach = function(list){  # attach a list to the global environment
    #     for(i in 1:length(names(list))){
    #         assign(names(list)[i], list[[i]], envir=.GlobalEnv)
    #     }
    # }
    # 
    # 
    # Optimize = function(objective, startpars, thres = 1e-6, max.step = 1000){
    #     
    #     Diff = function(x,y) abs(x-y)/(x+thres)
    #     step = 0;    value.old = diff = thres + 1
    #     while(diff > thres & value.old > thres & step < max.step){
    #         step = step + 1
    #         opt = optim(startpars,objective,control=list(maxit=max.step))
    #         diff = Diff(opt$value,value.old)
    #         value.old = opt$value
    #         startpars = opt$par
    #         if(MESSAGE & step %% 10 == 0){
    #             cat("This is the ", step, "th step. The optimum value is ",
    #                 opt$value," after ",opt$counts[1]," iterations \n",sep="")   
    #         }
    #     }
    #     
    #     opt = list(par = opt$par, convergence = (step < max.step), 
    #                value = opt$value)
    #     return(opt)
    #     
    # }
    # 
    
    
    
  