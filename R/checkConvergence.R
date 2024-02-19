checkConvergence <- function(Conv, ratio, L2iterationReset)
{
if (Conv < ratio)
    {
       warning("One EM start did not converge: Number of maximum iterations reached", call. = FALSE)
        #cat("\n EM did not converge: Number of maximum iterations reached \n", file = "EMwarnings.txt", append = TRUE) 
        invisible(FALSE)
        
    }else
    {
        if (L2iterationReset)
        {# EM terminated after converging, but convergence happend after a reset
           warning("One EM start did not converge: Iteration terminated after reset stagnation", call. = FALSE)
            #cat("\n EM did not converge: Iteration terminated after reset stagnation \n", file = "EMwarnings.txt", append = TRUE) 
            invisible(FALSE)
        }
        else
        {
            invisible(TRUE)
        }
        
    }
}