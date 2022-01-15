# Simple R helper functions for use with rpy2
#

# adds paths to R's environmental variables PATH and PKG_CONFIG_PATH
set_ev = function(add_path=NULL, add_pconfig=NULL) {

        res_path = res_pconfig = NULL

        if( !is.null(add_path) )
        {
            existing_path = Sys.getenv('PATH')
            if(existing_path=='') existing_path = NULL
            new_path = paste(existing_path, add_path, sep=':')
            res_path = Sys.setenv('PATH' = new_path)
        }
        
        if( !is.null(add_pconfig) )
        {
            existing_pconfig = Sys.getenv('PKG_CONFIG_PATH')
            if(existing_pconfig=='') existing_pconfig = NULL
            new_pconfig = paste(existing_pconfig, add_pconfig, sep=':')
            res_pconfig = Sys.setenv('PKG_CONFIG_PATH' = new_pconfig)
        } 

        return(list(PATH=res_path, PKG_CONFIG_PATH=res_pconfig))
    }

# reshape input point data vectors (list) into snapped grid configuration
get_sglist = function(vec_in, gsnap) {

    lapply(vec_in, \(x) x[gsnap[["sg"]][["ipvec"]]])
}

# call kriging function (on list of vectors) with garbage collection 
run_krig = function(Rsgdata, Rgsnap, Rpvario) {

    pc = pkern_cmean(Rsgdata[[1]], Rgsnap, Rpvario, pc=TRUE)

    zout = lapply(Rsgdata, \(x) {
        
        z = pkern_cmean(x, Rgsnap, Rpvario, pc=pc)
        gc()
        return(z)
        
    })

    gc()
    return(zout)
    
}