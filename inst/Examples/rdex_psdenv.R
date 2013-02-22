#RDEX#\dontrun{
##
## psd working environment
##
# Get some status information about the .psdenv environment
psd_envStatus()
#
# Get a list of all variables in .psdenv
psd_envList()
#
# Pull the variable "init" into .GlobalEnv
print(x <- psd_envGet("init"))
#
# pull the adaptive history into .GlobalEnv
get_adapt_history()
#RDEX#}
