#RDEX#\dontrun{
require(psd)
##
## psd working environment
##
# Get some status information about the psd working environment
psd_envStatus()
#
# Get a list of all variables
psd_envList()
#
# Pull the variable "init" into .GlobalEnv
print(x <- psd_envGet("init"))
#
# Pull the adaptive history into .GlobalEnv
get_adapt_history()
#RDEX#}
