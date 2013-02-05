\dontrun{
##
## rlpSpec working environment
##
# Get some status information about the .rlpenv environment
rlp_envStatus()
#
# Get a list of all variables in .rlpenv
rlp_envList()
#
# Pull the variable "init" into .GlobalEnv
print(x <- rlp_envGet("init"))
#
# pull the adaptive history into .GlobalEnv
get_adapt_history()
}
