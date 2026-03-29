.doDebugLogging <- function(){
  return(Sys.getenv('NIMBLE_DEGUG_LOGGING') == 1)
}