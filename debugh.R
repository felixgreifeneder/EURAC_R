sourceDebugging<-function(f){
 #Function to inject the code to
 theCode<-function(){}
 #Injection
 parse(text=c('{',readLines(f),'}'))->body(theCode)
 #Triggering debug
 debug(theCode)
 #Lift-off
 theCode()
}
sourceDebugging(X:/Workspaces/GrF/05_Code/R/netcdft2csv.R)