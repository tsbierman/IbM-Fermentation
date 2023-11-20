# The sole purpose of this file is to initialise the Java Virtual Machine (JVM)
# This only needs to run once, so it wont be necessary to comment when create_mat is run multiple thermodynamic_parameters

using JavaCall

JavaCall.addClassPath(string(pwd(), "\\lib\\bacteria\\shovingQuadTreekDist.jar"))
JavaCall.addClassPath(string(pwd(), "\\lib\\bacteria\\Results.java"))
JavaCall.init()
println(">>>>>>>>>>> DONE INITIALISING!")
