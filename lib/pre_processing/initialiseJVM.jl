# The sole purpose of this file is to initialise the Java Virtual Machine (JVM)
# This file only needs to run once per Julia session, to initialise the JVM
# Running it more than once gives both an error and a warning "JVM already initialised"

using JavaCall

JavaCall.addClassPath(string(pwd(), "\\lib\\bacteria\\shovingQuadTreekDist.jar"))
JavaCall.addClassPath(string(pwd(), "\\lib\\bacteria\\Results.java"))
JavaCall.init()
println(">>>>>>>>>>>>>>JVM INITIALISED!")
