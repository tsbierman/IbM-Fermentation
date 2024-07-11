"""
This file initialises the Java Virtual Machine (JVM).
This file only needs to un once per Julia session.
Running it more than once gives an error and a warning: "JVM already initialised"
This file is already run by the inclusion_file in general usage.
While developing, run this file by include("lib\\pre_processing\\initialiseJVM.jl")
"""

using JavaCall

JavaCall.addClassPath(string(pwd(), "\\lib\\bacteria\\shovingQuadTreekDist.jar"))
JavaCall.addClassPath(string(pwd(), "\\lib\\bacteria\\Results.java"))
JavaCall.init()
println(">>>>>>>>>>>>>>JVM INITIALISED!")
