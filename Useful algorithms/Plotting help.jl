a = [1, 2, 3, 4, 5, 6]

b = repeat(a, 6)
c = transpose(repeat(a, 6))


x = transpose(repeat(a, 1, 6))
y = transpose(repeat(a, inner=(1,6)))

tempX = nodeMidCoordinatesX[1:5]
x = repeat(tempX,inner=4)
x = repeat(tempX,4)
println(x)
