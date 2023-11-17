# pH testing
# Case 1, denominator fits 
Ct = 30
pH = 6.0
H = 10^(-pH)
k1 = 5.27e-4
k2 = 5.27e-7
theta = H^2 + k1 * H + k1*k2
B = H^2 * Ct/theta
B1 = H * Ct * k1/theta
B2 = k1*k2*Ct/theta
println("First case")
println("B: $(B)")
println("B1: $(B1)")
println("B2: $(B2)")

k3 = 0
theta2 = H^3 + H^2 * k1 + H * k1 * k2 + k1 * k2 * k3
B = H^3 * Ct/theta2
B1 = Ct * k1 * H^2/theta2
B2 = Ct * k1 * k2 * H/theta2
B3 = k1 * k2 * k3 * Ct/theta2

println("Second case")
println("B: $(B)")
println("B1: $(B1)")
println("B2: $(B2)")
println("B3: $(B3)")
