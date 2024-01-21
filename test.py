from scipy.optimize import fsolve
import numpy as np

# Define the inequality function
def inequality_function(k):
    return 1 - (1 - 3/k)**k - 0.1

# Initial guess for k
initial_guess = 10.0

# Use fsolve to find the root of the inequality function
result = fsolve(inequality_function, initial_guess)

# Result will contain the value of k that satisfies the inequality
print("Approximate value of k:", result[0])
k=150
print((1 - 3/k)**k)