import numpy as np

# Input values of r1, r2, and r3 from the user
while True:
    r1_input = input("Enter the position vector r1 (comma-separated values in km): ").split(',')
    if all(r1_input):
        break
    else:
        print("Please provide values for all components of r1.")

while True:
    r2_input = input("Enter the position vector r2 (comma-separated values in km): ").split(',')
    if all(r2_input):
        break
    else:
        print("Please provide values for all components of r2.")

while True:
    r3_input = input("Enter the position vector r3 (comma-separated values in km): ").split(',')
    if all(r3_input):
        break
    else:
        print("Please provide values for all components of r3.")

# Convert input strings to numpy arrays
r1 = np.array([float(x) for x in r1_input])
r2 = np.array([float(x) for x in r2_input])
r3 = np.array([float(x) for x in r3_input])

# Calculate vectors and magnitudes
r1_mag = np.linalg.norm(r1)
r2_mag = np.linalg.norm(r2)
r3_mag = np.linalg.norm(r3)

# Calculate delta terms
delta1 = np.cross(r2, r3)
delta2 = np.cross(r3, r1)
delta3 = np.cross(r1, r2)

delta1_mag = np.linalg.norm(delta1)
delta2_mag = np.linalg.norm(delta2)
delta3_mag = np.linalg.norm(delta3)

# Calculate N vector (in km^3)
N = r1_mag * delta1 + r2_mag * delta2 + r3_mag * delta3
N_mag = np.linalg.norm(N)

# Calculate D vector (in km^2)
D = delta1 + delta2 + delta3
D_mag = np.linalg.norm(D)

# Calculate S vector (in km^2)
S = r1 * (r2_mag - r3_mag) + r2 * (r3_mag - r1_mag) + r3 * (r1_mag - r2_mag)

# Print the values
print("Magnitude of r1:", r1_mag)
print("Magnitude of r2:", r2_mag)
print("Magnitude of r3:", r3_mag)

print("Delta 1:", delta1)
print("Delta 2:", delta2)
print("Delta 3:", delta3)

print("Magnitude of Delta 1:", delta1_mag)
print("Magnitude of Delta 2:", delta2_mag)
print("Magnitude of Delta 3:", delta3_mag)

print("N vector (km^3):", N)
print("Magnitude of N (km^3):", N_mag)

print("D vector (km^2):", D)
print("Magnitude of D (km^2):", D_mag)

print("S vector (km^2):", S)

# Define the gravitational parameter (mu)
mu = 398600  # Gravitational parameter of Earth in km^3/s^2

# Calculate v2
v2 = np.sqrt(mu / (N_mag * D_mag)) * ((np.cross(r2, D) / r2_mag) - S)

# Print the calculated v2
print("Velocity at time t2:", v2)

# Define a function to calculate classical orbital elements from the velocity vector
def classical_orbital_elements(r2, v2, mu):
    # Calculate specific angular momentum vector
    h = np.cross(r2, v2)
    h_mag = np.linalg.norm(h)
    
    # Calculate eccentricity vector
    e_vector = ((np.linalg.norm(v2)**2 - mu/np.linalg.norm(r2)) * r2 - np.dot(r2, v2) * v2) / mu
    e = np.linalg.norm(e_vector)
    
    # Calculate inclination (i)
    i = np.arccos(h[2] / h_mag)
    
    # Calculate the node vector (n)
    k = np.array([0, 0, 1])
    n = np.cross(k, h)
    n_mag = np.linalg.norm(n)
    
    # Calculate right ascension of the ascending node (Ω)
    Omega = np.arctan2(n[1], n[0])
    if Omega < 0:
        Omega += 2 * np.pi
    
    # Calculate argument of periapsis (ω)
    if n_mag != 0:
        omega = np.arccos(np.dot(n, e_vector) / (n_mag * e))
        if e_vector[2] < 0:
            omega = 2 * np.pi - omega
    else:
        omega = 0
    
    # Calculate true anomaly (ν)
    nu = np.arccos(np.dot(e_vector, r2) / (e * np.linalg.norm(r2)))
    if np.dot(r2, v2) < 0:
        nu = 2 * np.pi - nu
    
    # Calculate semi-major axis (a)
    E = np.linalg.norm(v2)**2 / 2 - mu / np.linalg.norm(r2)
    a = -mu / (2 * E)
    
    # Convert angles from radians to degrees
    i = np.degrees(i)
    Omega = np.degrees(Omega)
    omega = np.degrees(omega)
    nu = np.degrees(nu)
    
    return a, e, i, Omega, omega, nu

# Calculate classical orbital elements
a, e, i, Omega, omega, nu = classical_orbital_elements(r2, v2, mu)

# Print the classical orbital elements
print("Semi-major axis (a):", a, "km")
print("Eccentricity (e):", e)
print("Inclination (i):", i, "degrees")
print("Right Ascension of the Ascending Node (Ω):", Omega, "degrees")
print("Argument of Periapsis (ω):", omega, "degrees")
print("True Anomaly (ν):", nu, "degrees")
