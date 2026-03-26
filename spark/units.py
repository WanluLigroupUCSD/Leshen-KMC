"""Physical constants and unit conversions (SI base)."""

import numpy as np

# Fundamental constants
kB = 1.380649e-23          # Boltzmann constant [J/K]
h = 6.62607015e-34         # Planck constant [J·s]
hbar = h / (2 * np.pi)    # Reduced Planck constant [J·s]
NA = 6.02214076e23         # Avogadro number [1/mol]
R_gas = kB * NA            # Gas constant [J/(mol·K)]
e_charge = 1.602176634e-19 # Elementary charge [C]

# Energy
eV = e_charge              # 1 eV in Joules
kJ_per_mol = 1e3 / NA
kcal_per_mol = 4184.0 / NA

# Length
angstrom = 1e-10           # 1 Angstrom in meters

# Pressure
bar = 1e5                  # 1 bar in Pa
atm = 101325.0             # 1 atm in Pa

# Mass
umass = 1.66053906660e-27  # 1 amu in kg

# Math
pi = np.pi

# Atomic masses (amu)
ATOMIC_MASSES = {
    'H': 1.008, 'He': 4.003, 'Li': 6.941, 'C': 12.011,
    'N': 14.007, 'O': 15.999, 'F': 18.998, 'S': 32.06,
    'Cl': 35.45, 'Ar': 39.948, 'Fe': 55.845,
    'Mo': 95.95, 'W': 183.84, 'Pt': 195.08, 'Pd': 106.42,
}

# Molecular masses (amu)
MOLECULAR_MASSES = {
    'H2': 2.016, 'N2': 28.014, 'O2': 31.998,
    'CO': 28.010, 'CO2': 44.009, 'H2O': 18.015,
    'NH3': 17.031, 'N2H4': 32.045, 'NNH': 29.022,
    'HNNH': 30.030, 'N2H2': 30.030,
}
