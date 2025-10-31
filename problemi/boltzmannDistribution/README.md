# Boltzmann Distribution simulation

This project simulates the evolution of a N molecules system, with the same initial energy equal to $3/2 k_BT$:

- At each interaction, a randomly chosen molecule transfers a quanta of energy $\Delta U = (k_BT)/h$ to another randomly chosen molecule. $h\approx 10$.
- After each interaction, an histogram shows the energy distribution of the molecules.

The number of molecules and the number of interactions are user defined parameters. $k_B T$ and $h$ are fixed to 2.0 and 10.0 respectively (eventually can be also set by the user).

## Usage

```bash
boltzmannDistribution.out <num_molecules> <interactions [list comma separated]> 
```

for example:

```bash
boltzmannDistribution.out 1000 1000,10000,100000
```

## Plotting

Two possible plotting solutions are provided:

- Plain stdout based plotting (ASCII);
- ROOT based plotting (requires ROOT installed).