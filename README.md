# Electrical Circuit Simulator

An object-oriented electrical circuit simulator implemented in C++, supporting DC, transient, AC/frequency, and phase-domain analyses through a Modified Nodal Analysis (MNA) formulation.

The project combines a circuit simulation engine with a command-line interface and an SDL-based graphical interface for schematic construction, analysis, and result visualization.

## Overview

The simulator represents electrical circuits as interconnected nodes and elements and automatically constructs the corresponding system of equations using **Modified Nodal Analysis (MNA)**.

Different circuit elements contribute to the system through element-specific **stamping** routines. The resulting linear systems are solved numerically, with complex-valued systems used for AC analysis.

The project evolved in two phases, with the second phase extending the original circuit-analysis engine with frequency/phase analysis and a graphical schematic and plotting environment.

## Features

### Circuit Elements

The simulator supports:

* **Passive components**

  * Resistor
  * Capacitor
  * Inductor

* **Independent voltage sources**

  * DC
  * AC
  * Sinusoidal
  * Pulse

* **Independent current sources**

  * DC
  * AC
  * Sinusoidal
  * Pulse

* **Dependent sources**

  * Voltage-Controlled Voltage Source (VCVS)
  * Voltage-Controlled Current Source (VCCS)
  * Current-Controlled Voltage Source (CCVS)
  * Current-Controlled Current Source (CCCS)

### Circuit Analysis

The simulator provides several analysis modes:

* **DC Analysis**
* **DC Sweep**
* **Transient Analysis**
* **AC / Frequency Sweep**
* **Phase Sweep**

Voltage and current quantities can be selected as analysis outputs, and multiple quantities can be plotted when supported by the analysis mode.

### Numerical Methods

The circuit equations are formulated using **Modified Nodal Analysis (MNA)**.

Key parts of the numerical implementation include:

* Automatic node indexing
* Additional current unknowns required by MNA
* Element-specific matrix stamping
* Real-valued MNA systems for DC-related analyses
* Complex-valued MNA systems for AC analysis
* Gaussian elimination for solving linear systems
* Phasor-based voltage and current extraction in AC analysis

For AC analysis, the simulator constructs a complex MNA system as a function of angular frequency and solves it for each frequency point.

## Graphical Interface

The project includes an SDL-based graphical interface providing:

* Schematic construction
* Component placement
* Node editing and merging
* Circuit saving and loading
* Dedicated interfaces for different analysis modes
* Plotting of simulation results
* Interactive data cursor functionality

The graphical interface is designed to provide a visual workflow on top of the underlying circuit simulation engine.

## Command-Line Interface

In addition to the graphical interface, the simulator provides a command-line workflow for defining circuits and requesting analyses.

Analysis commands include operations corresponding to:

* DC sweep
* Transient analysis
* AC sweep
* Phase sweep
* Plotting simulation results
* Matrix inspection

This allows the numerical engine to be used independently of the graphical interface.

## Architecture

The project follows an object-oriented design centered around the circuit and element abstractions.

At the core are:

* `Circuit` — manages the circuit topology, node indexing, MNA matrices, and analysis state.
* `Node` — represents circuit nodes and their connectivity.
* `Element` — abstract base class for circuit elements.
* Derived element classes — implement the electrical behavior and matrix-stamping operations of individual components.

The interface and analysis layers are built around the simulation core:

```text
                 ┌──────────────────────┐
                 │   User Interface     │
                 │                      │
                 │  CLI / SDL GUI       │
                 └──────────┬───────────┘
                            │
                            ▼
                 ┌──────────────────────┐
                 │   Analysis Layer     │
                 │                      │
                 │ DC / DC Sweep        │
                 │ Transient            │
                 │ AC / Frequency Sweep │
                 │ Phase Sweep          │
                 └──────────┬───────────┘
                            │
                            ▼
                 ┌──────────────────────┐
                 │   Circuit Model      │
                 │                      │
                 │ Circuit / Node       │
                 │ Element hierarchy    │
                 └──────────┬───────────┘
                            │
                            ▼
                 ┌──────────────────────┐
                 │    MNA Formulation   │
                 │                      │
                 │ Matrix Stamping      │
                 │ Linear System Solve  │
                 └──────────────────────┘
```

## Project Evolution

### Phase 1

The first phase focused on building the circuit simulation engine and implementing the fundamental circuit-analysis functionality, including:

* Object-oriented circuit representation
* Circuit file parsing
* Passive and active circuit elements
* Dependent sources
* MNA formulation
* DC analysis
* DC sweep
* Transient analysis
* Initial frequency-domain functionality

### Phase 2

The second phase extended the simulator with a graphical environment and additional frequency-domain functionality, including:

* AC analysis and frequency sweep
* Phase sweep
* Interactive plotting
* SDL-based schematic editor
* Component placement and circuit construction
* Node editing and merging
* Schematic saving and loading
* GUI integration of the different analysis modes

Both phases are preserved in the Git history of this repository.

## Build

The project uses **CMake** for build configuration.

A typical build workflow is:

```bash
mkdir build
cd build
cmake ..
cmake --build .
```

The exact runtime dependencies and configuration may depend on the development environment and SDL installation.

## Technologies

* **C++**
* **CMake**
* **SDL**
* **Object-Oriented Programming**
* **Modified Nodal Analysis**
* **Gaussian Elimination**
* **Complex Arithmetic**
* **Numerical Circuit Simulation**

## Project Context

This project was developed as an object-oriented electrical circuit simulation project and subsequently extended with a graphical analysis and schematic-editing environment.

The repository consolidates the two development phases into a single project while preserving their development histories.

## Author

**Soheil Bolvardi**

Electrical Engineering & Mathematics Student
