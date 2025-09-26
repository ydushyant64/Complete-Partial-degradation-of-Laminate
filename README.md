# Complete-Partial-degradation-of-Laminate
Developed MATLAB codes for simulating complete and partial degradation of glass epoxy laminates under combined me- chanical and hygro-thermal loading conditions

# omechanical and Hygro-Thermal Degradation of Composite Laminates in MATLAB

This repository contains MATLAB scripts for analyzing the progressive failure of composite laminates under mechanical and hygro-thermal loads. The project implements two distinct degradation models—**Complete Degradation** and **Partial Degradation**—to simulate how a glass/epoxy laminate behaves as its individual layers (plies) begin to fail.

---

## The Big Idea: Why Does This Matter?

Composite materials, like the glass/epoxy laminates studied here, are incredibly strong and lightweight, making them essential in aerospace, automotive, and high-performance sports equipment. However, predicting how and when they will fail is complex. Unlike metals, which tend to fail all at once, composites fail progressively, one layer at a time.

This project tackles that complexity head-on. The goal was to build a computational tool in MATLAB that could:
-   Simulate the stresses on each individual layer of a composite laminate.
-   Predict which layer will fail first (an event known as "First Ply Failure").
-   Analyze how the load is redistributed among the remaining layers after a failure occurs.
-   Compare two different assumptions about what happens after a layer fails: Does it become completely useless (**Complete Degradation**), or does it only lose strength in a specific direction (**Partial Degradation**)?

---

## How It's Built: The Simulation Process

The analysis is built on the principles of **Classical Lamination Theory**. The MATLAB scripts automate the entire process, from defining the material properties to iteratively checking for failure.

#### 1. Defining the Material and Laminate
-   **Material:** The project uses standard properties for a **Glass/Epoxy** unidirectional lamina, a common composite material.
-   **Laminate Stacking Sequence:** The laminate is defined with a specific stacking sequence of `[0, 45, -45, 90, 90, -45, 45, 0]`, which is a common quasi-isotropic layup designed to provide strength in multiple directions.
-   **Loading Conditions:** The simulation applies both a mechanical load (uniaxial tension) and hygro-thermal loads (changes in temperature and moisture), which create internal stresses in the laminate.

#### 2. The Core of the Analysis: Stress Calculation
The heart of the simulation is a loop that calculates the stress on each ply. It constructs the **[ABD] matrix**, a cornerstone of lamination theory that relates the applied loads to the strains and curvatures of the laminate. From there, the stresses on each individual ply are calculated in its own fiber direction.

#### 3. Predicting Failure: The Maximum Stress Criterion
To determine if a ply has failed, the simulation uses the **Maximum Stress Failure Criterion**. This theory compares the calculated stresses in each direction (longitudinal, transverse, and shear) against the known strength limits of the material. The **stress ratio (SR)** is calculated for each mode, and if any SR value reaches or exceeds 1, that ply is considered to have failed.

---

## The Two Degradation Models

This is where the project gets interesting. The two MATLAB scripts explore different "what-if" scenarios for post-failure behavior:

#### a) `complete_degradation.m`
-   **The Assumption:** This model assumes a catastrophic failure. Once a ply fails in any mode (tension, compression, or shear), all of its material properties (E1, E2, G12) are instantly set to zero.
-   **The Outcome:** This is a conservative approach, representing a "worst-case" scenario. The load that was carried by the failed ply is immediately redistributed to the remaining plies, which can often cause a rapid, cascading failure of the entire laminate.

#### b) `partial_degradation.m`
-   **The Assumption:** This is a more realistic and nuanced model. It assumes that a ply only loses stiffness in the specific direction of its failure. For example:
    -   If it fails in the transverse direction, only E2 and G12 are set to zero, but it can still carry some load in the fiber direction (E1).
    -   If it fails in shear, only G12 is set to zero.
-   **The Outcome:** This approach allows for a more gradual and realistic simulation of damage progression. The laminate can sustain more load after the First Ply Failure, providing a better estimate of its ultimate strength.

---

## How to Use the Code

1.  Open either `complete_degradation.m` or `partial_degradation.m` in MATLAB.
2.  The material properties, laminate stacking sequence, and loading conditions are all defined at the beginning of the script. You can modify these to simulate different scenarios.
3.  Run the script.
4.  The MATLAB command window will output a step-by-step analysis of the failure progression, including the stress ratios in each ply for each iteration and the load at which each ply failure occurs.
