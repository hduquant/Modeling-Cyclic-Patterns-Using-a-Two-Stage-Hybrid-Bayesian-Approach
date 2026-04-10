# Real Data Example

This folder contains a real-data example for implementing the proposed method in R using either `rjags` or `rblimp`.

## Files

- **`data.csv`**  
  A real-data example dataset.

- **`real data example using rjags.R`**  
  R code for fitting the proposed method using the `rjags` package.

- **`real data example using rblimp.R`**  
  R code for fitting the proposed method using the `rblimp` package.

## Software Options

Two implementations of the proposed method are provided:

### 1. `rjags`
The script `real data example using rjags.R` uses the R package `rjags` to fit the proposed method.

### 2. `rblimp`
The script `real data example using rblimp.R` uses the R package `rblimp` to fit the proposed method.  
The `rblimp` package is based on the software **Blimp** and is generally more user-friendly than `rjags`. However, to use `rblimp`, you must install **Blimp** on your computer in advance.

## Requirements

Depending on which script you use, you will need:

- R
- The `rjags` package for the `rjags` example
- The `rblimp` package and **Blimp** software for the `rblimp` example

## Notes

The `rblimp` implementation may be easier for users who prefer a more user-friendly interface. The `rjags` implementation may be preferable for users who want to work directly within the JAGS framework.
