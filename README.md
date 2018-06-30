# Universe age

This program calculates the age of the universe by integrating Friedmann equation.
![friedmann_eq.png](https://raw.githubusercontent.com/shu-yusa/universe-age/master/friedmann_eq.png)

The following parameters are fixed.
* H<sub>0</sub> = 72 km/s/Mpc
* &Omega;<sub>R</sub> = 0
* a<sub>0</sub> = 1

## Parameter cases

* (&Omega;<sub>M</sub>, &Omega;<sub>K</sub>, &Omega;<sub>&Lambda;</sub>) = (1, 0, 0)
* (&Omega;<sub>M</sub>, &Omega;<sub>K</sub>, &Omega;<sub>&Lambda;</sub>) = (0.3, 0.7, 0)
* (&Omega;<sub>M</sub>, &Omega;<sub>K</sub>, &Omega;<sub>&Lambda;</sub>) = (2, -1, 0)
* (&Omega;<sub>M</sub>, &Omega;<sub>K</sub>, &Omega;<sub>&Lambda;</sub>) = (3, -2, 0)

## Requirements
Fortran compiler is required.
By using `gnuplot`, you can output graph for `a` vs `t`.
```bash
brew install gfortran
brew install gnuplot
```

## Compile and Run
```bash
gfortran universe.f90
./a.out
gnuplot < gnu
```

## Results
![graph.png](https://raw.githubusercontent.com/shu-yusa/universe-age/master/graph.png)

|  (&Omega;<sub>M</sub>, &Omega;<sub>K</sub>, &Omega;<sub>&Lambda;</sub>) |  Age  |
| ---- | ---- |
| (1, 0, 0)  |  9.07 Gyr  |
| (0.3, 0.7, 0) |  11.00 Gyr  |
| (2, -1, 0) |  13.11 Gyr    |
| (3, -2, 0) |  6.98 Gyr    |

