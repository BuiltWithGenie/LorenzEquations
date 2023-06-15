## Simulation of the Lorenz equations

Dashboard showing a real-time simulation of the Lorenz equations.

![preview](lorenz.gif)

## Installation

Clone the repository and install the dependencies:

First `cd` into the project directory then run:

```bash
$> julia --project -e 'using Pkg; Pkg.instantiate()'
```

Then run the app

```bash
$> julia --project
```

```julia
julia> using GenieFramework
julia> Genie.loadapp() # load app
julia> up() # start server
```

## Usage

Open your browser and navigate to `http://localhost:8000/`

---

Alternatively, register/login at <https://geniecloud.app> and import the app. You can then deploy it to the cloud.
