using OptimalDesign # this already imports Pumas
using PumasUtilities
using DataFramesMeta
using PharmaDatasets

### Warfarin - anticoagulant (https://en.wikipedia.org/wiki/Warfarin)
# data from Bauer et. al 
# Tutorial for $DESIGN in NONMEM: Clinical trial evaluation and optimization
# DOI: 10.1002/psp4.12713
datapath = dataset("warfarin", String)
df = CSV.read(datapath, DataFrame; missingstring=[".", "", "NA"])
pop = read_pumas(df)

## create a :route column for NCA
@rtransform! df :route = "ev"

## NCA Population
pop_nca = read_nca(
  df;
  observations=:dv)
## Mean Concentration vs Time Plot
summary_observations_vs_time(
  pop_nca,
  axis=(; xlabel = "Time (hr)",
          ylabel = "Warfarin Concentration (μg/mL)"))

### Step 1 - Model
# 1-cmt oral model
model = @model begin
  @param begin
    θ₁ ∈ RealDomain(lower=0.0, init=0.15)
    θ₂ ∈ RealDomain(lower=0.0, init=8.0)
    θ₃ ∈ RealDomain(lower=0.0, init=1.0)
    Ω  ∈ PSDDomain(3)
    σ  ∈ RealDomain(lower=0.0001, init=sqrt(0.01))
  end
  @random begin
    η ~ MvNormal(Ω)
  end
  @pre begin
    Tvcl = θ₁
    Tvv  = θ₂
    Tvka = θ₃
    CL   = Tvcl*exp(η[1])
    Vc   = Tvv*exp(η[2])
    Ka   = Tvka*exp(η[3])
  end
  @dynamics Depots1Central1
  @vars begin
    conc = Central / Vc
  end
  @derived begin
    dv ~ @. Normal(log(conc), σ)
  end
end

params = (
  θ₁ = 0.15,
  θ₂ = 8.0,
  θ₃ = 1.0,
  Ω  = [0.07 0.0 0.0;
        0.0 0.02 0.0;
        0.0 0.0 0.6],
  σ  = sqrt(0.01),
)

### Interlude - Date and time API in Julia
## `Date` constructors
# 25 Jan 2021
Date(2021, 1, 25)
Date("2021-01-25")

## `DateTime` constructors
# 25 Jan 2021 - 9:00 am
t0 = DateTime(2021, 1, 25, 9)  # also could be `DateTime("2021-01-25T09:00")`
# 25 Jan 2021 - 9:00 pm
tend = DateTime(2021, 1, 25, 21)

## `DateTime` calendar arithmetic with `Period` types
# Next second
t0 + Second(1)
# Next minute
t0 + Minute(1)
# Next hour
t0 + Hour(1)
# Next day
t0 + Day(1)
# Next week
t0 + Day(7)
t0 + Week(1)
# Next month
t0 + Month(1)
# Next year
t0 + Year(1)
# It also works for mix and match `Period` types
t0 + Second(1) + Minute(1) + Hour(1) + Day(1) + Week(1) + Month(1) + Year(1)

## Intervals for time bounds using the `..` operator
bound1 = t0..tend
# Careful with syntax, use `()`
t0..tend + Day(1)   # just the finishing bound day
(t0..tend) + Day(1) # both start and finishing bounds

### Fixed samples per time window
# You can pass a `Dict` with:
# - keys: Time bound interval
# - values: number of samples in that interval
bounds = Dict(
  (t0..tend) => 15,
  (t0..tend) + Day(1) => 10,
)


### Step 2 - Create a decision with constraints
## Decision variables
Nsubjects = length(pop)
## Decision
dec = decision(
  model, pop, params;            # Pumas API model, pop, params 
  type=:observation_times,     # only type supported currently
  bounds=[bound1],             # a `Vector` or `Dict` of time intervals
  # bounds=bounds,             # (if `Dict` adjust `N` argument properly)
  N=15,                        # total number of samples per subject
  subject_multiples=fill(3, Nsubjects), # the number of times each subject is counted (default 1)
  minimum_offset=Minute(30),   # enforce a minimum duration between any 2 samples
  model_time_unit=Hour(1),     # unit of time assumed in the model definition and dynamics model
)

### Step 3 - Optimize the design
# :aoptimal: minimizes the trace of the inverse of the expected information matrix.
# :doptimal: maximizes the (log) determinant of the expected information matrix.
# :toptimal: maximizes the trace of the expected information matrix.
@time result = design(
  dec;
  optimality=:doptimal,
  time_limit=400.0,             # time limit in seconds, maximum 20min
  nlp_tol=1e-4,                 # tolerance for the nonlinear solver
  verbose=true,                 # display progress and information?
  processors=Threads.nthreads() # multi-threaded and how many? (default 1)
)

### Inspect the result
# Optimal sample times for each subject
optimaltimes(result)
# Initial objective value
initvalue(result)
# Optimal objective value
optimalvalue(result)
# D-Efficieny of the FIM
# (det(FIM_2) / det(FIM_1))^(1/M) where:
# - FIM_2 is the Fisher information matrix (FIM) at the optimal design
# - FIM_1 is the FIM at the initial design
# - M is the number of model parameters
defficiency(result)
