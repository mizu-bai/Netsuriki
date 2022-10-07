module Translational

using Molecules
using Unitful
import PhysicalConstants.CODATA2018: R, k_B, N_A, h


@doc raw"""
    Translational.q(mol::Molecule, T::Unitful.Temperature=298.15u"K", p::Unitful.Pressure=1u"atm")::Float64

Compute **Translational Partition Function**.

# Detail

```math
q_{tr} = {\left( \frac{2 \pi m k_B T}{h^2} \right)}^{\frac{3}{2}} V
```

where ``m`` is the mass of molecule, ``k_B`` is **Boltzmann Constant**, ``h`` is **Plank Constant**,
``T`` is temperature, ``p`` is pressure, and ``V`` is the volume of 1 mol molecule in gas phase
under temperature ``T`` and pressure ``p``.

# Arguments

- `mol::Molecule`: a `Molecule` object
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
- `p::Unitful.Pressure`: Pressure, with default value of `1u"atm"`
"""
function q(mol::Molecule, T::Unitful.Temperature=298.15u"K", p::Unitful.Pressure=1u"atm")::Float64
    T = T |> u"K"
    p = p |> u"Pa"
    m = [atom.mass for atom in mol.atoms] |> sum |> m -> m * (1u"g/mol" |> u"kg/mol") / N_A
    V = T |> T -> 1.0u"mol" * R * T / p
    (2 * π * m * k_B * T / h^2)^(3/2) * V |> Unitful.NoUnits
end


@doc raw"""
    Translational.Am(q_tr::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Translational Molar Helmholtz Free Energy**.

# Detail

```math
\begin{align}
A_{m, tr} = & -k_B T \ln \frac{{(q_{tr})}^{N_A}}{N!} \\
          = & - N_A k_B T (\ln q_{tr} - \ln N_A + 1) \\
          = & - R T (\ln q_{tr} - \ln N_A + 1) 
\end{align}
```

# Arguments

- `q::Float64`: **Translational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Am(q::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    -R * T * (log(q) - log(N_A |> ustrip) + 1)
end


@doc raw"""
    Translational.Gm(q::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Translational Molar Gibbs Free Energy**.

# Detail

```math
G_{m, tr} = A_{m, tr} + pV = A_{m, tr} + RT
```

# Arguments

- `q::Float64`: **Translational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Gm(q::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    Am(q) + R * T
end

@doc raw"""
    Translational.Um(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Translational Molar Internal Thermal Energy**.

# Detail

```math
G_{m, tr} = \frac{3}{2} N k_B T = \frac{3}{2} R T
```

# Arguments

- `q::Float64`: **Translational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Um(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    3/2 * R * T
end

@doc raw"""
    Translational.Hm(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Translational Molar Enthalpy**.

# Detail

```math
H_{m, tr} = U_{m, tr} + p V = U_{m, tr} + R T = \frac{3}{2} R T + R T = \frac{5}{2} R T
```

# Arguments

- `q::Float64`: **Translational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Hm(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    5/2 * R * T
end


@doc raw"""
    Translational.Sm(q::Float64)::typeof(1.0u"J/mol/K")

Compute **Translational Entropy**.

# Detail

```math
S_{m, tr} = N_A k_B [ \ln q_{tr} - \ln N_A + \frac{5}{2} ] = R [ \ln q_{tr} - \ln N_A + \frac{5}{2} ]
```

# Arguments

- `q::Float64`: Electronic Partition Function.
"""
function Sm(q::Float64)::typeof(1.0u"J/mol/K")
    R * (log(q) - log(N_A |> ustrip) + 5/2) 
end


@doc raw"""
    Translational.CVm()::typeof(1.0u"J/mol/K")

Compute **Translational Molar Volume Heat Capacity**.

# Detail

```math
C_{Vm, tr} = \frac{3}{2} N_A k_B = \frac{3}{2} R
```
"""
function CVm()::typeof(1.0u"J/mol/K")
    3/2 * R
end

@doc raw"""
    Translational.Cpm()::typeof(1.0u"J/mol/K")

Compute **Translational Molar Pressure Heat Capacity**.

# Detail

```math
C_{pm, tr} = \frac{5}{2} N_A k_B = \frac{5}{2} R
```
"""
function Cpm()::typeof(1.0u"J/mol/K")
    5/2 * R
end

end