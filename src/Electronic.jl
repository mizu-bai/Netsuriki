module Electronic

using Molecules
using Unitful
import PhysicalConstants.CODATA2018: R


@doc raw"""
    Electronic.q(mol::Molecule)::Float64

Compute **Electronic Partition Function**, only ground state is considered.

# Detail

```math
q_{el} = g_{el}
```

where ``g_{el}`` is the degeneracy of the ground state energy level.

# Arguments

- `mol::Molecule`: A `Molecule` object
"""
function q(mol::Molecule)::Float64
    mol.multiplicity |> float |> Unitful.NoUnits
end


@doc raw"""
    Electronic.Am(q::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute Electronic Helmholtz Free Energy.

# Detail

```math
A_{m, el} = -N_A k_B \ln q_{el} = -RT \ln q_{el}
```

where ``R`` is **Molar Gas Constant**, ``T`` is temperature, and ``q_{el}`` is **Electronic Partition Function**.

# Arguments

- `q::Float64`: **Electronic Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Am(q::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    -R * T * log(q)
end


@doc raw"""
    Electronic.Gm(q::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute Electronic **Gibbs Free Energy**.

# Detail

```math
G_{m, el} = -N_A k_B \ln q_{el} = -RT \ln q_{el}
```

where ``R`` is **Molar Gas Constant**, ``T`` is temperature, and ``q_{el}`` is **Electronic Partition Function**.

# Arguments

- `q::Float64`: **Electronic Partition Function**
- `T::Unitful.Temperature`: Temperature, `298.15u"K"` by default
"""
function Gm(q::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    -R * T * log(q)
end


@doc raw"""
    Electronic.Um()::typeof(1.0u"J/mol")

Compute **Electronic Molar Internal Thermal**.

# Detail

For **Electronic Partition Function** has no contribution to **Molar Internal Thermal**, this function will always
return `0.0u"J/mol"`.
"""
function Um()::typeof(1.0u"J/mol")
    0.0u"J/mol"
end


@doc raw"""
    Electronic.Hm()::typeof(1.0u"J/mol")

Compute **Electronic Molar Enthalpy**.

# Detail

For **Electronic Partition Function** has no contribution to **Molar Enthalpy**, this function will always return
`0.0u"J/mol"`.
"""
function Hm()::typeof(1.0u"J/mol")
    0.0u"J/mol" 
end


@doc raw"""
    Electronic.Sm(q::Float64)::typeof(1.0u"J/mol/K")

Compute **Electronic Entropy**.

# Detail

```math
S_{el} = N_A k_B \ln q_{el} = -R \ln q_{el}
```

# Arguments

- `q::Float64`: Electronic Partition Function.
"""
function Sm(q::Float64)::typeof(1.0u"J/mol/K")
    R * log(q)
end


@doc raw"""
    Electronic.CVm()::typeof(1.0u"J/mol/K")

Compute **Electronic Molar Volume Heat Capacity**.

# Detail

For **Electronic Partition Function** has no contribution to **Molar Volume Heat Capacity**, this function will always
return `0.0u"J/mol/K"`.
"""
function CVm()::typeof(1.0u"J/mol/K")
    0.0u"J/mol/K"
end


@doc raw"""
    Electronic.Cpm()::typeof(1.0u"J/mol/K")

Compute **Electronic Molar Pressure Heat Capacity**.

# Detail

For **Electronic Partition Function** has no contribution to **Molar Pressure Heat Capacity**, this function will always
return `0.0u"J/mol/K"`.
"""
function Cpm()::typeof(1.0u"J/mol/K")
    0.0u"J/mol/K"
end

end
