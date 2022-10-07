module Rotational

using Molecules
using Molecules.Symmetry
using Unitful
import PhysicalConstants.CODATA2018: R, k_B, N_A, h, m_u

@doc raw"""
    Rotational.calcΘrot(I::typeof(1.0u"kg*m^2"))::typeof(1.0u"K")

Compute **Characteristic Rotational Temperature**.

# Detail

```math
\Theta_{rot} = \frac{h^2}{8 \pi^2 I k_B}
```

where ``k_B`` is **Boltzmann Constant**, ``h`` is **Plank Constant**, ``I`` is the moment of inertia.

# Arguments

- `I::typeof(1.0u"kg*m^2")`: The moment of inertia
"""
function calcΘrot(I::typeof(1.0u"kg*m^2"))::typeof(1.0u"K")
    h^2 / (8π^2 * I * k_B)
end


@doc raw"""
    Rotational.calcσrot(point_group::String)::Int

Compute **Rotational Symmetry Number**.

# Detail

Rotational symmetry number is computed according to the point group of molecule.

# Arguments

- `point_group::String`: Point group symbol.
"""
function calcσrot(point_group::String)::Int
    if point_group in ["C1", "Ci", "Cs", "Cinfv"]
        return 1
    elseif point_group == "Dinfh"
        return 2
    elseif point_group[1:1] == "C"
        return parse(Int, point_group[2:2])
    elseif point_group[1:1] == "D"
        return parse(Int, point_group[2:2]) * 2
    elseif point_group[1:1] == "S"
        return parse(Int, point_group[2:2]) ÷ 2
    elseif point_group[1:1] == "T"
        return 12
    elseif point_group[1:1] == "O"
        return 24
    elseif point_group[1:1] == "I"
        return 60
    else
        throw(DomainError(point_group, "no such point group"))
    end
end


@doc raw"""
    Rotational.q(mol::Molecule, T::Unitful.Temperature=298.15u"K")::Float64

Compute **Rotational Partition Function**

# Detail

- Single Atom

```math
q_{rot} = 0
```

- Linear Molecule

```math
q_{rot} = \frac{T}{\sigma_{rot} \Theta_{rot}}
````

- Nonlinear Molecule

```math
q_{rot} = \frac{\pi^\frac{1}{2} T^{\frac{3}{2}}}{\sigma_{rot}(\Theta_a \Theta_b \Theta_c)^{\frac{1}{2}}}
```

# Arguments

- `mol::Molecule`: A `Molecule` object
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function q(mol::Molecule, T::Unitful.Temperature=298.15u"K")::Float64
    T = T |> u"K"
    # if Atom
    if length(mol.atoms) == 1
        return 0.0
    end

    # here I_abc in amu * angstrom^2
    I_abc, _ =  mol.atoms |> calcmoit |> eigenmoit
    I_abc *= 1m_u * 1u"angstrom"^2 |> u"kg * m^2"
    Θ_abc = I_abc .|> calcΘrot
    point_group = mol.atoms |> find_point_group
    σ = point_group |> calcσrot

    if point_group in ["Cinfv", "Dinfh"]
        # Linear Molecule
        1 / σ * T / Θ_abc[2]
    else
        # Nonlinear Molecule
        π^(1/2) / σ * T^(3/2) / (Θ_abc |> prod)^(1/2)
    end |> Unitful.NoUnits
end


@doc raw"""
    Rotational.Am(q_rot::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Rotational Molar Helmholtz Free Energy**.

# Detail

```math
A_{m, rot} = - N_A k_B T \ln q_{rot} = - R T \ln q_{rot}
```

# Arguments

- `q::Float64`: **Rotational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Am(q::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    -R * T * log(q)
end


@doc raw"""
    Rotational.Gm(q::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Rotational Molar Gibbs Free Energy**.

# Detail

```math
G_{m, rot} = - N_A k_B T \ln q_{rot} = - R T \ln q_{rot}
```

# Arguments

- `q::Float64`: **Rotational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Gm(q::Float64, T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    -R * T * log(q)
end


@doc raw"""
    Rotational.Um(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Rotational Molar Internal Thermal Energy**

# Detail

```math
U_{m, rot} = \frac{3}{2} N_A k_B T = \frac{3}{2} R T
```

# Arguments

- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Um(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    3/2 * R * T
end

@doc raw"""
    Rotational.Hm(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Rotational Molar Enthalpy**

# Detail

```math
H_{m, rot} = \frac{3}{2} N_A k_B T = \frac{3}{2} R T
```

# Arguments

- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Hm(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    3/2 * R * T
end


@doc raw"""
    Rotational.Sm(q::Float64)::typeof(1.0u"J/mol/K")

Compute **Rotational Entropy**.

# Detail

```math
S_{m, rot} = N_A k_B (\frac{3}{2} + \ln q_{rot}) = R (\frac{3}{2} + \ln q_{rot})
```

# Arguments

- `q::Float64`: Electronic Partition Function.
"""
function Sm(q::Float64)::typeof(1.0u"J/mol/K")
    R * (3/2 + log(q))
end


@doc raw"""
    Rotational.CVm()::typeof(1.0u"J/mol/K")

Compute **Rotational Molar Volume Heat Capacity**.

# Detail

```math
C_{Vm, rot} = \frac{3}{2} N_A k_B = \frac{3}{2} R
```
"""
function CVm()::typeof(1.0u"J/mol/K")
    3/2 * R
end

@doc raw"""
    Rotational.Cpm()::typeof(1.0u"J/mol/K")

Compute **Rotational Molar Pressure Heat Capacity**.

# Detail

```math
C_{pm, rot} = \frac{3}{2} N_A k_B = \frac{3}{2} R
```
"""
function Cpm()::typeof(1.0u"J/mol/K")
    3/2 * R
end

end
