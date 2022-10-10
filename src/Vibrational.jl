module Vibrational

using Unitful
import PhysicalConstants.CODATA2018: R, k_B, N_A, h, c_0

@doc raw"""
    Vibrational.calcΘvib(ν̃)::typeof(0.0u"K")

Compute **Characteristic Vibrational Temperature**.

# Detail

```math
\nu = c_0 \tilde\nu \\
```

where ``c_0`` is the **Speed of Light**, ``\tilde\nu`` is wavenumber.

```math
\Theta_{rot} = \frac{h \nu}{k_B}
```

where ``k_B`` is **Boltzmann Constant**, ``h`` is **Plank Constant**, ``\nu`` is frequency.

# Arguments

- `ν̃`: Wavenumber, must in the unit of `u"1/cm"`
"""
function calcΘvib(ν̃)::typeof(0.0u"K")
    @assert unit(ν̃) == unit(1.0u"1/cm") "ν̃ should in the unit of 1/cm!"
    h * c_0 * ν̃ / k_B |> upreferred
end


@doc raw"""
    Vibrational.q(ν̃s; T::Unitful.Temperature=298.15u"K";)::Float64

Compute **Vibrational Partition Function**.

# Detail

For each vibration mode, here ``Θ_{vib, i}`` is characteristic vibration temperature for vibration mode ``i``

```math
q_{vib, i} = \frac{\exp [\frac{- \Theta_{vib, i}}{2T}]}{1 - \exp [\frac{- \Theta_{vib, i}}{T}]}
```

Overall Vibrational Partition Function

```math
q_{vib} = \Pi q_{vib, i}
```

Choose the first vibration energy level to be the zero of energy

```math
q_{vib, i} = \frac{1}{1 - \exp [\frac{- \Theta_{vib, i}}{T}]}
```

Overall Vibrational Partition Function

```math
q_{vib} = \prod q_{vib, i}
```

# Arguments

- `ν̃s`: Vector of wavenumbers, all elements must in the unit of `u"1/cm"`
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""

function q(ν̃s; T::Unitful.Temperature=298.15u"K")::Float64
    @assert unit(ν̃s[1]) == unit(1.0u"1/cm") "list of ν̃ should in the unit of 1/cm!"
    T = T |> u"K"
    map(ν̃s) do ν̃
        Θ_vib = calcΘvib(ν̃)
        1 / (1 - exp(-Θ_vib / T))
    end |> prod
end


@doc raw"""
    Vibrational.Am(q_rot::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Vibrational Molar Helmholtz Free Energy**.

# Detail

```math
A_{m, rot} = - N_A k_B T \ln q_{vib} = - R T \ln q_{vib}
```

# Arguments

- `q::Float64`: **Rotational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Am(q::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    -R * T * log(q)
end


@doc raw"""
    Vibrational.Gm(q::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Vibrational Molar Gibbs Free Energy**.

# Detail

```math
G_{m, vib} = - N_A k_B T \ln q_{vib} = - R T \ln q_{vib}
```

# Arguments

- `q::Float64`: **Rotational Partition Function**
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Gm(q::Float64; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    T = T |> u"K"
    Am(q, T=T)
end


@doc raw"""
    Vibrational.Um(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Vibrational Molar Internal Thermal Energy**

# Detail

```math
U_{m, vib} = R \sum \frac{\Theta_{vib, i}}{\exp \left[\frac{\Theta_{vib, i}}{T}\right] - 1}
```

# Arguments

- `ν̃s`: Vector of wavenumbers, all elements must in the unit of `u"1/cm"`
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Um(ν̃s; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    @assert unit(ν̃s[1]) == unit(1.0u"1/cm") "list of ν̃ should in the unit of 1/cm!"
    T = T |> u"K"
    map(ν̃s) do ν̃
        Θ_vib = calcΘvib(ν̃)
        R * Θ_vib * (1/2 + 1 / (exp(Θ_vib / T) - 1))
    end |> sum
end


@doc raw"""
    Vibrational.Hm(T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")

Compute **Vibrational Molar Enthalpy**

# Detail

```math
H_{m, vib} = R \sum \frac{\Theta_{vib, i}}{\exp \left[\frac{\Theta_{vib, i}}{T}\right] - 1}
```

# Arguments

- `ν̃s`: Vector of wavenumbers, all elements must in the unit of `u"1/cm"`
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Hm(ν̃s; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol")
    @assert unit(ν̃s[1]) == unit(1.0u"1/cm") "list of ν̃ should in the unit of 1/cm!"
    T = T |> u"K"
    Um(ν̃s, T=T)
end


@doc raw"""
    Vibrational.Sm(q::Float64)::typeof(1.0u"J/mol/K")

Compute **Vibrational Entropy**.

# Detail

```math
S_{m, vib} = R \sum \left\{ - \ln \left[ 1 - \exp [- \frac{\Theta_{vib, i}}{T}] \right] + \frac{\frac{\Theta_{vib,i}}{T}}{\exp [\frac{\Theta_{vib, i}}{T}] - 1} \right\}
```

# Arguments

- `ν̃s`: Vector of wavenumbers, all elements must in the unit of `u"1/cm"`
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Sm(ν̃s; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol/K")
    @assert unit(ν̃s[1]) == unit(1.0u"1/cm") "list of ν̃ should in the unit of 1/cm!" 
    T = T |> u"K"
    map(ν̃s) do ν̃
        Θ_vib = calcΘvib(ν̃)
        R * (Θ_vib / T / (exp(Θ_vib / T) - 1) - log(1 - exp(-Θ_vib / T)))
    end |> sum
end


@doc raw"""
    Vibrational.CVm(ν̃s; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol/K")

Compute **Vibrational Molar Volume Heat Capacity**.

# Detail

```math
C_{Vm, vib} = R \sum \frac{{(\frac{\Theta_{vib, i}}{T})}^2 \exp [ \frac{\Theta_{vib, i}}{T} ]}{{(\exp [ \frac{\Theta_{vib, i}}{T} ] - 1)}^2}
```

# Arguments

- `ν̃s`: Vector of wavenumbers, all elements must in the unit of `u"1/cm"`
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function CVm(ν̃s; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol/K")
    @assert unit(ν̃s[1]) == unit(1.0u"1/cm") "list of ν̃ should in the unit of 1/cm!"
    T = T |> u"K"
    map(ν̃s) do ν̃
        Θ_vib = calcΘvib(ν̃)
        R  * (Θ_vib / T)^2 * exp(Θ_vib / T) / (exp(Θ_vib / T) - 1)^2
    end |> sum
end


@doc raw"""
    Vibrational.Cpm(ν̃s; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol/K")

Compute **Vibrational Molar Pressure Heat Capacity**.

# Detail

```math
C_{pm, vib} = R \sum \frac{{(\frac{\Theta_{vib, i}}{T})}^2 \exp [ \frac{\Theta_{vib, i}}{T} ]}{{(\exp [ \frac{\Theta_{vib, i}}{T} ] - 1)}^2}
```

# Arguments

- `ν̃s`: Vector of wavenumbers, all elements must in the unit of `u"1/cm"`
- `T::Unitful.Temperature`: Temperature, with default value of `298.15u"K"`
"""
function Cpm(ν̃s; T::Unitful.Temperature=298.15u"K")::typeof(1.0u"J/mol/K")
    @assert unit(ν̃s[1]) == unit(1.0u"1/cm") "list of ν̃ should in the unit of 1/cm!"
    T = T |> u"K"
    CVm(ν̃s, T=T)
end

end
