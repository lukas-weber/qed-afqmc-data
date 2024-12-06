
function digon(w, h)
    R = (h^2 + w^2) / 4w
    α = asin(h / 2R)
    digon = BezierPath([
        EllipticalArc(Point(-R + w / 2, 0), R, R, 0, -α, α),
        EllipticalArc(Point(R - w / 2, 0), R, R, 0, π - α, π + α),
        ClosePath(),
    ])
end

const molcolors =
    Dict("CH₄" => "#51585d", "NH₃" => "#588cc2", "H₂O" => "#ca5662", "FH" => "#9db974")

molnamemap(name) =
    Dict("CH₄" => L"CH$_4$", "H₂O" => L"H$_2$O", "NH₃" => L"NH$_3$", "FH" => L"HF$$")[name]

function basisnamemap(name)
    if name == "cbs"
        return "aug-cc-pV∞Z"
    end
    return replace(name, r"v.z" => uppercase)
end


function load_data(filename)
    data = JSON.parsefile("../data/$filename.json")

    df = DataFrame(Tables.dictrowtable(data))

    function combine_errors(mean, err)
        if isnothing(mean)
            return missing
        end

        if iszero(err)
            return mean
        end

        return measurement(mean, err)
    end

    for column in names(df)
        if endswith(column, "Error")
            prefix = column[1:end-length("Error")]
            select!(df, Not(column), [prefix, column] => ByRow(combine_errors) => prefix)
        end
    end

    return df
end

function unstack_obs(df, obsname)
    return unstack(
        df,
        [:ε, :molname, :basis, :quadrupole_enabled, :gauge],
        :method,
        obsname,
        renamecols = m -> "$obsname$(m[5:end])",
        allowmissing = true,
    )
end

function extrapolate_basis(
    df::DataFrame;
    observables = [:CorrEnergyQMC, :CorrEnergyCCSD, :CorrEnergyCCSDT],
    params = [:molname, :ε],
)
    dfinf = combine(
        groupby(df, params),
        [] => Returns("aug-cc-pv∞z") => :basis,
        ([:basis, obs] => extrapolate_basis => obs for obs in observables)...,
    )
    return dfinf
end

function extrapolate_basis(basis::AbstractVector, observable::AbstractVector)
    xdict = Dict('d' => 2, 't' => 3, 'q' => 4, '5' => 5)

    @. model(x, p) = p[1] + p[2] * x^(-3)
    x = getindex.(basis, 10)

    obst = observable[findfirst(x .== 't')]
    obsq = observable[findfirst(x .== 'q')]

    return (obsq * 4^3 - obst * 3^3) / (4^3 - 3^3)
end
