using CairoMakie
using FileIO
using LaTeXStrings
using DataFrames
using Statistics
using Measurements
using ColorSchemes
using JSON

include("makie.jl")
include("utils.jl")

function fig_gauges(df)
    df = filter(r -> r.molname == "H₂" && !r.quadrupole_enabled, df)
    dfref = select(
        filter(r -> r.method == "QED-HF" && r.gauge == "Coulomb", df),
        :ε,
        :basis,
        :Energy => :EnergyHF,
    )
    df = filter(r -> r.method in ("QED-QMC", "QED-ED"), df)
    df = leftjoin(unstack_obs(df, :Energy), dfref, on = [:ε, :basis])
    disallowmissing!(df, [:EnergyED, :EnergyHF])
    df.EnergyQMC .= float.(df.EnergyQMC)

    dftmp = copy(df)

    dftmp.CorrEnergyQMC = dftmp.EnergyQMC - dftmp.EnergyHF
    dftmp.CorrEnergyED = dftmp.EnergyED - dftmp.EnergyHF
    dfinf = extrapolate_basis(
        filter(:basis => startswith("aug-"), dftmp);
        observables = [:CorrEnergyQMC, :CorrEnergyED],
        params = [:molname, :ε, :gauge, :quadrupole_enabled],
    )
    dfinf = select(
        leftjoin(
            dfinf,
            select(filter(:basis => ==("aug-cc-pvqz"), dfref), Not(:basis)),
            on = :ε,
        ),
        Not(:CorrEnergyQMC, :CorrEnergyED),
        [:CorrEnergyQMC, :EnergyHF] => (+) => :EnergyQMC,
        [:CorrEnergyED, :EnergyHF] => (+) => :EnergyED,
    )

    df = vcat(df, dfinf)


    E0 = filter(:ε => ==(0.0), dfinf).EnergyED[1]

    df = filter(:ε => >(0.0), df)
    disallowmissing!(df, :EnergyQMC)

    @show dfinf.EnergyED[4] - dfinf.EnergyED[3]
    @show dfinf.EnergyQMC[4] - dfinf.EnergyQMC[3]

    fig = Figure(size = (380, 250))

    g = 0.3
    basismap = Dict(
        "cc-pvdz" => 1 - g,
        "cc-pvtz" => 3 - g,
        "cc-pvqz" => 5 - g,
        "cc-pv5z" => 7 - g,
        "aug-cc-pvdz" => 1 + g,
        "aug-cc-pvtz" => 3 + g,
        "aug-cc-pvqz" => 5 + g,
        "aug-cc-pv5z" => 7 + g,
        "aug-cc-pv∞z" => 8 + g,
    )
    basismap = filter(x -> endswith(x[1], "∞z") || in(x[1], df.basis), basismap)
    bases = collect(keys(basismap))


    dfdipole = filter(:gauge => ==("Dipole"), df)
    dfcoulomb = filter(:gauge => ==("Coulomb"), df)


    ax = Axis(
        fig[1, 1],
        ylabel = L"$E$ (Ha)",
        xticks = (getindex.(Ref(basismap), bases), basisnamemap.(bases)),
        xticklabelrotation = -π / 4,
        ygridvisible = false,
        xgridvisible = false,
    )

    inset_ax = Axis(
        fig[1, 1],
        width = Relative(0.5),
        height = Relative(0.5),
        halign = 0.3,
        valign = 1.1,
        aspect = DataAspect(),
    )

    image!(inset_ax, rotr90(load("../raw/h2_cavity.png")))
    hidedecorations!(inset_ax)
    hidespines!(inset_ax)
    Acolor = "#88ccff"
    x = 1160
    arrows!(
        inset_ax,
        [x],
        [410],
        [0],
        [200],
        space = :relative,
        linewidth = 2,
        color = Acolor,
    )
    text!(inset_ax, x + 30, 452, text = "𝛜", color = Acolor)

    colors = Makie.wong_colors()

    basisx(basis) = getindex.(Ref(basismap), basis)

    shift = 0.08
    edplots = [
        scatter!(
            ax,
            basisx.(dfcoulomb.basis) .- shift,
            dfcoulomb.EnergyED,
            color = :white,
            strokecolor = colors[1],
            strokewidth = 1,
            label = "ED Coulomb",
        ),
        scatter!(
            ax,
            basisx.(dfdipole.basis) .+ shift,
            dfdipole.EnergyED,
            color = :white,
            strokecolor = colors[2],
            marker = :rect,
            strokewidth = 1,
            label = "ED Dipole",
        ),
    ]
    xlims!(ax, 0.4, 10.1)
    qmcplots = [
        measplot!(
            ax,
            basisx.(dfcoulomb.basis) .- shift,
            dfcoulomb.EnergyQMC,
            errorsinfront = true,
            label = "QMC Coulomb",
        ),
        measplot!(
            ax,
            basisx.(dfdipole.basis) .+ shift,
            dfdipole.EnergyQMC,
            color = colors[2],
            marker = :rect,
            label = "QMC dipole",
            errorsinfront = true,
        ),
    ]

    @show mean(dfcoulomb.EnergyED - dfcoulomb.EnergyQMC)
    @show mean(dfdipole.EnergyED - dfdipole.EnergyQMC)

    emptyelem = [PolyElement(color = :transparent)]
    axislegend(
        ax,
        [edplots, qmcplots, [emptyelem, emptyelem]],
        [repeat([""], 4), repeat([""], 4), ["Coulomb", "Dipole"]],
        ["FCI", "QMC", ""],
        position = :rt,
        orientation = :horizontal,
        nbanks = 2,
        groupgap = 5,
        titlegap = 3,
        titlesize = 12,
        patchlabelgap = 0,
        padding = 0,
        rowgap = 0,
        patchsize = (1, 14),
        margin = (0, 3, 0, 0),
    )

    text!(
        ax,
        7.9,
        -1.163,
        text = L"H$_2$, $|𝛜| = 0.1\,\text{a.u.}$",
        align = (:center, :center),
    )

    if E0 !== nothing
        Ω = 0.3
        hlines!(ax, E0, color = :black)
        text!(ax, 4.0, E0 + 0.0009; text = "no cavity", align = (:center, :center))
    end
    return fig
end

function fig_mols(df)
    df = filter(
        r ->
            r.molname in ["CH₄", "NH₃", "H₂O", "FH"] &&
                !r.quadrupole_enabled &&
                (r.method, r.gauge) in
                [("QED-QMC", "Coulomb"), ("QED-HF", "Coulomb"), ("QED-CC", "Dipole")] &&
                startswith(r.basis, "aug"),
        df,
    )

    df.gauge .= missing
    df = unstack_obs(df, :Energy)

    ccsdt_shift = DataFrame(JSON.parsefile("../data/ccsdt_shifts.json"))

    df = leftjoin(df, ccsdt_shift, on = [:molname, :basis])

    fig = Figure(size = (380, 450), figure_padding = 7)
    ax0 = Axis(fig[1, 1:2], aspect = DataAspect())

    image!(ax0, rotr90(load("../raw/mol_cavity.png")))
    hidedecorations!(ax0)
    hidespines!(ax0)

    Acolor = "#88ccff"
    x = 1650
    arrows!(ax0, [x], [430], [0], [200], space = :relative, linewidth = 2, color = Acolor)
    text!(ax0, x + 19, 498, text = "𝛜", color = Acolor)

    ax1 = Axis(
        fig[2, 1],
        xlabel = L"|𝛜|~\mathrm{(a.u.)}",
        ylabel = L"E_\text{c}\,\mathrm{(Ha)}",
    )

    df.CorrEnergyQMC = df.EnergyQMC - df.EnergyHF
    df.CorrEnergyCCSD = df.EnergyCC - df.EnergyHF
    df.CorrEnergyCCSDT = df.EnergyCC - df.EnergyHF + df.ccsdt

    df = vcat(df, extrapolate_basis(df), cols = :union)

    plot_mols_bases(ax1, df)
    ax2 = Axis(fig[2, 2], xlabel = L"|𝛜|~\mathrm{(a.u.)}", ylabel = L"E-E_\mathrm{HF}")
    text!(ax0, 0.025, 0.98, space = :relative, text = "(a)", align = (:left, :top))

    plot_mols(ax2, df)
    linkyaxes!(ax1, ax2)
    ylims!(ax2, -0.44, -0.161)
    text!(ax1, 0.025, 0.98, space = :relative, text = L"(b) CH$_4$", align = (:left, :top))
    text!(ax2, 0.025, 0.98, space = :relative, text = "(c)", align = (:left, :top))
    hideydecorations!(ax2, ticks = false)
    rowgap!(fig.layout, 1, -20)
    return fig
end


function plot_mols_bases(ax, df)

    qmcplots = []
    ccsdplots = []
    bases = []

    dfc = filter(r -> r.molname == "CH₄", df)

    marker_cycle = [digon(0.45, 0.8), :utriangle, :rect, :circle]

    colorscheme = [colorschemes[:BuPu][[0.4, 0.5, 0.6]]..., molcolors["CH₄"]]
    for (i, dfbasis) in enumerate(groupby(dfc, :basis))
        push!(
            qmcplots,
            scatter!(
                ax,
                dfbasis.ε,
                Measurements.value.(dfbasis.CorrEnergyQMC),
                color = colorscheme[i],
                marker = marker_cycle[i],
            ),
        )
        push!(
            ccsdplots,
            scatter!(
                ax,
                dfbasis.ε,
                (dfbasis.CorrEnergyCCSD),
                color = :white,
                strokecolor = colorscheme[i],
                strokewidth = 1,
                marker = marker_cycle[i],
            ),
        )
        push!(bases, basisnamemap(dfbasis.basis[1]))
    end


    dfinf = filter(:basis => ==("aug-cc-pv∞z"), dfc)
    push!(
        ccsdplots,
        scatter!(
            ax,
            dfinf.ε,
            dfinf.CorrEnergyCCSDT,
            label = "aug-cc-pv∞z",
            color = :white,
            strokecolor = molcolors["CH₄"],
            strokewidth = 1,
            markersize = 6,
            marker = :circle,
        ),
    )
    push!(bases, L"+(T$_0$)")
    push!(qmcplots, PolyElement(color = :white))

    axislegend(
        ax,
        [ccsdplots, qmcplots],
        [repeat([""], 5), bases],
        ["CCSD", "QMC"],
        orientation = :horizontal,
        nbanks = 5,
        position = :lb,
        titlesize = 12,
        titlehalign = :left,
        groupgap = 3,
        gridshalign = :right,
        padding = 0,
        patchlabelgap = 0,
        titlegap = 3,
        rowgap = 0,
    )
end

function plot_mols(ax, df)
    qmcplots = []
    ccsdplots = []
    mols = LaTeXString[]

    dfinf = filter(:basis => ==("aug-cc-pv∞z"), df)

    @show mean(abs.(dfinf.CorrEnergyQMC - dfinf.CorrEnergyCCSD))
    @show mean(abs.(dfinf.CorrEnergyQMC - dfinf.CorrEnergyCCSDT))

    for (i, dfmol) in enumerate(groupby(dfinf, :molname))
        push!(
            qmcplots,
            scatter!(
                ax,
                dfmol.ε,
                Measurements.value.(dfmol.CorrEnergyQMC),
                color = molcolors[dfmol.molname[1]],
                marker = marker_cycle[i],
            ),
        )
        push!(
            ccsdplots,
            scatter!(
                ax,
                dfmol.ε,
                dfmol.CorrEnergyCCSDT,
                color = :white,
                strokecolor = molcolors[dfmol.molname[1]],
                strokewidth = 1,
                marker = marker_cycle[i],
                markersize = 6,
            ),
        )
        push!(mols, molnamemap(dfmol.molname[1]))
    end
    axislegend(
        ax,
        [ccsdplots, qmcplots],
        [repeat([""], 4), mols],
        ["CCSD(T₀)", "QMC"],
        orientation = :horizontal,
        nbanks = 4,
        groupgap = 3,
        titlegap = 3,
        titlesize = 12,
        patchlabelgap = 0,
        gridshalign = :right,
        padding = 0,
        rowgap = 0,
        valign = :top,
        margin = (0, 3, 0, 0),
    )

end

function fig_fixed_basis(df)
    df = filter(
        r ->
            !r.quadrupole_enabled &&
                r.molname != "H₂" &&
                (r.gauge, r.method) in
                [("Coulomb", "QED-HF"), ("Dipole", "QED-QMC"), ("Dipole", "QED-CC")] &&
                r.basis == "aug-cc-pvdz",
        df,
    )
    df.gauge .= missing
    df = unstack_obs(df, :Energy)

    ccsdt_shift = DataFrame(JSON.parsefile("../data/ccsdt_shifts.json"))
    df = leftjoin(df, ccsdt_shift, on = [:molname, :basis])

    fig = Figure(size = (380, 250))

    ax = Axis(fig[1, 1], xlabel = L"$|𝛜|$ (a.u.)", ylabel = L"$E_\text{c}$ (Ha)")

    qmcplots = []
    ccsdplots = []
    molnames = []
    for (i, dfmol) in Iterators.reverse(enumerate(groupby(df, :molname)))
        push!(
            qmcplots,
            scatter!(
                ax,
                dfmol.ε,
                Measurements.value.(dfmol.EnergyQMC) - dfmol.EnergyHF,
                color = molcolors[dfmol.molname[1]],
                label = molnamemap(dfmol.molname[1]),
                marker = marker_cycle[i],
                strokewidth = 1,
            ),
        )
        push!(
            ccsdplots,
            scatter!(
                ax,
                dfmol.ε,
                dfmol.EnergyCC + dfmol.ccsdt - dfmol.EnergyHF,
                color = :transparent,
                strokecolor = molcolors[dfmol.molname[1]],
                strokewidth = 1,
                marker = marker_cycle[i],
                markersize = 6,
            ),
        )
        push!(molnames, molnamemap(dfmol.molname[1]))
    end
    @show mean(abs.(df.EnergyQMC - df.EnergyCC - df.ccsdt))
    text!(ax, 0.0, -0.283, text = "aug-cc-pVDZ\nDipole gauge", align = (:left, :center))
    axislegend(
        ax,
        [ccsdplots, qmcplots],
        [repeat([""], 4), molnames],
        ["CCSD(T₀)", "QMC"],
        orientation = :horizontal,
        nbanks = 4,
        groupgap = 3,
        titlegap = 3,
        titlesize = 12,
        patchlabelgap = 0,
        gridshalign = :right,
        padding = 0,
        rowgap = 0,
        valign = :top,
        margin = (0, 3, 0, 0),
    )
    return fig
end


function fig_photon_number(df)
    df = filter(
        r ->
            r.molname == "FH" &&
                !r.quadrupole_enabled &&
                r.method in ["QED-QMC", "QED-CC"] &&
                startswith(r.basis, "aug") &&
                r.ε > 0,
        df,
    )

    df = unstack_obs(df, :PhotonNumber)

    fig = Figure(size = (380, 250))

    ax = Axis(fig[1, 1], xlabel = L"$|𝛜|$ (a.u.)", ylabel = L"n_\mathrm{ph}")
    colorscheme = [colorschemes[:BuPu][[0.4, 0.5, 0.6]]..., molcolors["FH"]]
    inset_ax = Axis(
        fig[1, 1],
        width = Relative(0.25),
        height = Relative(0.3),
        halign = 0.45,
        valign = 0.9,
        ylabel = L"n_\mathrm{ph}",
        xticks = ([2, 3, 4], ["D", "T", "Q"]),
        yticks = [0.025, 0.030],
        xlabel = "aug-cc-pV•Z",
    )

    text!(
        ax,
        0.03,
        0.98,
        space = :relative,
        text = "aug-cc-pVQZ\nHF",
        align = (:left, :top),
    )

    plot_photon_number_convergence(inset_ax, df)

    dfqz = filter(:basis => ==("aug-cc-pvqz"), df)
    dfcoulomb = dropmissing(filter(:gauge => ==("Coulomb"), dfqz), :PhotonNumberQMC)
    dfdipole = dropmissing(filter(:gauge => ==("Dipole"), dfqz))
    scatter!(
        ax,
        dfdipole.ε,
        dfdipole.PhotonNumberCC,
        color = :white,
        strokecolor = molcolors["FH"],
        marker = :rect,
        strokewidth = 1,
        label = "CCSD",
    )
    measplot!(
        ax,
        dfcoulomb.ε,
        dfcoulomb.PhotonNumberQMC,
        marker = :rect,
        color = molcolors["FH"],
        label = "QMC Coulomb",
    )
    measplot!(
        ax,
        dfdipole.ε,
        dfdipole.PhotonNumberQMC,
        marker = :rect,
        color = molcolors["FH"],
        label = "QMC Dipole",
        strokecolor = :black,
        strokewidth = 1,
    )

    @show mean(
        (dfdipole.PhotonNumberCC - dfcoulomb.PhotonNumberQMC) ./ dfcoulomb.PhotonNumberQMC,
    )
    @show mean(
        (dfdipole.PhotonNumberQMC - dfcoulomb.PhotonNumberQMC) ./ dfcoulomb.PhotonNumberQMC,
    )

    axislegend(ax, position = :rb)
    return fig
end

function plot_photon_number_convergence(ax, df)
    dfcoulomb =
        dropmissing(filter(r -> r.ε == 0.1 && r.gauge == "Coulomb", df), :PhotonNumberQMC)
    dfdipole = dropmissing(filter(r -> r.ε == 0.1 && r.gauge == "Dipole", df))
    marker_cycle = [digon(0.45, 0.8), :utriangle, :rect]
    colors = molcolors["FH"]

    xs = Dict('d' => 2, 't' => 3, 'q' => 4)

    scatter!(
        ax,
        getindex.(Ref(xs), getindex.(dfdipole.basis, 10)),
        dfdipole.PhotonNumberCC,
        color = :white,
        strokecolor = colors,
        strokewidth = 1,
        marker = marker_cycle,
    )
    measplot!(
        ax,
        getindex.(Ref(xs), getindex.(dfcoulomb.basis, 10)),
        dfcoulomb.PhotonNumberQMC,
        color = colors,
        marker = marker_cycle,
    )
    measplot!(
        ax,
        getindex.(Ref(xs), getindex.(dfdipole.basis, 10)),
        dfdipole.PhotonNumberQMC,
        color = colors,
        strokewidth = 1,
        marker = marker_cycle,
    )

    text!(
        ax,
        0.5,
        0.86,
        text = L"$|𝛜|=0.1$ a.u.",
        space = :relative,
        align = (:center, :center),
    )
    ylims!(ax, 0.023, 0.0350)
end

function fig_quadrupole(df)
    dfccsd = filter(:method => ==("QED-CC"), df)

    df = filter(:ε => ==(0.1), dfccsd)
    df0 = filter(:ε => ==(0.0), dfccsd)

    df = leftjoin(
        df,
        select(df0, :molname, :basis, :quadrupole_enabled, :Energy => :Energy0),
        on = [:molname, :basis, :quadrupole_enabled],
    )

    fig = Figure(size = (380, 250), figure_padding = 3)

    xs = Dict("aug-cc-pvdz" => 2, "aug-cc-pvtz" => 3, "aug-cc-pvqz" => 4, "cbs" => 5)
    ax = Axis(
        fig[1, 1],
        ylabel = L"E_{|𝛜|=0.1}-E_{𝛜=0}",
        xlabel = "aug-cc-pV•Z",
        xticks = ([2, 3, 4, 5], ["D", "T", "Q", "∞"]),
    )


    for dfmol in groupby(df, :molname)
        dfquad = subset(dfmol, :quadrupole_enabled)
        dfnoquad = subset(dfmol, :quadrupole_enabled => ByRow(!))
        dfcbs = filter(:basis => ==("cbs"), dfnoquad)

        text!(
            ax,
            5,
            dfcbs.Energy[1] - dfcbs.Energy0[1] + 0.0002;
            text = molnamemap(dfcbs.molname[1]),
            align = (:center, :bottom),
        )
        color = molcolors[dfmol.molname[1]]
        hlines!(ax, dfcbs.Energy - dfcbs.Energy0, color = color)
        scatter!(
            ax,
            getindex.(Ref(xs), dfquad.basis),
            dfquad.Energy - dfquad.Energy0,
            color = color,
        )
        scatter!(
            ax,
            getindex.(Ref(xs), dfnoquad.basis),
            dfnoquad.Energy - dfnoquad.Energy0,
            color = :white,
            strokecolor = color,
            strokewidth = 1,
        )
    end

    elements = [
        MarkerElement(color = :black, marker = :circle),
        MarkerElement(
            color = :white,
            marker = :circle,
            strokewidth = 1,
            strokecolor = :black,
        ),
    ]
    axislegend(
        ax,
        elements,
        [L"P\textbf{d}^2P", L"(P\textbf{d}P)^2"],
        position = (0.01, 0.01),
        padding = 0,
        orientation = :horizontal,
    )

    ylims!(ax, 0.0095, 0.0195)
    fig
end



function fig_nonint(dfnonint)
    fig = Figure(size = (380, 250), figure_padding = 3)
    ax = Axis(fig[1, 1], ylabel = L"$E_\mathrm{c}$ (Ha)", xlabel = L"$|𝛜|$ (a.u.)")

    df = filter(
        r ->
            (r.method, r.gauge) in
            [("QED-CC", "Dipole"), ("QED-QMC", "Dipole"), ("QED-HF", "Coulomb")] &&
                r.basis == "aug-cc-pvqz",
        dfnonint,
    )
    df.gauge .= missing

    df = unstack_obs(df, :Energy)
    disallowmissing(df, [:EnergyHF, :EnergyQMC, :EnergyCC])


    measplot!(
        ax,
        df.ε,
        df.EnergyQMC - df.EnergyHF,
        marker = :rect,
        color = molcolors["FH"],
        label = "QMC",
    )
    scatter!(
        ax,
        df.ε,
        df.EnergyCC - df.EnergyHF,
        marker = :xcross,
        markersize = 6,
        color = :black,
        label = "CCSD-2",
    )

    axislegend(ax, padding = 0)
    text!(
        ax,
        -0.001,
        0.0035,
        text = "HF, aug-cc-pVQZ\nDipole gauge, no el-el repulsion",
        align = (:left, :bottom),
        justification = :left,
        fontsize = 12,
    )
    inset_ax = Axis(
        fig[1, 1],
        width = Relative(0.25),
        height = Relative(0.3),
        halign = 0.35,
        valign = 0.3,
        ylabel = L"$E_\mathrm{QMC}-E_\mathrm{CCSD‑2}$ (μHa)",
        xlabel = "|𝛜| (a.u.)",
    )

    measplot!(
        inset_ax,
        df.ε,
        1e6 * (df.EnergyQMC - df.EnergyCC),
        color = molcolors["FH"],
        marker = :rect,
    )
    xlims!(ax, -0.0035, 0.105)
    ylims!(ax, -0.068, 0.016)

    fig
end

function fig_ccsd1(df)
    df = filter(r -> r.method in ["QED-CC", "QED-CC-1"] && !r.quadrupole_enabled, df)
    df = unstack_obs(df, :Energy)

    fig = Figure(size = (380, 250))
    ax = Axis(
        fig[1, 1],
        ylabel = L"$E_\mathrm{CCSD‑2}-E_\mathrm{CCSD‑1}$ (mHa)",
        xlabel = L"$|𝛜|$ (a.u.)",
    )

    marker_cycle = Dict(
        "aug-cc-pvdz" => digon(0.45, 0.8),
        "aug-cc-pvtz" => :utriangle,
        "aug-cc-pvqz" => :rect,
        "cbs" => :circle,
    )

    scatter!(
        ax,
        df.ε,
        1e3 * (df.EnergyCC - df.var"EnergyCC-1"),
        strokecolor = getindex.(Ref(molcolors), df.molname),
        strokewidth = 1,
        color = :white,
        marker = getindex.(Ref(marker_cycle), df.basis),
    )

    molnames = unique(df.molname)
    bases = unique(df.basis)

    molelems = [PolyElement(; color = molcolors[molname]) for molname in molnames]

    basiselems = [
        MarkerElement(; strokewidth = 1, marker = marker_cycle[basis], color = :white)
        for basis in bases
    ]

    axislegend(
        ax,
        vcat(molelems, basiselems),
        vcat(molnamemap.(molnames), basisnamemap.(bases)),
        nbanks = 4,
        orientation = :horizontal,
        position = :lb,
    )

    fig
end

function plots()
    pt_per_unit = 0.625

    df = load_data("molecules")
    df_nonint = load_data("molecules_no_coulomb")

    with_theme(qed_theme()) do
        save("../figs/mols.pdf", fig_mols(df); pt_per_unit)
        save("../figs/nph.pdf", fig_photon_number(df); pt_per_unit)
        save("../figs/quadrupole.pdf", fig_quadrupole(df); pt_per_unit)
        save("../figs/h2_gauges.pdf", fig_gauges(df); pt_per_unit)
        save("../figs/dipole_fixed_basis.pdf", fig_fixed_basis(df); pt_per_unit)
        save("../figs/hf_nonint.pdf", fig_nonint(df_nonint); pt_per_unit)
        save("../figs/ccsd1.pdf", fig_ccsd1(df); pt_per_unit)
    end
end

function (@main)(args)
    plots()
    return 0
end
