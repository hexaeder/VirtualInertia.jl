using Graphs
using Printf
using MetaGraphs
using GraphMakie
using Makie
using GLMakie
using Makie.Colors
using Makie.ColorSchemes
using NetworkDynamics
using GraphMakie.NetworkLayout: spring

export inspect_solution

function inspect_solution(sol, network)
    fig = Figure(resolution = (1200, 1200))
    ####
    #### selector grid to control plot
    ####
    symgrid = fig[1,1] = GridLayout(tellwidth=false, tellheight=true)
    buttons = symgrid[1,1:5] = [Button(fig, label="_ω"),
                                Button(fig, label="_u_arg"),
                                Button(fig, label="_u_mag"),
                                Button(fig, label="_P"),
                                Button(fig, label="_rocof")]
    symbox = symgrid[1,6] = Textbox(fig, width=150)
    reltoggle = symgrid[1,7] = Toggle(fig)
    symgrid[1,8] = Label(fig, "relativ to u0")

    rel_to_u0 = reltoggle.active

    set_nsym! = function(sym)
        symbox.displayed_string[] = sym
        symbox.stored_string[] = sym
    end

    for i in 1:5
        on(buttons[i].clicks) do n
            sym = buttons[i].label[]
            set_nsym!(sym)
        end
    end

    nstatesym = Observable(:_ω)
    on(symbox.stored_string) do s
        nstatesym[] = Symbol(s)
    end

    set_nsym!("_ω")

    #####
    ##### Selectors for sel_nodes and sel_edges
    #####
    selgrid = fig[2,1] = GridLayout(tellwidth=false, tellheight=true)

    sel_nodes = Observable(Set())
    sel_edges = Observable(Set())
    nselectors = selgrid[1,1:3] = [Label(fig, "Selected Nodes:"),
                                       Textbox(fig; width=300, placeholder="insert comma separated indices"),
                                       Button(fig, label="open node plot")]
    on(nselectors[2].stored_string) do s
        if occursin(r"^\s*$", s)
            sel_nodes[] = Set{Int}()
        else
            try
                parts = split(s, ',')
                ints = parse.(Int, parts)
                sel_nodes[] = Set(ints)
            catch e
                @warn "Parsing error of $s"
            end
        end
    end
    on(sel_nodes) do set
        str = prod(string.(sort(collect(set))) .* ", ")[begin:end-2]
        nselectors[2].displayed_string[] = str *" "
    end

    #####
    ##### Graphplot
    #####
    gpgrid = fig[3,1] = GridLayout()

    ##### Bottom elements
    bottom_sg = SliderGrid(fig[4,1],
                           (label = "Node color scaling", range=Base.range(-10, 0, length=100), startvalue=0, format=x->@sprintf("%.2E", 10^x)),
                           (label = "Time", range=Base.range(sol.t[begin], sol.t[end], length=1000), startvalue=sol.t[begin], format=x->@sprintf("%.4f s", x)))

    # return bottom_sg
    tslider = bottom_sg.sliders[2]
    t = Observable(0.0)
    connect!(t, tslider.value)

    # hacky way to add another slider
    tinterval_slider = bottom_sg.layout[3, 2] = IntervalSlider(fig; range=tslider.range, tellheight=true)
    bottom_sg.layout[3, 1] = Label(fig, @lift(@sprintf("%.4f s", $(tinterval_slider.interval)[1])); tellheight=true, halign=:right)
    bottom_sg.layout[3, 3] = Label(fig, @lift(@sprintf("%.4f s", $(tinterval_slider.interval)[2])); tellheight=true, halign=:left)

    # on(tinterval_slider.interval) do interval
    #     if t[] < interval[1]
    #         set_close_to!(tslider, interval[1])
    #     elseif t[] > interval[2]
    #         set_close_to!(tslider, interval[2])
    #     end
    # end

    # t = tslider.value
    ncolorscale=@lift 10^$(bottom_sg.sliders[1].value)

    ##### Color range stuff
    maxrange = lift(nstatesym, rel_to_u0; ignore_equal_values=true) do nstatesym, rel
        extremas = Float64[]
        for i in 1:nv(network)
            try
                x = timeseries(sol, i, nstatesym).x
                if rel
                    x .= x .- x[begin]
                end
                append!(extremas, Base.extrema(x))
            catch
                continue
            end
        end
        (min, max) = isempty(extremas) ? (0,1) : Base.extrema(extremas)
        if min > 0 && max >0
            min = 0
        elseif min<0 && max>0
            m = Base.max(-min, max)
            min = -m
            max =  m
        end
        (min, max)
    end

    ncolorscheme = @lift if ($maxrange)[1] < 0
        # ColorScheme([colorant"blue", colorant"gray50", colorant"red"])
        ColorSchemes.coolwarm
    else
        ColorSchemes.thermal
    end

    ncolorrange = lift(ncolorscale, maxrange; ignore_equal_values=true) do ncolorscale, maxrange
        ncolorscale .* maxrange
    end

    ## Graphplot
    gpax = Axis(gpgrid[1,1])
    args = gparguments(sol, network; t, ncolorscheme, nstatesym, ncolorrange, sel_nodes, rel_to_u0)
    graphplot!(gpax,network; args...)
    hidespines!(gpax)
    hidedecorations!(gpax)

    Colorbar(gpgrid[2,1]; colormap=ncolorscheme, colorrange=ncolorrange, vertical=false, label=@lift("node colors: "*string($nstatesym)), flipaxis=false)

    HOVER_DEFAULT = "Hover node/edge to see info!"
    hover_text = Observable{String}(HOVER_DEFAULT)
    gpgrid[1, 1] = Label(fig, hover_text, tellwidth=false, tellheight=false, justification=:left, halign=:left, valign=:top)

    # delete other interactions on scene
    deregister_interaction!(gpax, :rectanglezoom)

    edgeclick = EdgeClickHandler() do idx, event, axis
        sel_edges[] = idx ∈ sel_edges[] ? delete!(sel_edges[], idx) : push!(sel_edges[], idx)
    end
    register_interaction!(gpax, :eclick, edgeclick)

    nodeclick = NodeClickHandler() do idx, event, axis
        sel_nodes[] = idx ∈ sel_nodes[] ? delete!(sel_nodes[], idx) : push!(sel_nodes[], idx)
    end
    register_interaction!(gpax, :nclick, nodeclick)

    nodehover = NodeHoverHandler() do state, idx, event, axis
        if state
            string = """Node $idx """
        else
            string = HOVER_DEFAULT
        end
        hover_text[] = string
    end
    register_interaction!(gpax, :nhover, nodehover)

    edgehover = EdgeHoverHandler() do state, idx, event, axis
        if state
            string = """Edge $idx """
        else
            string = HOVER_DEFAULT
        end
        hover_text[] = string
    end
    register_interaction!(gpax, :ehover, edgehover)

    #####
    ##### Keyboard interaction for time
    #####
    register_keyboard_interaction!(fig, tslider)

    #####
    ##### Open node plots in new windows
    #####
    on(nselectors[3].clicks) do n
        f = nodeplot_window(sol, tslider, sel_nodes; tlims=tinterval_slider.interval)
        sc = display(GLMakie.Screen(), f)
    end

    fig
end

function nodeplot_window(sol, tslider, sel_nodes; tlims=Observable((sol.t[begin], sol.t[end])))
    fig = Figure(resolution=(1000, 800))

    symgrid = fig[1,1] = GridLayout(tellwidth=false, tellheight=true)
    buttons = symgrid[1,1:5] = [Button(fig, label="_ω"),
                                Button(fig, label="_u_arg"),
                                Button(fig, label="_u_mag"),
                                Button(fig, label="_P"),
                                Button(fig, label="_rocof")]
    symbox = symgrid[1,6] = Textbox(fig, width=150)

    for i in 1:5
        on(buttons[i].clicks) do n
            lab = buttons[i].label[]
            symbox.displayed_string[] = lab
            symbox.stored_string[] = lab
        end
    end

    sym = Observable(:_ω)
    on(symbox.stored_string) do s
        sym[] = Symbol(s)
    end

    ## add menus
    menugrid = fig[2,1] = GridLayout(tellwidth=false, tellheight=true)
    menus = menugrid[1,1:6] = states_dropdown(fig, sol, sel_nodes)
    for menu in filter(m -> m isa Menu, menus)
        on(menu.selection) do sel
            if sel isa String
                symbox.displayed_string[] = sel
                symbox.stored_string[] = sel
            end
        end
    end

    ax = Axis(fig[3, 1])

    plots = Dict{Int, Lines}()

    on(sym; update=true) do sym
        empty!(ax)
        empty!(plots)
        Makie.vlines!(ax, tslider.value; color=:black)
    end

    on(tlims, update=true) do lims
        xlims!(ax, lims)
    end

    legend = nothing
    onany(sel_nodes, sym) do selected, sym
        added   = setdiff(selected, keys(plots))
        removed = setdiff(keys(plots), selected)
        for i in added
            try
                ts = timeseries(sol, i, sym)
                p = lines!(ax, ts; label=string(i), linewidth=3)
                plots[i] = p
            catch e
                # variable not found
            end
        end
        for i in removed
            p = plots[i]
            delete!(plots, i)
            delete!(ax, p)
        end
        if !(isempty(added) && isempty(removed))
            if !isnothing(legend)
                # legend.width[] = 0
                # legend.height[] = 0
                # legend.tellwidth[]=false
                # delete!(legend)
            end
            # isnothing(legend) || delete!(legend)
            if !isempty(plots)
                # legend = Legend(fig[2,2], ax, ":"*string(sym)*" at nodes:"; bgcolor=:white)
            end
            # display(GLMakie.Screen(fig.scene), fig)
        end
    end

    symbox.displayed_string[] = string(sym[])
    symbox.stored_string[] = string(sym[])

    # arrows to move time
    register_keyboard_interaction!(fig, tslider)

    # click to set time
    set_time_interaction = (event::MouseEvent, axis) -> begin
        if event.type === MouseEventTypes.leftclick
            pos = mouseposition(axis.scene)[1]
            set_close_to!(tslider, pos)
            return true
        end
        return false
    end
    register_interaction!(set_time_interaction, ax, :set_time)

    fig
end

function states_dropdown(fig, sol, sel_nodes)
    DEFAULT = "(no option)"
    states = Observable(String[DEFAULT])
    rem_states = Observable(String[DEFAULT])
    meas_states = string.(blockstates(sol, 1; print=false).meas_states)
    on(sel_nodes; update=true) do idxs
        if !isempty(idxs)
            nt = _common_states(sol, idxs)
            states[] = isempty(nt.states) ? [DEFAULT] : string.(nt.states)
            rem_states[] = isempty(nt.rem_states) ? [DEFAULT] : string.(nt.rem_states)
        else
            states[] = [DEFAULT]
            rem_states[] = [DEFAULT]
        end
    end
    l1 = Label(fig, "Common:"; tellwidth=true)
    m1 = Menu(fig; options=meas_states)
    l2 = Label(fig, "States:"; tellwidth=true)
    m2 = Menu(fig; options=states)
    l3 = Label(fig, "Removed:"; tellwidth=true)
    m3 = Menu(fig; options=rem_states)
    return [l1, m1, l2, m2, l3, m3]
end

function register_keyboard_interaction!(fig, tslider)
    steps = 1
    on(fig.scene.events.keyboardbutton) do e
        if e.key == Keyboard.left_shift || e.key == Keyboard.right_shift
            if e.action == Keyboard.press
                steps = 10
            elseif e.action == Keyboard.release
                steps = 1
            end
        elseif e.action == Keyboard.press || e.action == Keyboard.repeat
            if e.key == Keyboard.left
                mv_slider(tslider, -steps)
                return true
            elseif e.key == Keyboard.right
                mv_slider(tslider, steps)
                return true
            end
        end
        return false
    end
end

function mv_slider(slider, steps=1)
    idx = slider.selected_index[]
    set_close_to!(slider, slider.range[][clamp(idx+steps, firstindex(slider.range[]), lastindex(slider.range[]))])
end

function gparguments(sol::ODESolution,
                     network::MetaGraph;
                     t::Observable,
                     nstatesym,
                     ncolorrange,
                     ncolorscheme,
                     rel_to_u0,
                     sel_nodes = Observable(Set{Int}()))

    NV = nv(network)

    markdict = Dict(:load => :rect, :gen => :circle, :syncon => :circle);
    node_marker = [markdict[k] for k in get_prop(network, 1:NV, :type)];

    u0statevec = Observable(Vector{Float64}(undef, NV))
    onany(nstatesym, rel_to_u0) do nstatesym, rel
        if rel
            for i in 1:NV
                u0statevec[][i] = getstate(sol, sol.t[begin], sol.prob.p, i, nstatesym; err=false)
            end
        end
        notify(u0statevec)
    end
    notify(nstatesym) # trigger call to fill u0statevec

    statevec = Observable(Vector{Float64}(undef, NV))
    onany(t, nstatesym, u0statevec) do t, nstatesym, u0statevec
        for i in 1:NV
            statevec[][i] = getstate(sol, t, sol.prob.p, i, nstatesym; err=false)
            if rel_to_u0[]
                statevec[][i] -= u0statevec[i]
            end
        end
        notify(statevec)
    end

    node_color = Observable(Vector{RGB{Float64}}(undef, NV))
    onany(statevec, ncolorrange, ncolorscheme) do statevec, range, scheme
        for i in 1:NV
            node_color[][i] = isnan(statevec[i]) ? RGB(0,0,0) : get(scheme, statevec[i], range)
        end
        notify(node_color)
    end
    notify(t)

    SMALL = 30
    BIG = 50
    node_size = Observable(fill(SMALL, NV))
    onany(sel_nodes) do selected
        fill!(node_size[], SMALL)
        for sel in selected
            node_size[][sel] = BIG
        end
        notify(node_size)
    end
    notify(sel_nodes)

    return (;layout=read_pos_or_spring,
            node_marker,
            node_color,
            node_size,
            node_attr=(;colorrange=ncolorrange, colormap=ncolorscheme));
end

function read_pos_or_spring(g::MetaGraph)
    pos = get_prop(g, 1:nv(g), :pos)
    if any(ismissing.(pos))
        return spring(g)
    else
        return pos
    end
end
