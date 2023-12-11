using GLMakie; GLMakie.activate!() #framerate=1.0, scalefactor=1.0)
using SorptionModels
using MembraneBase


START_SLICE = 1
RECTANGLE_MARGIN = 0.01


COLORMAP_CONTOUR = :turbo
COLORMAP_HEATMAP = :rainbow1


function transform_range(original_index, original_range, new_range)
    percent_complete = (original_index - original_range[1]) / (original_range[end] - original_range[1])
    data_coordinate = percent_complete * (new_range[end] - new_range[1]) + new_range[1]
    return data_coordinate
end

"""
Find min and max values in an Array, 3x faster than extrema
"""
function minmax(x)
    a = b = x[1]
    for i in 2:length(x)
        if x[i] > b 
            b = x[i]
        elseif x[i] < a 
            a = x[i]
        end
    end
    a, b
end

function format_axis_spine!(axis::Axis, spinecolor, spinethickness)
    axis.bottomspinecolor = spinecolor
    axis.topspinecolor = spinecolor
    axis.leftspinecolor = spinecolor
    axis.rightspinecolor = spinecolor
    axis.spinewidth = spinethickness
end

function get_dualmode_sinfs(isotherms)
    dual_mode_models = [fit_model(DualMode(), isotherm) for isotherm in isotherms]
    dual_mode_infinite_dilution_solubilities = infinite_dilution_solubility.(dual_mode_models)
    return dual_mode_infinite_dilution_solubilities
end

function get_nelf_sinfs(nelf_params, bulk_phase_characteristic_params, isotherms)
    pol_densities = polymer_density.(isotherms)

    # generate the equations of state necessary to predict properties
    # bulk_phase_eos = [SL(params...) for params in bulk_phase_characteristic_params]
    # polymer_phase_eos = [SL([[nelf_param, bulk_param] for (nelf_param, bulk_param) in zip(nelf_params, bulk_params)]...) for bulk_params in bulk_phase_characteristic_params]
    
    polymer_phase_eos = [SorptionModels.construct_binary_sl_eosmodel(nelf_params, bulk_params, 0) for bulk_params in bulk_phase_characteristic_params]

    # generate the nelf models for each isotherm
    nelf_models = [NELFModel(polymer_phase, pol_dens) for (polymer_phase, pol_dens) in zip(polymer_phase_eos, pol_densities)]
    
    # predict infinite dilution solubilities and parity them (dual mode is assumped to basically be an empirical representation of the isotherm)
    nelf_infinite_dilution_solubilities = [infinite_dilution_solubility(nelf_model, t)[1] for (nelf_model, t) in zip(nelf_models, temperature.(isotherms))]
    return nelf_infinite_dilution_solubilities
end

function sanchez_lacome_parameter_fit_error_explorer(isotherms, isotherm_names, bulk_phase_params, pstars, tstars, rhostars; errgrid = missing, fit_kij=false)
    if ismissing(errgrid)
        errgrid, _, _ = SorptionModels.nelf_characteristic_parameter_error_map(isotherms, bulk_phase_params, rhostars=rhostars, pstars=pstars, tstars=tstars; fit_kij) 
    end
    data = log.(errgrid)

    ln_dm_sinfs = log.(get_dualmode_sinfs(isotherms))

    nx, ny, nz = size(data)
    minerr, maxerr = minmax(data)
    gminx, gminy, gminz = Tuple(argmin(data))
    
    if !ismissing(pstars)
        p_range = pstars  
        t_range = tstars
        ρ_range = rhostars
    else
        p_range = range(-0.5, nx - 0.5, length = nx)  
        t_range = range(-0.5, ny - 0.5, length = ny)
        ρ_range = range(-0.5, nz - 0.5, length = nz)
    end

    constant_ρ_local_minima_indices = [Tuple(argmin(@view data[:, :, idx])) for idx in 1:nz]  # (x, y) coords
    constant_t_local_minima_indices = [Tuple(argmin(@view data[:, idx, :])) for idx in 1:ny]  # (x, z) coords
    constant_ρ_local_minima = [(transform_range(idx, 1:nx, p_range), transform_range(idy, 1:ny, t_range)) for (idx, idy) in constant_ρ_local_minima_indices]
    constant_t_local_minima = [(transform_range(idx, 1:nx, p_range), transform_range(idz, 1:nz, ρ_range)) for (idx, idz) in constant_t_local_minima_indices]


    function z_rectangle_indicator_positions(height::Int64)
        cartesian_height = transform_range(height, 1:nz, ρ_range)
        return [
            (p_range[1]-RECTANGLE_MARGIN, t_range[1]-RECTANGLE_MARGIN, cartesian_height), (p_range[end]+RECTANGLE_MARGIN, t_range[1]-RECTANGLE_MARGIN, cartesian_height),
            (p_range[end]+RECTANGLE_MARGIN, t_range[1]-RECTANGLE_MARGIN, cartesian_height), (p_range[end]+RECTANGLE_MARGIN, t_range[end]+RECTANGLE_MARGIN, cartesian_height), 
            (p_range[end]+RECTANGLE_MARGIN, t_range[end]+RECTANGLE_MARGIN, cartesian_height), (p_range[1]-RECTANGLE_MARGIN, t_range[end]+RECTANGLE_MARGIN, cartesian_height), 
            (p_range[1]-RECTANGLE_MARGIN, t_range[end]+RECTANGLE_MARGIN, cartesian_height), (p_range[1]-RECTANGLE_MARGIN, t_range[1]-RECTANGLE_MARGIN, cartesian_height)
        ]
    end

    function y_rectangle_indicator_positions(width::Int64)
        cartesian_width = transform_range(width, 1:ny, t_range)
        return [
            (p_range[1]-RECTANGLE_MARGIN, cartesian_width, ρ_range[1]-RECTANGLE_MARGIN), (p_range[end]+RECTANGLE_MARGIN, cartesian_width, ρ_range[1]-RECTANGLE_MARGIN),
            (p_range[end]+RECTANGLE_MARGIN, cartesian_width, ρ_range[1]-RECTANGLE_MARGIN), (p_range[end]+RECTANGLE_MARGIN, cartesian_width, ρ_range[end]+RECTANGLE_MARGIN), 
            (p_range[end]+RECTANGLE_MARGIN, cartesian_width, ρ_range[end]+RECTANGLE_MARGIN), (p_range[1]-RECTANGLE_MARGIN, cartesian_width, ρ_range[end]+RECTANGLE_MARGIN), 
            (p_range[1]-RECTANGLE_MARGIN, cartesian_width, ρ_range[end]+RECTANGLE_MARGIN), (p_range[1]-RECTANGLE_MARGIN, cartesian_width, ρ_range[1]-RECTANGLE_MARGIN)
        ]
    end

    function yz_line_indicator_positions(height::Int64, width::Int64)
        cartesian_width = transform_range(width, 1:ny, t_range)
        cartesian_height = transform_range(height, 1:nz, ρ_range)
        return [
            (p_range[1]-RECTANGLE_MARGIN * 2, cartesian_width, cartesian_height), (p_range[end]+RECTANGLE_MARGIN * 2, cartesian_width, cartesian_height),
        ]
    end

    # define the layout
    fig = Figure( backgroundcolor=:lightsteelblue1)
    fig.layout[1, 1] = main_view_and_control_layout = GridLayout()

    main_view_and_control_layout[1, 1] = main_view_axis = Axis3(main_view_and_control_layout[1, 1], viewmode = :fitzoom, azimuth=-0.25π, aspect = (1, 1, 1),
        xlabel = "P*", ylabel = "T*", zlabel="ρ*") # is fig the real parent here?
    fig[1, 2] = sl_ρ = Slider(main_view_and_control_layout[1, 2], range = 1:nz, startvalue = START_SLICE, horizontal=false)
    fig[2, 1] = sl_t = Slider(main_view_and_control_layout[2, 1], range = 1:ny, startvalue = START_SLICE, horizontal=true)
    fig[2, 2] = reset_button = Button(main_view_and_control_layout[2, 2], label = "Global Minima")
    
    main_view_and_control_layout[1, 2] = constant_ρstar_axis = Axis(fig, aspect=nx/ny, xlabel = "P*", ylabel = "T*")
    main_view_and_control_layout[2, 1] = constant_tstar_axis = Axis(fig, aspect=nx/nz, xlabel = "P*", ylabel = "ρ*")
    main_view_and_control_layout[2, 2] = lines_layout = GridLayout()
    
    lines_layout[1, 1] = constant_ρt_axis = Axis(fig, xlabel = "P*", ylabel = "log(RSS)",
        backgroundcolor = :gray16,)
    lines_layout[1, 2] = sinf_axis = Axis(fig, xlabel = "DualMode Sinf", ylabel = "NELF Sinf", aspect=1)
    limits!(sinf_axis, (0, maximum(ln_dm_sinfs)+2), (0, maximum(ln_dm_sinfs)+2))

    
    # format the axes
    ρt_axis_margins_constant_y = (maxerr - minerr) * 0.05
    ρt_axis_margins_constant_x = (p_range[end] - p_range[1]) * 0.1
    GLMakie.ylims!(constant_ρt_axis, (minerr - ρt_axis_margins_constant_y, maxerr + ρt_axis_margins_constant_y))
    GLMakie.xlims!(constant_ρt_axis, (p_range[1] - ρt_axis_margins_constant_x, p_range[end] + ρt_axis_margins_constant_x))
    format_axis_spine!(constant_ρstar_axis, :red, 4)
    format_axis_spine!(constant_tstar_axis, :blue, 4)
    format_axis_spine!(constant_ρt_axis, :black, 4)
    # backgroundcolor = :gray90

    # fill data with interaction-able stuff
    GLMakie.contour!(main_view_axis, p_range, t_range, ρ_range, data, colormap=COLORMAP_CONTOUR, nan_color = :black, levels=10)
   
    constant_ρ_slice(val::Integer) = @view data[:, :, val]
    constant_t_slice(val::Integer) = @view data[:, val, :]
    constant_ρt_slice(ρ_idx::Integer, t_idx::Integer) = @view data[:, t_idx, ρ_idx]

    get_constant_ρ_minima(idx) = constant_ρ_local_minima[idx]
    get_constant_t_minima(idx) = constant_t_local_minima[idx]
    function get_constant_ρt_minima(ρ_idx, t_idx) 
        minimum_index = argmin(constant_ρt_slice(ρ_idx, t_idx))
        minimum_p = transform_range(minimum_index, 1:nx, p_range)
        return (minimum_p, data[minimum_index, t_idx, ρ_idx])
    end
    
    chosen_position = lift(sl_ρ.value, sl_t.value) do ρ_idx, t_idx
        p_value, _ = get_constant_ρt_minima(ρ_idx, t_idx)
        t_value = transform_range(t_idx, 1:ny, t_range)
        ρ_value = transform_range(ρ_idx, 1:nz, ρ_range)
        return p_value, t_value, ρ_value
    end
    GLMakie.scatter!(main_view_axis, chosen_position, markersize=20, color = :magenta,  overdraw=true)

    chosen_position_text = lift(chosen_position) do pos
        p_value = round(pos[1], digits = 0)
        t_value = round(pos[2], digits = 0)
        ρ_value = round(pos[3], digits = 2)
        return "P⋆ = $p_value\nT⋆ = $t_value\nρ⋆ = $ρ_value"
    end

    ln_nelf_sinf_values = lift(chosen_position) do pos
        nelf_params = [pos..., 100000]
        return log.(get_nelf_sinfs(nelf_params, bulk_phase_params, isotherms))
    end

    constant_ρ_data = lift(constant_ρ_slice, sl_ρ.value)
    constant_t_data = lift(constant_t_slice, sl_t.value)
    constant_ρt_data = lift(constant_ρt_slice, sl_ρ.value, sl_t.value)

    constant_ρ_minima = lift(get_constant_ρ_minima, sl_ρ.value)
    constant_t_minima = lift(get_constant_t_minima, sl_t.value)
    constant_ρt_minima = lift(get_constant_ρt_minima, sl_ρ.value, sl_t.value)



    GLMakie.heatmap!(constant_ρstar_axis, p_range, t_range, constant_ρ_data, colormap=COLORMAP_HEATMAP)
    GLMakie.heatmap!(constant_tstar_axis, p_range, ρ_range, constant_t_data, colormap=COLORMAP_HEATMAP)
    GLMakie.contour!(constant_ρstar_axis, p_range, t_range, constant_ρ_data, color=:white)
    GLMakie.contour!(constant_tstar_axis, p_range, ρ_range, constant_t_data, color=:white)
    GLMakie.lines!(constant_ρt_axis, p_range, constant_ρt_data, color=:purple, linewidth=8) # todo why is there not an axis here?

    GLMakie.scatter!(constant_ρstar_axis, constant_ρ_local_minima, color=:gray16)
    GLMakie.scatter!(constant_ρstar_axis, constant_ρ_minima, color=:white)
    GLMakie.scatter!(constant_tstar_axis, constant_t_local_minima, color=:gray16)
    GLMakie.scatter!(constant_tstar_axis, constant_t_minima, color=:white)
    
    GLMakie.scatter!(constant_ρt_axis, constant_ρt_minima, markersize = 20, color = :magenta)
    text!(constant_ρt_axis, 
        chosen_position_text, position = constant_ρt_minima, 
        align = (:center, :center), justification = :center, 
        offset = (0, 50), color = :white
    )

    GLMakie.scatter!(sinf_axis, ln_dm_sinfs, ln_nelf_sinf_values)
    sinf_textpositions = lift(ln_nelf_sinf_values) do ln_nelf_sinf_vals
        collect(zip(ln_dm_sinfs, ln_nelf_sinf_vals))
    end
    lines!(sinf_axis, [0, maximum(ln_dm_sinfs)], [0, maximum(ln_dm_sinfs)])
    text!(sinf_axis, 
        isotherm_names, position = sinf_textpositions, 
        align = (:left, :center), justification = :center, 
        offset = (15, 0), color = :black,  fontsize = 12
    )
    
    # create the slider visualizers 
    z_rect_positions = lift(z_rectangle_indicator_positions, sl_ρ.value)
    y_rect_positions = lift(y_rectangle_indicator_positions, sl_t.value)
    x_line_positions = lift(yz_line_indicator_positions, sl_ρ.value, sl_t.value)

    linesegments!(main_view_axis, z_rect_positions, linewidth=5, color=:red)
    linesegments!(main_view_axis, y_rect_positions, linewidth=5, color=:blue)
    linesegments!(main_view_axis, x_line_positions, linewidth=5, color=:magenta)

    # add the local minimum button behavior
    function reset_sliders_to_global_minima(_=missing)
        set_close_to!(sl_ρ, gminz)
        set_close_to!(sl_t, gminy)
    end
    on(reset_sliders_to_global_minima, reset_button.clicks)


    return fig
end


# view parameter space navigator
function parameter_space_navigator(isotherms, bulk_phase_characteristic_params, iso_names; fit_kij=false)
    rhostars=1.6:0.03:1.85
    pstars=300:15:1200
    tstars=300:15:900

    errgrid = SorptionModels.nelf_characteristic_parameter_error_map(
        isotherms, bulk_phase_characteristic_params, rhostars=rhostars, pstars=pstars, tstars=tstars; 
        nan_on_failure=true, verbose=true
    )

    for item in errgrid
        if isnan(item)
            @error "NaNs found in grid!"
        end
    end
    fig = sanchez_lacome_parameter_fit_error_explorer(isotherms, iso_names, bulk_phase_characteristic_params, pstars, tstars, rhostars, errgrid = @view errgrid[:, :, :]; fit_kij)
    return fig
end


tpbo_ch4_5c = IsothermData(; 
    partial_pressures_mpa = [0, 0.051022404, 0.097858902, 0.172285293, 0.272711361, 0.371012386, 0.475068139], 
    concentrations_cc = [0, 7.37368523, 12.0433614, 17.76552202, 23.28373709, 27.50367509, 31.07457011],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_ch4_20c = IsothermData(; 
    partial_pressures_mpa = [0, 0.054409635, 0.101117776, 0.17570818, 0.275518841, 0.373215869], 
    concentrations_cc = [0, 5.084334416, 8.57780684, 12.92845101, 17.55678063, 21.17207162],
    temperature_k = 293.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_ch4_35c = IsothermData(; 
    partial_pressures_mpa = [0, 0.05676287, 0.103720596, 0.177868877, 0.268361442, 0.371119351, 0.478013248], 
    concentrations_cc = [0, 3.741224553, 6.311976164, 9.748565324, 13.23714075, 16.47955269, 19.49981169],
    temperature_k = 308.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_5c = IsothermData(; 
    partial_pressures_mpa = [0, 0.012474405, 0.043391439, 0.11294484, 0.260641173, 0.563818055, 1.023901167, 1.633212939, 2.105369806, 2.615898101], 
    concentrations_cc = [0, 15.19436518, 33.5283133, 54.43691918, 76.32730082, 99.6986546, 121.6846495, 142.8977897, 157.6456881, 174.186204],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_20c = IsothermData(; 
    partial_pressures_mpa = [0, 0.018051448, 0.043568055, 0.07237145, 0.134854673, 0.239969739, 0.386544681], 
    concentrations_cc = [0, 11.64079279, 22.18344653, 30.5524273, 43.09494749, 56.42684262, 67.28267947],
    temperature_k = 293.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_35c = IsothermData(; 
    partial_pressures_mpa = [0, 0.031852099, 0.066294896, 0.104825903, 0.142384952, 0.18000105, 0.217614795, 0.29451087, 0.380689171], 
    concentrations_cc = [0, 11.7485186, 19.85668992, 26.48815715, 31.5333868, 36.0129803, 39.61060375, 45.73820821, 51.29191881],
    temperature_k = 308.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_co2_50c = IsothermData(; 
    partial_pressures_mpa = [0, 0.023513454, 0.050773712, 0.080001807, 0.144376557, 0.249710838, 0.396483131], 
    concentrations_cc = [0, 5.93630284, 11.36554572, 15.98552528, 23.62447856, 33.06426088, 42.47173453],
    temperature_k = 323.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_n2_5c = IsothermData(; 
    partial_pressures_mpa = [0, 0.068906986, 0.181377336, 0.424951374, 0.731306858, 1.064696014, 1.413103086], 
    concentrations_cc = [0, 2.252715738, 5.601581157, 11.31054253, 16.87930294, 21.39238669, 25.17075548],
    temperature_k = 278.15, rho_pol_g_cm3 = 1.3937
    )
tpbo_n2_50c = IsothermData(; 
    partial_pressures_mpa = [0, 0.269930265, 0.705742173, 1.060688385, 1.42192415, 1.813024602, 2.228349107], 
    concentrations_cc = [0, 2.435212223, 5.677614879, 8.139676474, 10.6450967, 12.90356804, 14.82380991],
    temperature_k = 323.15, rho_pol_g_cm3 = 1.3937
    )


isotherms = [tpbo_ch4_5c, tpbo_ch4_20c, tpbo_ch4_35c, tpbo_co2_5c, tpbo_co2_20c, tpbo_co2_35c, tpbo_co2_50c, tpbo_n2_5c, tpbo_n2_50c]
    
char_co2 = [630, 300, 1.515, 44]
char_ch4 = [250, 215, 0.500, 16.04]
char_n2 = [160, 145, 0.943, 28.01]
bulk_phase_char_params = [char_ch4, char_ch4, char_ch4, char_co2, char_co2, char_co2, char_co2, char_n2, char_n2]


iso_names = ["ch4_5c", "ch4_20c", "ch4_35c", "co2_5c", "co2_20c", "co2_35c", "co2_50c", "n2_5c", "n2_50c"]
parameter_space_navigator(isotherms, bulk_phase_char_params, iso_names; fit_kij=true)