XLSX.openxlsx("./Results/results_continuous.xlsx", mode="w") do xf
    sheet = xf[1]
    XLSX.rename!(sheet, "Costs")
    sheet["A1"] = "Nominal cost"
    sheet["B1"] = round(objective_value.(model))
    sheet["A3"] = "Recourse cost"
    sheet["B3"] = round(objective_value(wc))
    sheet["A5"] = "Computation time"
    sheet["B5"] = round(solve_time(model))

    XLSX.addsheet!(xf, "Electricity")
    sheet = xf[2]
    sheet["A1"] = "IL"
    sheet["A2"] = IL
    sheet["A3"] = "EC_nom"
    sheet["A4"] = EC_nom
    sheet["A6"] = "EC nominal"
    for i in I
        sheet["A$(2i+5)"] = "Compressor $i"
        sheet["A$(2i+6)"] = round.(EC[i], digits = 1)
    end
    sheet["A17"] = "EC recourse"
    for i in I
        sheet["A$(2i+16)"] = "Compressor $i"
        sheet["A$(2i+17)"] = round.(EC[i], digits = 1)
    end

    XLSX.addsheet!(xf, "Inventory")
    sheet4 = xf[3]
    sheet4["A2"] = "IV"
    sheet4["A3"] = "Storage tank 1"
    sheet4["A4"] = round.(IV[:C1], digits = 1)
    sheet4["A5"] = "Storage tank 2"
    sheet4["A6"] = round.(IV[tank2_position], digits = 1)
    sheet4["A7"] = "Output gas"
    sheet4["A8"] = round.(IV[:C6], digits = 1)
    sheet4["A10"] = "IV recourse"
    sheet4["A11"] = "Storage tank 1"
    sheet4["A12"] = round.(IVtilde[:C1], digits = 1)
    sheet4["A13"] = "Storage tank 2"
    sheet4["A14"] = round.(IVtilde[tank2_position], digits = 1)
    sheet4["A15"] = "Output gas"
    sheet4["A16"] = round.(IVtilde[:C6], digits = 1)
end

p1 = bar(IL, label = "", 
    linecolor = :match, 
    # color = :springgreen4,
    color = :green, 
    ylims = (0, 450), 
    ylabel = "Electricity requirement [MWh]", 
    legend = :topleft, 
    left_margin = 2Plots.mm,  
    size = (700,370),
    yguidefont = font(13, "times bold"),
    xguidefont = font(11, "times"),
    # title = "(a) Continuous recourse only",
    # titlefont = font(13, "times bold"),
    framestyle = :box,
    # titlelocation = :left,
    tickfont = font(9, "times"),
    # guidefontweight = :bold,
    # legendfontsize = 9
)
# 
plot!(p1, sum(EC[i][:] for i in I), 
    label = "", 
    color = :blue, 
    linewidth = 3,
    ylims = (0, 450), 
    ylabel = "Electricity requirement [MWh]", 
    legend = :topleft,
    )
axis2 = twinx()
plot!(axis2, alphaEC[1:tfin], 
    ylabel = "Price [\$/MWh]", 
    color = :black,  
    linewidth = 2.5,
    label = "",
    right_margin = 5Plots.mm, 
    legend = :topright, 
    yguidefont = font(13, "times bold"),
    xguidefont = font(11, "times"),
    guidefontcolor = :black,
    tickfont = font(9, "times"),
    ylims = (0,120),
    framestyle = :box,
)
plot!(axis2, alphaIL[1:tfin], 
    color = :red, 
    line = (3, :dash), 
    label = "",  
    ylabel = "Price [\$/MWh]", 
    right_margin = 5Plots.mm,
)
xlims!(1, tfin, widen = true)
xticks!(0:6:tfin)           
xlabel!("Time [h]")

# # p2 = bar(ILm, label = "", 
# #     linecolor = :match, 
# #     # color = :springgreen4,
# #     color = :green, 
# #     ylims = (0, 450), 
# #     ylabel = "Electricity requirement [MWh]", 
# #     yguidefont = font(13, "times bold"),
# #     xguidefont = font(11, "times"),
# #     guidefontcolor = :black,
# #     tickfont = font(9, "times"),
# #     legend = :topleft, 
# #     # left_margin = 5Plots.mm,
# #     left_margin = 2Plots.mm, 
# #     title = "(b) Mixed-integer recourse",
# #     titlefont = font(13, "times bold"),
# #     framestyle = :box,
# # )
# # # 
# # plot!(p2, sum(ECm[i][:] for i in I), 
# #     label = "", 
# #     color = :blue, 
# #     linewidth = 3,
# #     ylims = (0, 450), 
# #     ylabel = "Electricity requirement [MWh]", 
# #     legend = :topleft, 
# #     # framestyle = :box,
# #     )
# # axis2 = twinx()
# # plot!(axis2, α_EC[1:tfin], 
# #     ylabel = "Price [\$/MWh]", 
# #     color = :black,  
# #     linewidth = 2.5,
# #     yguidefont = font(13, "times bold"),
# #     xguidefont = font(11, "times"),
# #     framestyle = :box,
# #     # guidefontcolor = :black,
# #     tickfont = font(9, "times"),
# #     label = "",
# #     right_margin = 5Plots.mm, 
# #     legend = :topright, 
# #     ylims = (0,140),
# # )
# # plot!(axis2, α_IL[1:tfin], 
# #     color = :red, 
# #     line = (3, :dash), 
# #     label = "",  
# #     ylabel = "Price [\$/MWh]", 
# #     right_margin = 5Plots.mm,
# #     # framestyle = :box,
# # )
# # xlims!(1, tfin, widen = true)
# # xticks!(0:6:tfin)           
# # xlabel!("Time [h]")

p3 = plot(1,
    color = :black,
    linewidth = 2.5,
    label = "Electricity price",
    legendfontsize = 11,
    framestyle = :none,
    foreground_color_legend = :white,
    legend = :topright
)
plot!(p3, 1,
    color = :red,
    line = (3, :dash),
    label = "Interruptible load price",
    bottom_margin = 2Plots.mm,
    legendfont = font(11,"times"),
)
axis2 = twinx()
bar!(axis2, 1,
    color = :green,
    linecolor = :match,
    label = "Interruptible load provided",
    foreground_color_legend = nothing,
    framestyle = :none,
)
plot!(axis2, 1,
    color = :blue,
    linewidth = 3,
    label = "Nominal electricity consumption",
    legend = :topleft,
    legendfont = font(11,"times"),
)
# plot!(axis2, zeros(tfin),
#     color = :white,
#     linewidth = 6,
#     label = "",
# )

# plt = plot(p1, p2, p3, layout = grid(3,1, heights = (3.95/8,3.95/8, 0.1/8)))
# plt = plot(p1, p2, p3, layout = grid(3,1, heights = (3.47/7,3.47/7, 0.06/7)))
plt = plot(p1, p3, layout = grid(2,1, heights = (0.94,0.06)))