function numberstotex(d)
    join(map(n -> begin
                 if divrem(n, 1)[2] == 0
                     string(Int(n))
                 else
                     string(n)
                 end
             end, d'[:]), " & ")
end
function datatotextable(raw, rearranged)
    s = "\\begin{tabular}{l|c|cccccc||cccccc}\n"
    s *= " \\multicolumn{2}{c|}{} & \\multicolumn{2}{c}{A} &
    \\multicolumn{2}{c}{B} & \\multicolumn{2}{c||}{C} &
    \\multicolumn{2}{c}{K\$_1\$}  & \\multicolumn{2}{c}{K\$_2\$} &
    \\multicolumn{2}{c}{K\$_3\$} \\\\"
    s *= "bird & protocoll & P & M & P & M & P & M & f\$_1\$ & f\$_2\$  & f\$_1\$  & f\$_2\$  & f\$_1\$  & f\$_2\$
          \\\\\\midrule\n"
    for (k, (c, f, d)) in raw
        s *= "$k & $(join(map(string, c)))-$(join(map(string, f))) & $(numberstotex(d)) & $(numberstotex(rearranged[k].cached)) \\\\\n"
     end
     s *= "\\end{tabular}"
end
function datatotextable2(raw, rearranged)
    s = "\\begin{tabular}{l|c|ccc||ccc}\n"
    s *= "bird & protocoll & A & B & C & K\$_1\$  & K\$_2\$ & K\$_3\$ \\\\\\midrule\n"
    for (k, (c, f, d)) in raw
        s *= "$k & $(join(map(string, c)))-$(join(map(string, f))) & $(numberstotex(d)) & $(numberstotex(rearranged[k].cached)) \\\\\n"
     end
     s *= "\\end{tabular}"
end

function plotindividualdata(data, title;
                            isaxis = true, isticks = true,
                            y_axis_line_style = "{<->}",
                            ylabel = "f\$_2\$\\hspace{1.2cm} f\$_1\$",
                            xticklabels = "")
    if length(data) == 3
        append!(data, zeros(3)) # this is a hack to deal with exp 2
    end
    @pgf p = Plot(Coordinates([1:3; 1:3], [data[1:3]; -data[4:end]]))
    if isaxis
        a1 = @pgf Axis({ybar, bar_width = "10pt", title = title,
                        width = "3.5cm", ymin = - maximum(data[4:end]) - 1,
                        ymax = maximum(data[1:3]) + 1,
                        axis_x_line = "middle", axis_y_line = "left",
                        xmin = .5, xmax = 3.5, xtick = [1, 2, 3],
                        xticklabels = xticklabels,
                        xticklabel_shift = "2mm",
                        y_axis_line_style = y_axis_line_style,
                        x_axis_line_style = "{-}",
                        ytick = "\\empty",
                        ylabel = ylabel,
                        ylabel_style = {yshift = "0pt"},
                        yticklabels = ["f\$_2\$", "f\$_1\$"],
                        yticklabel = "{\\scriptsize
                                   \\ifdim\\tick pt < 0pt
                                   \\pgfmathparse{-1*\\tick}
                                   \\pgfmathprintnumber{\\pgfmathresult}
                                   \\else
                                   \\pgfmathprintnumber{\\tick}
                                   \\fi}" }, p)
        if !isticks
            a2["xtick"] = "\\empty"
        end
        TikzPicture(a1)
    else
        p
    end
end

function plotallindividualdata(data; exp = 1, dir, xlabelat = "",
                               titles = Dict(), kwargs...)
    for (k, (c, f, d)) in data
        td = TikzDocument(plotindividualdata(d, haskey(titles, k) ? titles[k] : string(k); xticklabels =
                                             (string(k) == xlabelat ? ["K\$_1\$", "K\$_2\$", "K\$_3\$"] : ""),
                                             kwargs...))
        PGFPlotsX.savetex(joinpath(dir, "data$(replace(string(k), " " => ""))-$exp.tex"), td,
                          include_preamble = false)
    end
end

function plotcomparison(y, models; title = "", ymax = 1, legend = true,
                            legendentries = ["bird independent", "bird dependent",
                                            "food independent", "food dependent"],
                            legendtoname = "modelcomparisonlegend",
                            altpattern = "crosshatch",
                            colors = ["red!70", "blue!70"],
                            exp2 = false
                           )
    plots = Any[]
    coords = [0, 1, 1.6, 2.2, 2.8]
    idx1 = findall(x -> x.bird_indep && x.food_indep, models)
    push!(plots, @pgf Plot({fill = colors[2], draw = colors[2]},
                           Coordinates(coords[1:length(idx1)], y[idx1])))
    idx2 = findall(x -> !x.bird_indep && x.food_indep, models)
    push!(plots, @pgf Plot({fill = colors[1], draw = colors[1]},
                           Coordinates(coords[1:length(idx2)], y[idx2])))
    if !exp2
    idx3 = findall(x -> x.bird_indep && !x.food_indep, models)
    push!(plots, @pgf Plot({pattern = altpattern, draw = colors[2],
                            pattern_color = colors[2]},
                           Coordinates(coords[1:length(idx3)], y[idx3])))
    idx4 = findall(x -> !x.bird_indep && !x.food_indep, models)
    push!(plots, @pgf Plot({pattern = altpattern, draw = colors[1],
                            pattern_color = colors[1]},
                           Coordinates(coords[1:length(idx3)], y[idx4])))
    end
    a = @pgf Axis({ybar, bar_width = "14pt", xmin = -.5, xmax = 3.2, width = "\\textwidth",
                   height = "7cm", line_width = "0pt", xtick = coords, ymax = ymax,
                   major_tick_length = 0,
                          xticklabels = ["comp. independent", "comp. dependent",
                                         (exp2 ? "" : "\\hspace{10mm}") * " CCH",
                                         (exp2 ? "" : "\\hspace{9mm}") * " FPH 1",
                                         (exp2 ? "" : "\\hspace{9mm}") * " FPH 2"],
                          ymode = "log", log_origin = "infty",
                          title = title}, plots...)
    if legend
        prepend!(a.contents, ["\\addlegendimage{color = $(colors[2])}
                              \\addlegendentry[color = black]{$(legendentries[1])}
                         \\addlegendimage{color = $(colors[1])}
                         \\addlegendentry[color = black]{$(legendentries[2]) }" *
                         (exp2 ? "" : "\\addlegendimage{fill = black, line width = 0pt}
                         \\addlegendentry[color = black]{$(legendentries[3]) },
                         \\addlegendimage{pattern = $altpattern, line width = 0pt}
                         \\addlegendentry[color = black]{$(legendentries[4]) }")])
        a["legend_to_name"] = legendtoname
        a["legend_columns"] = "-1"
        a["legend_style"] = "draw = none, nodes = {inner sep = 8pt}"
    end
    a
end

