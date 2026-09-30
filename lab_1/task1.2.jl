using Plots

gr()

struct ADC
	M::Int
	Vmin::Float64
	Vmax::Float64
end

const adc_1 = ADC(16, -5.0, 5.0)
const adc_2 = ADC(3, 0.0, 3.0)

levels(adc::ADC) = 2 ^ adc.M
step(adc::ADC) = (adc.Vmax - adc.Vmin) / levels(adc)

function quantize(signal, adc::ADC)
	# формула 1*
	x_D = floor.(Int, signal ./ (adc.Vmax - adc.Vmin) .* levels(adc))
	x_D_star = min.(max.(x_D, 0), 2^adc.M - 1)
	x_D_star
end

function restore(code, adc::ADC)
	# формула 2
	code .* (adc.Vmax - adc.Vmin) ./ levels(adc)
end

function make_plots(adc::ADC, name::String)
	t = -1.0:0.01:4.0
	x = collect(t)
	code = quantize(x, adc)
	xq = restore(code, adc)
	error = x .- xq

	level_numbers = 0:levels(adc)-1
	level_step = max(1, ceil(Int, levels(adc) / 200))
	shown_levels = restore(collect(level_numbers[1:level_step:end]), adc)

	characteristic = plot(t, x,
		label="x(t) = t", linewidth=2,
		xlabel="t", ylabel="Амплитуда",
		title="Характеристика АЦП: M=$(adc.M), L=$(levels(adc))")
	plot!(characteristic, t, xq,
		label="xq(t)", linewidth=2, seriestype=:steppost)
	hline!(characteristic, shown_levels,
		label="Уровни квантования", color=:gray, linestyle=:dot, alpha=0.35)
	savefig(characteristic, joinpath(@__DIR__, "out1.2", name * "_characteristic.png"))

	error_plot = plot(t, error,
		label="eq(t) = x(t) - xq(t)", linewidth=2,
		xlabel="t", ylabel="Ошибка квантования",
		title="Шум квантования: M=$(adc.M), ΔU=$(step(adc))")
	hline!(error_plot, [0.0], label="0", color=:black, linestyle=:dash)
	savefig(error_plot, joinpath(@__DIR__, "out1.2", name * "_error.png"))

	println(name)
	println("  L = ", levels(adc))
	println("  ΔU = ", step(adc))
	println("  min(eq) = ", minimum(error), ", max(eq) = ", maximum(error))
	println("  первые коды: ", code[1:min(10, length(code))])
end

mkpath(joinpath(@__DIR__, "out1.2"))
make_plots(adc_1, "adc1")
make_plots(adc_2, "adc2")
