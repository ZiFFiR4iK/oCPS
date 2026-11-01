using Plots

gr()

struct ADC
	fs::Float64
	M::Int
	Vmin::Float64
	Vmax::Float64
end

const adc_1 = ADC(1000.0, 16, -5.0, 5.0)
const adc_2 = ADC(10.0, 3, 0.0, 3.0)

# число уровней квантования L = 2^M
levels(adc::ADC) = 2 ^ adc.M

# кодирование (1) с ограничением (1*)
function adc_literal(signal, adc::ADC)
	codes = floor.(Int, signal ./ (adc.Vmax - adc.Vmin) .* levels(adc))
	clamp.(codes, 0, levels(adc) - 1) # ограничение от 0 до L-1
end

# восстановление сигнала по (2)
function reconstruct_literal(codes, adc::ADC)
	codes .* (adc.Vmax - adc.Vmin) ./ levels(adc)
end

# построение графика
function save_adc_plot(name, time, analog, adc::ADC; title_text=name,
		analog_seriestype=:path) # последнее - имя графика
	codes = adc_literal(analog, adc)
	quantized = reconstruct_literal(codes, adc)

	plot(time, analog, label="xA(t)", seriestype=analog_seriestype, linewidth=2,
		xlabel="t, s", ylabel="Амплитуда",
		 title="$title_text, fs=$(Int(adc.fs)) Гц, M=$(adc.M)")
	plot!(time, quantized, label="x[n] по (1),(2)", seriestype=:sticks,
		marker=:circle, linewidth=1)
	savefig(joinpath(@__DIR__, "out1.1", name * ".png"))
	println(name, ": первые коды по (1): ", codes[1:min(10, length(codes))])
end

# задаём один ненулевой цифровой отсчёт
function digital_delta(adc::ADC)
	time = -1.0:1 / adc.fs:1.0
	codes = zeros(Int, length(time))
	codes[findfirst(==(0.0), time)] = levels(adc) - 1
	time, codes, reconstruct_literal(codes, adc)
end

# построение дельта
function save_delta_plot(adc::ADC, name)
	time, codes, signal = digital_delta(adc)
	analog = zeros(length(time))
	plot(time, analog, label="xA(t)", linewidth=2,
		 xlabel="t, s", ylabel="Амплитуда", ylims=(-0.1, 1.0),
		 title="Ручной цифровой дельта-импульс, M=$(adc.M)")
	plot!(time, signal, seriestype=:sticks, label="цифровой δ[n]", marker=:circle,
		 xlabel="t, s", ylabel="Амплитуда",
		 title="Ручной цифровой дельта-импульс, M=$(adc.M)")
	savefig(joinpath(@__DIR__, "out1.1", name * ".png"))
	println(name, ": код дельта-импульса = ", maximum(codes))
end

function run()
	# интервал -1:1 с частотой
	time_1 = -1.0:1 / adc_1.fs:1.0
	time_2 = -1.0:1 / adc_2.fs:1.0

	signals_1 = Dict(
		"constant_1" => Float64.(time_1 .>= 0.0),
		"constant_2" => 2 .* Float64.(time_1 .>= 0.0),
		"exp_decay" => 2 .* exp.(-0.5 .* time_1),
		"exp_growth" => 4 .* exp.(0.5 .* time_1),
		"harmonic_1" => 0.5 .* sin.(2pi .* time_1),
		"harmonic_2" => 5 .* sin.(2pi .* 100 .* time_1),
		"harmonic_3" => sin.(2pi .* 0.1 .* time_1 .+ pi / 4),
		"polyharmonic" => cos.(2pi .* time_1) .+
			0.5 .* cos.(2pi .* 3 .* time_1) .+
			0.2 .* cos.(2pi .* 5 .* time_1),
		"window_rectangular" => Float64.(abs.(time_1) .<= 1.0),
		"window_hann" => ifelse.(abs.(time_1) .<= 1.0,
			0.5 .* (1 .+ cos.(2pi .* time_1 ./ 2)), 0.0),
		"window_hamming" => ifelse.(abs.(time_1) .<= 1.0,
			0.54 .+ 0.46 .* cos.(2pi .* time_1 ./ 2), 0.0)
	)

	for (name, signal) in signals_1
		save_adc_plot("adc1_" * name, time_1, signal, adc_1;
			title_text=name, analog_seriestype=name in ("constant_1", "constant_2") ? :steppost : :path)
	end

	# второй
	signals_2 = Dict(
		"constant_1" => Float64.(time_2 .>= 0.0),
		"constant_2" => 2 .* Float64.(time_2 .>= 0.0),
		"exp_decay" => 2 .* exp.(-0.5 .* time_2),
		"exp_growth" => 4 .* exp.(0.5 .* time_2),
		"harmonic_1" => 0.5 .* sin.(2pi .* time_2),
		"harmonic_2" => 5 .* sinpi.(200 .* time_2),
		"harmonic_3" => sin.(2pi .* 0.1 .* time_2 .+ pi / 4),
		"polyharmonic" => cos.(2pi .* time_2) .+
			0.5 .* cos.(2pi .* 3 .* time_2) .+
			0.2 .* cos.(2pi .* 5 .* time_2),
		"window_rectangular" => Float64.(abs.(time_2) .<= 1.0),
		"window_hann" => ifelse.(abs.(time_2) .<= 1.0,
			0.5 .* (1 .+ cos.(2pi .* time_2 ./ 2)), 0.0),
		"window_hamming" => ifelse.(abs.(time_2) .<= 1.0,
			0.54 .+ 0.46 .* cos.(2pi .* time_2 ./ 2), 0.0)
	)

	for (name, signal) in signals_2
		save_adc_plot("adc2_" * name, time_2, signal, adc_2;
			title_text=name, analog_seriestype=name in ("constant_1", "constant_2") ? :steppost : :path)
	end

	save_delta_plot(adc_1, "adc1_delta_manual")
	save_delta_plot(adc_2, "adc2_delta_manual")
end

run()