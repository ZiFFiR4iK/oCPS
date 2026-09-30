using Plots
using LinearAlgebra

gr()

const output_dir = joinpath(@__DIR__, "out2.2")
const b7 = [-0.026, -0.6, -0.017, -0.01, -0.022, -0.015, -0.03, -0.023,
	-0.041, -0.036, -0.059, -0.059, -0.094, -0.117, -0.218, -0.629,
	0.629, 0.218, 0.117, 0.094, 0.059, 0.059, 0.036, 0.041,
	0.023, 0.03, 0.015, 0.022, 0.01, 0.017, 0.6, 0.026]
const b8 = [0.0646, -0.1235, 0.1949, -0.1949, 0.1235, -0.0646]
const a8_complex = ComplexF64[1.0000 + 0.0000im, 1.4576 - 0.1880im,
	2.3598 - 0.2122im, 1.8195 - 0.1694im, 1.0772 - 0.0556im,
	0.3218 + 0.0000im]
const a9 = [1.0000, 1.4277, 2.1448, 1.5880, 0.9051, 0.2682]

function fir_system(x, coefficients)
	y = zeros(promote_type(eltype(x), eltype(coefficients)), length(x))
	for n in eachindex(x)
		for i in eachindex(coefficients)
			input_index = n - i + 1
			if input_index >= firstindex(x)
				y[n] += coefficients[i] * x[input_index]
			end
		end
	end
	y
end

function f1(x)
	x
end

function f2(x)
	fir_system(x, fill(1 / 5, 6))
end

function f3(x)
	-1 / 4 .* x
end

function f4(x)
	y = zeros(eltype(x), length(x))
	for n in eachindex(x)
		input_index = n - 4
		if input_index >= firstindex(x)
			y[n] = x[input_index]
		end
	end
	y
end

function f5(x)
	y = zeros(eltype(x), length(x))
	for n in eachindex(x)
		y[n] = x[n]
		if n > firstindex(x)
			y[n] += 2 * y[n - 1]
		end
	end
	y
end

function f6(x)
	0.5 .* abs.(x)
end

function recursive_system(x, b, a)
	y = zeros(promote_type(eltype(x), eltype(b), eltype(a)), length(x))
	for n in eachindex(x)
		for i in eachindex(b)
			input_index = n - i + 1
			if input_index >= firstindex(x)
				y[n] += b[i] * x[input_index]
			end
		end
		for i in 2:length(a)
			output_index = n - i + 1
			if output_index >= firstindex(y)
				y[n] -= a[i] * y[output_index]
			end
		end
	end
	y
end

f7(x) = fir_system(x, b7)
f8(x) = recursive_system(x, b8, a8_complex)
f9(x) = recursive_system(x, b8, a9)

function save_impulse_plot(name, n, y; title_text=name)
	if eltype(y) <: Complex
		plot(n, real.(y), seriestype=:sticks, marker=:circle, label="real",
			xlabel="n", ylabel="Amplitude", title=title_text)
		plot!(n, imag.(y), seriestype=:sticks, marker=:circle, label="imaginary")
	else
		plot(n, y, seriestype=:sticks, marker=:circle, label="h[n]",
			xlabel="n", ylabel="Amplitude", title=title_text)
	end
	savefig(joinpath(output_dir, name * "_impulse.png"))
end

function save_system_plot(name, n, input, output; title_text=name)
	plot(n, input, label="input", xlabel="n", ylabel="Amplitude", title=title_text)
	plot!(n, real.(output), label="output")
	savefig(joinpath(output_dir, name * ".png"))
end

function rms(signal)
	sqrt(sum(abs2, signal) / length(signal))
end

function pole_magnitudes(a)
	order = length(a) - 1
	companion = zeros(ComplexF64, order, order)
	companion[1, :] .= -ComplexF64.(a[2:end])
	for row in 2:order
		companion[row, row - 1] = 1
	end
	abs.(eigvals(companion))
end

mkpath(output_dir)

n_short = 0:15
n_long = 0:40
impulse_short = [Float64.(n .== 0) for n in n_short]
impulse_long = [Float64.(n .== 0) for n in n_long]

systems = [
	("F1", f1, n_short, impulse_short, true, true, "stable"),
	("F2", f2, n_short, impulse_short, true, true, "stable"),
	("F3", f3, n_short, impulse_short, true, true, "stable"),
	("F4", f4, n_short, impulse_short, true, true, "stable"),
	("F5", f5, n_long, impulse_long, true, true, "unstable"),
	("F6", f6, n_short, impulse_short, false, true, "stable"),
	("F7", f7, n_long, impulse_long, true, true, "stable"),
	("F8", f8, n_long, impulse_long, true, true, "unstable (pole outside unit circle)"),
	("F9", f9, n_long, impulse_long, true, true, "stable")
]

println("Свойства систем (нулевые начальные условия):")
for (name, system, n, input, is_linear, is_time_invariant, stability) in systems
	output = system(input)
	save_impulse_plot(name, n, output; title_text="$name: impulse response")
	println("$name: linear=$(is_linear), time-invariant=$(is_time_invariant), stability=$stability")
end

println("F8 pole magnitudes: ", pole_magnitudes(a8_complex))
println("F9 pole magnitudes: ", pole_magnitudes(a9))
println("Импульсная характеристика имеет практический смысл для ЛТИ-систем F1-F5 и F7-F9.")
println("Для F6 она не описывает систему полностью, так как F6 нелинейна.")

input_1 = sin.(2pi * 0.01 .* n_long)
input_2 = sin.(2pi * 0.2 .* n_long)
save_system_plot("F2_low_frequency", n_long, input_1, f2(input_1);
	title_text="F2: input x1[n]")
save_system_plot("F2_high_frequency", n_long, input_2, f2(input_2);
	title_text="F2: input x2[n]")
save_system_plot("F7_low_frequency", n_long, input_1, f7(input_1);
	title_text="F7: input x1[n]")
save_system_plot("F7_high_frequency", n_long, input_2, f7(input_2);
	title_text="F7: input x2[n]")

println("Изменение характеристик входных сигналов:")
for (name, system) in [("F2", f2), ("F7", f7)]
	for (input_name, input) in [("x1", input_1), ("x2", input_2)]
		output = system(input)
		println("$name/$input_name: RMS input=$(round(rms(input), digits=4)), " *
			"RMS output=$(round(rms(output), digits=4)), " *
			"max |output|=$(round(maximum(abs.(output)), digits=4))")
	end
end
println("Частота входных синусоид сохраняется, а амплитуда и фаза изменяются фильтрами F2 и F7.")