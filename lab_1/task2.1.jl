using Plots

gr()

const Fs = 1000
const n = -4000:4000
const t = n ./ Fs
const output_dir = joinpath(@__DIR__, "out2.1")

# формула 6
function convolution(x, h)
	result = zeros(Float64, length(x) + length(h) - 1)
	for output_index in eachindex(result)
		for k in eachindex(x)
			h_index = output_index - k + 1
			if 1 <= h_index <= length(h)
				result[output_index] += x[k] * h[h_index]
			end
		end
	end
	result
end

# формула корреляции
function correlation(x, h)
	lags = -(length(h) - 1):(length(x) - 1)
	result = zeros(Float64, length(lags))
	for (result_index, m) in enumerate(lags)
		for n in eachindex(x)
			h_index = n - m
			if 1 <= h_index <= length(h)
				result[result_index] += x[n] * h[h_index]
			end
		end
	end
	result
end

function save_convolution(name, x, h, title_text)
	y = convolution(x, h)
	time = ((first(n) + first(n)):((last(n) + last(n)))) ./ Fs
	plot(time, y, label="x * h", xlabel="t, s", ylabel="Amplitude",
		title=title_text)
	savefig(joinpath(output_dir, name * ".png"))
	println(name, ": длина результата = ", length(y))
end

function save_correlation(name, x, h, title_text)
	y = correlation(x, h)
	lags = (-(length(h) - 1):(length(x) - 1)) ./ Fs
	plot(lags, y, label="x star h", xlabel="lag, s", ylabel="Correlation",
		title=title_text)
	savefig(joinpath(output_dir, name * ".png"))
	println(name, ": длина результата = ", length(y), ", максимум = ", maximum(y))
end

mkpath(output_dir)

x1 = Float64.((0 .<= t) .& (t .<= 2))
x2 = ifelse.(0 .<= t .< 1, t, ifelse.(1 .<= t .<= 2, 2 .- t, 0.0))
x3 = ifelse.(0 .<= t .<= 2, sin.(2pi .* t), 0.0)

println("Количество отсчётов: ", length(t))
println("Ненулевые отсчёты: ", sum(x1 .!= 0), ", ", sum(x2 .!= 0), ", ", sum(x3 .!= 0))

save_convolution("x1_x2_convolution", x1, x2, "x1 * x2")
save_convolution("x1_x3_convolution", x1, x3, "x1 * x3")
save_convolution("x2_x3_convolution", x2, x3, "x2 * x3")
save_convolution("x1_x1_convolution", x1, x1, "x1 * x1")
save_convolution("x3_x3_convolution", x3, x3, "x3 * x3")

save_correlation("x1_x2_correlation", x1, x2, "x1 star x2")
save_correlation("x3_x3_correlation", x3, x3, "x3 star x3")

function read_signal(filename)
	parse.(Float64, readlines(filename))
end

signal_1 = read_signal(joinpath(@__DIR__, "signal_1.csv"))
signal_2 = read_signal(joinpath(@__DIR__, "signal_2.csv"))

plot(1:length(signal_1), signal_1, label="signal 1", xlabel="Sample number",
	ylabel="Amplitude", title="Signal 1")
savefig(joinpath(output_dir, "signal_1.png"))

plot(1:length(signal_2), signal_2, label="signal 2", xlabel="Sample number",
	ylabel="Amplitude", title="Signal 2")
savefig(joinpath(output_dir, "signal_2.png"))

csv_correlation = correlation(signal_1, signal_2)
csv_lags = -(length(signal_2) - 1):(length(signal_1) - 1)
maximum_index = argmax(csv_correlation)
maximum_lag = csv_lags[maximum_index]

plot(csv_lags, csv_correlation, label="signal 1 star signal 2",
	xlabel="Lag, samples", ylabel="Correlation", title="CSV signal correlation")
scatter!([maximum_lag], [csv_correlation[maximum_index]],
	label="maximum", marker=:star, markersize=8, color=:red)
vline!([maximum_lag], label="maximum lag", linestyle=:dash, color=:red)
savefig(joinpath(output_dir, "signals_correlation.png"))

println("Длина корреляции CSV: ", length(csv_correlation))
println("Максимум корреляции CSV: ", csv_correlation[maximum_index])
println("Положение максимума: ", maximum_lag, " отсчётов")