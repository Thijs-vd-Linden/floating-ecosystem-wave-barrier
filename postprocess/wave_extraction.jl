# wave_extraction.jl
#
# Numerical helpers for extracting wave amplitudes from time-domain
# probe signals: fitting a complex amplitude to a sampled time series,
# and decomposing a free-surface signal into incident and reflected
# wave components from three probe positions.

"Least-squares fit of a complex wave amplitude A from a real time series
`signal` sampled at times `t_ss`, assuming signal(t) ≈ Re(A * exp(-iωt))."
function fit_amplitude(signal, t_ss, ω)
    A = [cos.(ω .* t_ss) sin.(ω .* t_ss)]
    c = A \ signal
    return complex(c[1], -c[2])
end

"Decompose a complex free-surface amplitude measured at three probe
positions `x_probes` into incident (A) and reflected (B) wave components,
given wavenumber k. Solves the 3-probe system for A, B in
κ(x) = A·exp(ikx) + B·exp(-ikx)."
function decompose_wave(x_probes, κ_probes, k)
    M = [exp(im*k*x_probes[1])  exp(-im*k*x_probes[1]);
         exp(im*k*x_probes[2])  exp(-im*k*x_probes[2]);
         exp(im*k*x_probes[3])  exp(-im*k*x_probes[3])]
    AB = M \ κ_probes
    return AB[1], AB[2]
end