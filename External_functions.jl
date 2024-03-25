# Ref.: https://www.juliapackages.com/p/ipopt
function my_callback(
    alg_mod::Cint,
    iter_count::Cint,
    obj_value::Float64,
    inf_pr::Float64,
    inf_du::Float64,
    mu::Float64,
    d_norm::Float64,
    regularization_size::Float64,
    alpha_du::Float64,
    alpha_pr::Float64,
    ls_trials::Cint,
 )
    append!(objHist, obj_value)
    append!(gHist, inf_pr)
    # return `true` to keep going, or `false` to terminate the optimization.
    return true
 end