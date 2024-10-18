using Catalyst
using Random
using BlackBoxOptim
using Flux
using CSV
using DataFrames
using JLD2
using Sobol
include("train_NN.jl")
include("utils_GTM.jl")
include("constants.jl")
result_dir = joinpath(root_path, "result/synthetic")
model_type = "GTM"
d = 5
intensity = 1
logranges = [  1.0 15.0
                0.1 10.0
                1.0 15.0
                0.01 10.0
                0.1  400.0
             ]

function loss_hellinger_map_baye(x::AbstractVector, model,hist_yobs, tt=tt)
    bufs = zeros(Threads.nthreads())
    Threads.@threads for i in 1:length(tt)
        ps = x
        # pred = pred_pdf(model, ps, 0:length(hist_yobs))
        pred = pred_pdf_infe(model, ps, 0:length(hist_yobs))
        bufs[Threads.threadid()] += hellinger2(pred, hist_yobs)
    end
    global AllSolutions
    global AllLoss
    push!(AllSolutions, deepcopy(x))
    loss_result =sum(bufs)
    push!(AllLoss, loss_result)
    return loss_result
    
end

@load joinpath(DATA_DIR, "test_$model_type.jld2") params_arr_test sim_results_pro_test
@load joinpath(MODELWEIGHT_DIR, "model_stats_prob.jld2") model

estimates = []
s = SobolSeq(logranges[:,1] ,logranges[:,2])
ps_train = [ Sobol.next!(s) for i in 1:100 ]
start_point_num = 3
params_num = 10
params_name = ["kon","ron","koff","roff","mu"]

for i in eachindex(sim_results_pro_test[1:params_num])
    estimates = []
    AllSolutions = Vector{Float64}[]
    AllLoss = Vector{Float64}()
    hist_yobs = sim_results_pro_test[i]
    tt = [0,0]
    yobs = []
    p = Progress(1)
    Threads.@threads for j in 1:start_point_num
        opt_result = bboptimize(p -> loss_hellinger_map_baye(p, model,hist_yobs,tt),ps_train[j]; SearchRange = [ tuple(logranges[i,:]...) for i in 1:d ], TraceMode=:silent)
        push!(estimates, best_candidate(opt_result))
        ProgressMeter.next!(p)
    end
    AllLoss=exp.(-AllLoss)
    AllSolutions_matrix = reduce(vcat,transpose.(AllSolutions))
    paramsArr_matrix = reduce(vcat,transpose.(params_arr_test[1:params_num]))

    AllSolution_df = DataFrame(AllSolutions_matrix,:auto)
    AllSolution_df[:, :tauOn] = AllSolution_df[!,1]./AllSolution_df[!,2]
    AllSolution_df[:, :tauOff] = AllSolution_df[!,3]./AllSolution_df[!,4]
    AllSolution_df[:, :tauOnRatio] = AllSolution_df[:, :tauOn] ./(AllSolution_df[:, :tauOn]+AllSolution_df[:, :tauOff])
    AllSolution_df[:, :bf] = 1 ./ (AllSolution_df[:, :tauOn] .+ AllSolution_df[:, :tauOff])
    AllSolution_df[:,:bs] = AllSolution_df[:, :tauOn] .* AllSolution_df[!,5]
    AllSolution_df[:,:loss] = AllLoss
    CSV.write( joinpath(result_dir, "posteriorDist/AllSolution_$i.csv"),AllSolution_df)
end
