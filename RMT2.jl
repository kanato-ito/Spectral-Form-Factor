#ランダム行列を用いたSFFの計算#
using LinearAlgebra
#using Random
using Plots
using Statistics
using Distributions
using SparseArrays
using ProgressMeter
#using ProgressBars
#using BenchmarkTools
using FLoops

#複素エルミートランダム行列の生成
function rand_M(L)
    M=zeros(ComplexF64,L,L)
    for i in 1:L
        for j in 1:L
            M[i,j]=Distributions.rand(Normal(0,1/sqrt(L)))+im*Distributions.rand(Normal(0,1/sqrt(L)))
        end
    end
    return (M+M')/2
end

#|Z(β+it)|^2
function SFF_u(e,β,t)
    Z_t=sum(exp.(-(β+im*t)*e))
    return Z_t*Z_t'
end
#Z(β)
function get_Z(e,β)
    Z=sum(exp.(-β*e))
    return Z
end

############# SET UP ###############
L=4096;         #行列のサイズ
samples=1200;   #サンプリング数
β=5;            #逆温度
###################################

ts = exp10.(range(-1,stop=5,length=10000));   #10sample, 1000lengthで1分くらい
g_u=zeros(Float64,(samples,length(ts)));
g_d=zeros(Float64,samples);
println("L=$L, ","samples=$samples")

function calc_g()
    prog=Progress(samples,1);
    @floop for  i=1:samples
        M=rand_M(L);
        e=eigvals!(M);
        g_d[i]=get_Z(e,β);
        @inbounds for (j,t) in enumerate(ts)
            g_u[i,j]=SFF_u(e,β,t);
        end
        #print("$i,")
        next!(prog)
    end
end
calc_g()

g_u_mean=mean(g_u,dims=1);#<Z(β+it)*Z(β-it)>
g_d_mean=(mean(g_d))^2;#<Z(β)>^2
gs_mean=g_u_mean/g_d_mean;#g(t;β) = <|Z(β,t)|^2> / <Z(β)>^2

#プロット
plot(ts,gs_mean[:],xaxis=:log,yaxis=:log,xlims=(0.1,10^5),ylims=(3*10^-6,1),xticks=[0.1,1,10,10^2,10^3,10^4,10^5],yticks=[10^-5,10^-4,10^-3,10^-2,10^-1,1],xlabel="t",ylabel="g(t)",title="GUE, L=1024, 1200 samples, β=5, g(t)",minorticks=true,xgrid=false,ygrid=false)
#plot!(ts,gs_mean_2[:],xaxis=:log,yaxis=:log,xlims=(0.1,10^5),ylims=(3*10^-6,1),xticks=[0.1,1,10,10^2,10^3,10^4,10^5],yticks=[10^-5,10^-4,10^-3,10^-2,10^-1,1],xlabel="t",ylabel="g(t)",title="GUE, L=1024, 1200 samples, β=5, g(t)",minorticks=true,xgrid=false,ygrid=false)


################## (ts, gs_mean[:]) プロット用データの保存 ##################
#txtファイルにデータを保存
out=open("data_L=4096_1200samples.txt","w") # ~.txtを作成して開く
println(out,"ts, gs_mean")
Base.print_array(out,hcat(ts,gs_mean[:])) #(ts, gs_mean)を2列にして保存
close(out) #ファイルを閉じる

gs=zeros(samples,length(ts));
for  i in 1:samples
    for j in 1:length(ts)
        gs[i,j]=g_u[i,j]/g_d[i]^2
    end
end

gs_mean_2=mean(gs,dims=1)

out2=open("data_L=4096_1200samples_another.txt","w") # ~.txtを作成して開く
println(out2,"ts, gs_mean")
Base.print_array(out2,hcat(ts,gs_mean_2[:])) #(ts, gs_mean)を2列にして保存
close(out2) #ファイルを閉じる