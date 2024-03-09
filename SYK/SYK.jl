using LinearAlgebra
using Random
using Plots
using Statistics
using Distributions
using SparseArrays
using ProgressMeter
#using ProgressBars
#using BenchmarkTools
using FLoops


cr=[0 1
    0 0];
an=[0 0
    1 0];
id=Matrix(I,2,2);
id2=[-1 0
    0 1];
    function Table(M,n)
      A=[]
      for i=1:n
          push!(A,M)
      end
      return A
  end
  function K1(n)
      A=push!(Table(id,n-1),cr)
      for i=1:Nc-n
          A=push!(A,Table(id2,Nc-n)[i])
      end
      return A
  end    
  function K2(n)
      A=push!(Table(id,n-1),an)
      for i=1:Nc-n
          A=push!(A,Table(id2,Nc-n)[i])
      end
      return A
  end    
  function KroneckerProduct1(n)
      A=K1(n)[1]
      for i=1:Nc-1
          A=kron(A,K1(n)[i+1])
      end
      return A
  end        
  function KroneckerProduct2(n)
      A=K2(n)[1]
      for i=1:Nc-1
          A=kron(A,K2(n)[i+1])
      end
      return A
  end        
  function c(n)
      c=sparse(KroneckerProduct1(n))
      return c
  end
  function c_d(n)
      cd=sparse(KroneckerProduct2(n))
      return cd
  end
  function make_ψ()
      ψ=[]
      for n in 1:Nc
          push!(ψ,1/√2*(c(n)+c_d(n)))
      end
      for n in 1:Nc
          push!(ψ,1/√2*(-im*c(n)+im*c_d(n)))
      end
      return ψ
  end

#H_sq, Eig_sq, sum_eq = [], [], [] ベンチマークのプロット用
#for N in 4:2:24                   ベンチマークのプロット用

#フェルミオンの数#
N=24;
Nc=div(N,2);
  ψ=make_ψ();
  function make_ψ_2()
     ψ_2=[ψ[i]*ψ[j] for i in 1:N, j in 1:N]
      return ψ_2
  end
  ψ_2=make_ψ_2();
q=4;
J=1;
β=5;
function Hamiltonian()
    js=zeros(Float64,(N,N,N,N));
    js=rand!(Normal(0,sqrt(J^2*factorial((q-1))/(N^(q-1)))),js);
    H=Array(sum(sum(ψ_2[i,j]*sum(sum(js[i,j,k,l]*ψ_2[k,l] for l in k+1:N) for k in j+1:N-1) for j in i+1:N-2) for i in 1:N-3));
    return H
end
#SFF#
function SFF_u(e,β,t)
    #Z=sum(exp.(-β*e))
    Z_t=sum(exp.(-(β+im*t)*e))
    return Z_t*Z_t'
end
function get_Z(e,β)
    Z=sum(exp.(-β*e))
    return Z
end


samples=1200;
ts = exp10.(range(-1,stop=5,length=10000));   #10sample, 1000lengthで1分くらい#10000で十分？（mathematicaの結果と同じくらいに見える）
g_u=zeros(Float64,(samples,length(ts)));
g_d=zeros(Float64,samples);
gs=zeros(samples,length(ts));

#=function test()
    for  i=1:samples
        A=zeros(ComplexF32,(2^Nc,2^Nc))#目盛消費を抑えるため精度を落とす-> N=28で500s->370sになった
   A.=Hamiltonian()
   e=eigvals(A)
    Z=sum(exp.(-β*e))
   for (j,t) in enumerate(ts)
       g=SFF(e,β,t,Z)
       gs[i,j]=g
   end
end
end
=#
function test2()
    prog=Progress(samples,1);
   @floop for  i=1:samples
   A=Hamiltonian()
   e=eigvals!(A);
        g_d[i]=get_Z(e,β);
        @inbounds for (j,t) in enumerate(ts)
            g_u[i,j]=SFF_u(e,β,t);
        end
next!(prog)
end
end
println("N=$N, ","samples=$samples")
#=@time test()
gs_mean=mean(gs,dims=1)
plot(ts,gs_mean[:],xaxis=:log,yaxis=:log,xlabel="t",ylabel="g(t)",title="SFF in the SYK")
=#
@time test2()
g_u_mean=mean(g_u,dims=1);#<Z(β+it)*Z(β-it)>
g_d_mean=(mean(g_d));#<Z(β)>^2
gs_mean=g_u_mean/g_d_mean;#g(t;β) = <|Z(β,t)|^2> / <Z(β)>^2
plot(ts,gs_mean[:],xaxis=:log,yaxis=:log,xlabel="t",ylabel="g(t)",title="SFF in the SYK(N=$N, $samples samples)")

out=open("data_N=$N"*"_$samples"*"samples.txt","w") # ~.txtを作成して開く
println(out,"ts, gs_mean")
Base.print_array(out,hcat(ts,gs_mean[:])) #(ts, gs_mean)を2列にして保存
close(out) #ファイルを閉じる
#=

gs=zeros(samples,length(ts));
for  i in 1:samples
    for j in 1:length(ts)
        gs[i,j]=g_u[i,j]/g_d[i]^2
    end
end

gs_mean_2=mean(gs,dims=1)
print("N=$N","_pp")

out2=open("data_N=$N _1200samples_another.txt","w") # ~.txtを作成して開く
println(out2,"ts, gs_mean")
Base.print_array(out2,hcat(ts,gs_mean_2[:])) #(ts, gs_mean)を2列にして保存
close(out2) #ファイルを閉じる
=#

#サンプル保存#

g_u_mean=g_u_mean[:]
pushfirst!(g_u_mean,g_d_mean)
out3=open("data_N=$N"*"_$samples"*"samples_gs_part3.txt","w") # ~.txtを作成して開く
println(out3,"1行目はZ(β), 2行目からはZ(β+it)*Z(β-it)")
Base.print_array(out3,g_u_mean) #(ts, gs_mean)を2列にして保存
close(out3)

################## (ts, gs_mean[:]) プロット用データの保存 ##################
#txtファイルにデータを保存
#=out=open("data_N=16.txt","w") # ~.txtを作成して開く
println(out,"ts, gs_mean")
Base.print_array(out,hcat(ts,gs_mean[:])) #(ts, gs_mean)を2列にして保存
close(out) #ファイルを閉じる
=#

########### txtファイルからプロット用データを読み込む ###########
#=dataorg = DelimitedFiles.readdlm("data_N=16.txt",skipstart=1) # datorgにデータを読み込み(*1行目はラベルなのでskipstart=1で飛ばす)
ts = dataorg[:,1]
gs_mean = dataorg[:,2]
plot(ts,gs_mean[:],xaxis=:log,yaxis=:log,xlabel="t",ylabel="g(t)",title="SFF in the SYK(N=$N, $samples samples)")
=#

################## ベンチマークのプロット ##################
#push!(H_sq, a.time)
#push!(Eig_sq, b.time)
#push!(sum_eq,a.time+b.time)
#end

#plot(4:2:24,H_sq,label="create H")
#plot!(4:2:24,Eig_sq,label="diagonalization")
#p=plot!(4:2:24,sum_eq,xlabel="N",ylabel="t[s]",title="Benchmarks of Diagonalization",label="sum")
#pにベンチマークのプロットを保存
#savefig(p,".pdf")
