using LinearAlgebra
using Plots
using CSV
using DataFrames

x_L = 0.0  # left boundary x-coordinate
x_R = 10.0  # right boundary x-coordinate
N = 50  # # of segments composing system
cps = 20  # cells per segment
lam_order = 1  # Polynomial order of lambda
# 0 = constant
# 1 = linear
# 2 = quadratic

# Material 0 Data
lam00 = 101.0 / 20.0  # lambda_0 (x=0)
lam0T = 99.0 / 10.0  # lambda_0 (x=T)
lam0_MaxMult = 2  # Only needed/used if lam_order==2
lam0m = lam0_MaxMult * max(lam00, lam0T)
c0 = 0.0

XS00 = 10.0 / 99.0
XS0T = XS00

XSs00 = c0 * XS00
XSs0T = XSs00

# Material 1 Data
lam10 = 101.0 / 20.0
lam1T = 11.0 / 10.0
lam1_MaxMult = 0.5
lam1m = lam1_MaxMult * min(lam10, lam1T)
c1 = 0.0

XS10 = 100.0 / 11.0
XS1T = XS10

XSs10 = c1 * XS10
XSs1T = XSs10

# Formation of x vectors
Mn = N+1  # # of segments bounds points
mT = cps * N
cT = mT + 1
xn = range(x_L, stop=x_R, length=N+1)  # vectof of segment bound values
x_temp = zeros(cps+1, N)
x_temp[:, 1] = collect(range(x_L, stop=xn[2], length=cps+1))
for i = 2:N
    x_temp[:, i] = collect(range(xn[i], stop=xn[i+1], length=cps+1))
end
x = zeros(mT+1, 1)
x[1:cps+1] = x_temp[:,1]
for i = 2:N
    x[((i-1)*cps+2):(i*cps+1)] = x_temp[2:(cps+1),i]
end

# Data vectors
XS0temp = collect(range(XS00, stop=XS0T, length=N+1))
XS1temp = collect(range(XS10, stop=XS1T, length=N+1))
XSs0temp = collect(range(XSs00, stop=XSs0T, length=N+1))
XSs1temp = collect(range(XSs10, stop=XSs1T, length=N+1))

XS0 = zeros(N, 1)
XS1 = zeros(N, 1)
lam0n = zeros(N, 1)
lam1n = zeros(N, 1)
lamcn = zeros(N, 1)
XSs0 = zeros(N, 1)
XSs1 = zeros(N, 1)
XSR0 = zeros(N, 1)
XSR1 = zeros(N, 1)

if (lam_order == 0)
    for n = 1:N
        XS0[n] = (XS00 + XS0T) / 2.0
        XS1[n] = (XS10 + XS1T) / 2.0
        XSs0[n] = (XSs00 + XSs0T) / 2.0
        XSs1[n] = (XSs10 + XSs1T) / 2.0
        lam0n[n] = (lam00 + lam0T) / 2.0
        lam1n[n] = (lam10 + lam1T) / 2.0
        XSR0[n] = XS0[n] - XSs0[n] / 2.0 + 1.0 / lam0n[n]
        XSR1[n] = XS1[n] - XSs1[n] / 2.0 + 1.0 / lam1n[n]
        lamcn[n] = 1.0 ./ (1.0./lam0n[n] + 1.0./lam1n[n])
    end
elseif (lam_order == 1)
    lam0ntemp = collect(range(lam00, stop=lam0T, length=N+1))
    lam1ntemp = collect(range(lam10, stop=lam1T, length=N+1))
    lamcntemp = 1.0 ./ (1.0 ./ lam0ntemp + 1.0 ./ lam1ntemp)

    for n = 1:N
        XS0[n] = (XS0temp[n+1] + XS0temp[n]) / 2.0
        XS1[n] = (XS1temp[n+1] + XS1temp[n]) / 2.0
        XSs0[n] = (XSs0temp[n+1] + XSs0temp[n]) / 2.0
        XSs1[n] = (XSs1temp[n+1] + XSs1temp[n]) / 2.0
        lam0n[n] = (lam0ntemp[n+1] + lam0ntemp[n]) / 2.0
        lam1n[n] = (lam1ntemp[n+1] + lam1ntemp[n]) / 2.0
        XSR0[n] = XS0[n] - XSs0[n] / 2.0 + 1.0 / lam0n[n]
        XSR1[n] = XS1[n] - XSs1[n] / 2.0 + 1.0 / lam1n[n]
        lamcn[n] = (lamcntemp[n+1] + lamcntemp[n]) / 2.0
    end
elseif (lam_order == 2)
    lam_ep = [lam00 lam0T;
              lam10 lam1T]
    lam_m = [lam0m;
             lam1m]
    dl_0m = zeros(2)
    dl_T0 = zeros(2)
    a = zeros(2)
    b = zeros(2)
    c = zeros(2)
    lam = zeros(2, mT + 1)
    lamn = zeros(2, N + 1)
    for i = 1:2
        dl_0m[i] = lam_ep[i, 1] - lam_m[i]
        dl_T0[i] = lam_ep[i, 2] - lam_ep[i, 1]

        b[i] = -2.0 * dl_0m[i] / x_R * (sqrt(1.0 + dl_T0[i]/dl_0m[i]) + 1.0)
        a[i] = 1.0 / x_R^2 * (dl_T0[i] - b[i] * x_R)
        c[i] = lam_ep[i, 1]

        lam[i, :] = a[i] .* x.^2 .+ b[i] .* x .+ c[i]
        lamn[i, :] = a[i] .* xn.^2 .+ b[i] .* xn .+ c[i]
    end

    lam0ntemp = lamn[1, :]
    lam1ntemp = lamn[2, :]
    lamcntemp = 1.0 ./ (1.0 ./ lam0ntemp + 1.0 ./ lam1ntemp)

    for n = 1:N
        # This assumes cross-sections are at most linearly dependent
        XS0[n] = (XS0temp[n+1] + XS0temp[n]) / 2.0
        XS1[n] = (XS1temp[n+1] + XS1temp[n]) / 2.0
        XSs0[n] = (XSs0temp[n+1] + XSs0temp[n]) / 2.0
        XSs1[n] = (XSs1temp[n+1] + XSs1temp[n]) / 2.0

        dxn = xn[n+1] - xn[n]
        dxn2 = xn[n+1]^2 - xn[n]^2
        dxn3 = xn[n+1]^3 - xn[n]^3
        lam0n[n] = (a[1] * dxn3/3.0 + b[1] * dxn2/2.0 + c[1] * dxn) / dxn
        lam1n[n] = (a[2] * dxn3/3.0 + b[2] * dxn2/2.0 + c[2] * dxn)/dxn

        XSR0[n] = XS0[n] - XSs0[n] / 2.0 + 1.0 / lam0n[n]
        XSR1[n] = XS1[n] - XSs1[n] / 2.0 + 1.0 / lam1n[n]

        lamcn[n] = lam0n[n] * lam1n[n] / (lam0n[n] + lam1n[n])
    end
else
    println("No lambda order available for lam_order = ", lam_order)
    exit()
end

Theta_n = zeros(4, Mn)
Theta = zeros(4, cT)
p_n = zeros(2, Mn)
p = zeros(2, cT)
psi_n = zeros(4, Mn)
psi = zeros(4, cT)

# Boundary conditions

if (lam_order == 0)
    # Use the calculated constant values
    p_n[1, 1] = lam0n[1] / (lam0n[1] + lam1n[1])
else
    # Use the provided left-bound values
    p_n[1, 1] = lam00 / (lam00 + lam10)
end
p_n[2, 1] = 1 - p_n[1, 1]

# Important: Order of Theta_i^{+/-,n} is Different from Report
# Theta = < Theta_0+, Theta_1+, Theta_0-, Theta_1- >

Theta[1, 1] = 1.0  # Theta_0+
Theta[2, 1] = 1.0  # Theta_1+
Theta[3, cT] = 0.0  # Theta_0-
Theta[4, cT] = 0.0  # Theta_1-

M = zeros(4, 4, N)
k = zeros(4, N)
xi0p = zeros(4, N)  # xi_0+
xi0m = zeros(4, N)  # xi_1+
xi1p = zeros(4, N)  # xi_0-
xi1m = zeros(4, N)  # xi_1-

# Obtain eigenvalues and eigenvectors:
for n = 1:N
    M[:, :, n] = [ -XSR0[n] 1.0/lam1n[n] XSs0[n]/2.0 0.0;  # xi_0+
                   1.0/lam0n[n] -XSR1[n] 0.0 XSs1[n]/2.0;  # xi_1+
                   -XSs0[n]/2.0 0.0 XSR0[n] -1.0/lam1n[n];  # xi_0-
                   0.0 -XSs1[n]/2.0 -1.0/lam0n[n] XSR1[n]]  # xi_1-
    d, xixi = eigen(M[:, :, n])
    xixi = transpose(xixi)
    xi0p[:, n] = xixi[:, 1]
    xi1p[:, n] = xixi[:, 2]
    xi0m[:, n] = xixi[:, 3]
    xi1m[:, n] = xixi[:, 4]
    k[:, n] = d
end

C = zeros(4*N, 4*N)
CC = zeros(4, 8, N-1)
B = zeros(4*N, 1)

# Equations for the left boundary
en = [exp(x_L*k[1,1]) exp(x_L*k[2,1]) exp(x_L*k[3,1]) exp(x_L*k[4,1])]
C[1:2, 1:4] = [xi0p[1,1]*en[1] xi0p[2,1]*en[2] xi0p[3,1]*en[3] xi0p[4,1]*en[4];
               xi1p[1,1]*en[1] xi1p[2,1]*en[2] xi1p[3,1]*en[3] xi1p[4,1]*en[4]]
B[1:2] = [Theta[1,1] Theta[2,1]]
# Equations for the right boundary
eN = [exp(x_R*k[1,N]) exp(x_R*k[2,N]) exp(x_R*k[3,N]) exp(x_R*k[4,N])]
C[4*N-1:4*N,4*N-3:4*N] = [xi0m[1,N]*eN[1] xi0m[2,N]*eN[2] xi0m[3,N]*eN[3] xi0m[4,N]*eN[4]
                          xi1m[1,N]*eN[1] xi1m[2,N]*eN[2] xi1m[3,N]*eN[3] xi1m[4,N]*eN[4]]
# The rest of the internal segment boundaries (4N-2 of them, each w/ 4 equations)
for n = 2:N
    enmo = [exp(xn[n]*k[1,n-1]) exp(xn[n]*k[2,n-1]) exp(xn[n]*k[3,n-1]) exp(xn[n]*k[4,n-1])]
    global en = [exp(xn[n]*k[1,n]) exp(xn[n]*k[2,n]) exp(xn[n]*k[3,n]) exp(xn[n]*k[4,n])]

    CC[1,:,n-1] = [xi0p[1,n-1]*enmo[1] xi0p[2,n-1]*enmo[2] xi0p[3,n-1]*enmo[3] xi0p[4,n-1]*enmo[4] -xi0p[1,n]*en[1] -xi0p[2,n]*en[2] -xi0p[3,n]*en[3] -xi0p[4,n]*en[4]]

    CC[2,:,n-1] = [xi1p[1,n-1]*enmo[1] xi1p[2,n-1]*enmo[2] xi1p[3,n-1]*enmo[3] xi1p[4,n-1]*enmo[4] -xi1p[1,n]*en[1] -xi1p[2,n]*en[2] -xi1p[3,n]*en[3] -xi1p[4,n]*en[4]]

    CC[3,:,n-1] = [xi0m[1,n-1]*enmo[1] xi0m[2,n-1]*enmo[2] xi0m[3,n-1]*enmo[3] xi0m[4,n-1]*enmo[4] -xi0m[1,n]*en[1] -xi0m[2,n]*en[2] -xi0m[3,n]*en[3] -xi0m[4,n]*en[4]]

    CC[4,:,n-1] = [xi1m[1,n-1]*enmo[1] xi1m[2,n-1]*enmo[2] xi1m[3,n-1]*enmo[3] xi1m[4,n-1]*enmo[4] -xi1m[1,n]*en[1] -xi1m[2,n]*en[2] -xi1m[3,n]*en[3] -xi1m[4,n]*en[4]]
end

# Arrange CC's sheets into proper coefficient matrix positions
for n = 1:N-1
    C[4*n-1:4*n-1+3, 4*n-3:4*n-3+7] = CC[:, :, n]
end

# Obtain a_j^n coefficients from boundary condition matrix
aa = C\B

# Rearrange a_j^n coefficients for easier index-referencing later:
a = zeros(4, N)
for n = 1:N
    a[:,n] = aa[4*n-3:4*n-3+3]
end

jj = 0
for i = 1:cT
    # Determine which segment we are in
    for j = 1:N
        if x[i] == xn[j]
            global jj += 1
            break
        end
    end

    axi0p = [a[1,jj]*xi0p[1,jj] a[2,jj]*xi0p[2,jj] a[3,jj]*xi0p[3,jj] a[4,jj]*xi0p[4,jj]]
    axi1p = [a[1,jj]*xi1p[1,jj] a[2,jj]*xi1p[2,jj] a[3,jj]*xi1p[3,jj] a[4,jj]*xi1p[4,jj]]
    axi0m = [a[1,jj]*xi0m[1,jj] a[2,jj]*xi0m[2,jj] a[3,jj]*xi0m[3,jj] a[4,jj]*xi0m[4,jj]]
    axi1m = [a[1,jj]*xi1m[1,jj] a[2,jj]*xi1m[2,jj] a[3,jj]*xi1m[3,jj] a[4,jj]*xi1m[4,jj]]
    global en = [exp(x[i]*k[1,jj]) exp(x[i]*k[2,jj]) exp(x[i]*k[3,jj]) exp(x[i]*k[4,jj])]

    Theta[1,i] = axi0p[1]*en[1] + axi0p[2]*en[2] + axi0p[3]*en[3] + axi0p[4]*en[4]  # Theta_0+
    Theta[2,i] = axi1p[1]*en[1] + axi1p[2]*en[2] + axi1p[3]*en[3] + axi1p[4]*en[4]  # Theta_1+
    Theta[3,i] = axi0m[1]*en[1] + axi0m[2]*en[2] + axi0m[3]*en[3] + axi0m[4]*en[4]  # Theta_0-
    Theta[4,i] = axi1m[1]*en[1] + axi1m[2]*en[2] + axi1m[3]*en[3] + axi1m[4]*en[4]  # Theta_1-
end

# Calculation of p_i^n and p_i (x):
for n = 2:Mn
    dxn = xn[n] - xn[n-1]

    # Calculation of p_i^n, i=0,1
    global a = dxn/lamcn[n-1]
    global b = p_n[2,n-1]*exp(-a)
    global c = lamcn[n-1]/lam0n[n-1] * (1-exp(-a))
    p_n[2,n] = b + c
    p_n[1,n] = 1.0 - p_n[2,n]
end
p[1,1] = p_n[1,1]
p[2,1] = p_n[2,1]
n = 1
for i = 2:cT
    # Calculation of p_i(x), i=0,1
    ae = (x[i] - xn[n])/lamcn[n]
    be = p_n[2,n]*exp(-ae)
    ce = lamcn[n]/lam0n[n]*(1 - exp(-ae))
    p[2,i] = be + ce
    p[1,i] = 1.0 - p[2,i]

    for k = n+1:Mn
        if x[i] == xn[k]
            global n += 1
            break
        end
    end
end
p[1,cT] = p_n[1,Mn]
p[2,cT] = p_n[2,Mn]

for i = 1:cT
    psi[1,i] = Theta[1,i]/p[1,i]  # psi_0+
    psi[2,i] = Theta[2,i]/p[2,i]  # psi_1+
    psi[3,i] = Theta[3,i]/p[1,i]  # psi_0-
    psi[4,i] = Theta[4,i]/p[2,i]  # psi_1-
end

flux = zeros(2, cT)
for i = 1:cT
    flux[1,i] = (psi[1,i] + psi[3,i]) / 2.0  # flux_0
    flux[2,i] = (psi[2,i] + psi[4,i]) / 2.0  # flux_1
end

println("Calculation done, writing output...")

"""
plot(x, transpose(flux), title="Flux", label=["Material 0" "Material 1"], yaxis=:log10)
ylabel!("Flux")
xlabel!("x")
ylims!((1e-3, 1e1))
savefig("output.png")

plot(x, transpose(p), reuse=false, title="Volume Fraction", label=["Material 0" "Material 1"])
ylabel!("Volume Fraction")
xlabel!("x")
ylims!((0.0, 1.0))
savefig("volfrac.png")
"""

tabular = DataFrame(x=vec(x), flux_0=vec(transpose(flux[1,:])), flux_1=vec(transpose(flux[2,:])), vf_0=vec(transpose(p[1,:])), vf_1=vec(transpose(p[2,:])))

CSV.write("LP.csv", tabular)

println("Done")
