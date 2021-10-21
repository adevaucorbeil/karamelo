using CUDA, StaticArrays, Printf


macro krun(ex...)
    len = ex[1]
    call = ex[2]

    args = call.args[2:end]

    @gensym kernel config threads blocks
    code = quote
        local $kernel = @cuda launch=false $call
        local $config = launch_configuration($kernel.fun)
        local $threads = min($len, $config.threads)
        local $blocks = cld($len, $threads)
        print("Launch with ")
        print($threads)
        print(" threads and ")
        print($blocks)
        println(" blocks")
        $kernel($(args...); threads=$threads, blocks=$blocks)
    end

    return esc(code)
end

quality     = 1  # Use a larger value for higher-res simulations
n_particles = 9000*quality^2
n_grid      = 128*quality
dx          = 1.0 / n_grid
inv_dx      = float(n_grid)
dt          = 1e-4 / quality
p_vol       = (dx*0.5)^2
p_rho       = 1
p_mass      = p_vol*p_rho
E           = 5e3
nu          = 0.2  # Young's modulus and Poisson's ratio
mu_0        = E / (2*(1 + nu))
lambda_0    = E*nu / ((1 + nu)*(1 - 2*nu))  # Lame parameters

x = CuArray{SVector{2, Float64}}(undef, n_particles) # position
v = CuArray{SVector{2, Float64}}(undef, n_particles)  # velocity
F = CuArray{SMatrix{2,2, Float64, 4}}(undef, n_particles)
C = CuArray{SMatrix{2,2, Float64, 4}}(undef, n_particles) # affine velocity field
Jp = CUDA.zeros(Float64, n_particles)  # plastic deformation
material = CUDA.ones(Int64, n_particles)  # material

grid_vx = CUDA.zeros(Float64, n_grid, n_grid) # position along x
grid_vy = CUDA.zeros(Float64, n_grid, n_grid) # position along y
grid_m = CUDA.zeros(Float64, n_grid, n_grid)  # grid node mass

gravity_x = 0.0
gravity_y = -1.0
attractor_strength = 0.0
attractor_pos_x = 0
attractor_pos_y = 0


#I = CUDA.ones(Float64, 2, 2)


t = 2*3.141592653589793*CUDA.rand(n_particles)
u = CUDA.rand(n_particles) + CUDA.rand(n_particles)


function reset(x, v, F, Jp, C, t, u, n_particles)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    
    if (i <= n_particles)
        v0 = 2
        group_size = n_particles÷2

        solid = (i-1)÷group_size
        r = 0.0
        if u[i] > 1
            r = 2-u[i]
        else
            r = u[i]
        end
        r = r*0.2
        x[i] = SVector{2}(r*CUDA.cos(t[i]) + 0.25 + 0.5*solid, r*CUDA.sin(t[i]) + 0.25 + 0.5*solid)
        ##material[i] = 1#i // group_size  # 0: fluid 1: jelly 2: snow
        if (solid == 0)
            v[i] = SVector{2}(v0, v0)
        else
            v[i] = SVector{2}(-v0, -v0)
        end
        
        F[i] = SMatrix{2,2}(1.0, 0.0, 0.0, 1.0)
        Jp[i] = 1
        C[i] = SMatrix{2,2}(0.0, 0.0, 0.0, 0.0)
    end
    return nothing
end

function reset_grid(grid_vx, grid_vy, grid_m, n_grid)
    ii = (blockIdx().x - 1)*blockDim().x + threadIdx().x

    i = (ii - 1)%n_grid + 1
    j = (ii - 1)÷n_grid + 1

    if (i <= n_grid && j <= n_grid)
        grid_vx[i,j] = 0.0
        grid_vy[i,j] = 0.0
        grid_m[i,j] = 0.0
    end
    return nothing
end


function P2G(x, v, F, Jp, C, material, n_particles, grid_vx, grid_vy, grid_m, n_grid, dt, inv_dx, mu_0, lambda_0, p_mass, p_vol)
    p = (blockIdx().x - 1)*blockDim().x + threadIdx().x

    if (p <= n_particles)
        base = SVector{2, Int}((x[p][1]*inv_dx - 0.5)÷1, (x[p][2]*inv_dx - 0.5)÷1)
        #@cuprintf("base=[%ld, %ld]\n", base[1], base[2])
        fx = x[p]*inv_dx - base
        ## Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        wx = SVector{3}(0.5*(1.5 - fx[1])^2, 0.75 - (fx[1] - 1)^2, 0.5*(fx[1] - 0.5)^2)
        wy = SVector{3}(0.5*(1.5 - fx[2])^2, 0.75 - (fx[2] - 1)^2, 0.5*(fx[2] - 0.5)^2)

        F[p] = (SMatrix{2,2}(1.0, 0.0, 0.0, 1.0) + dt*C[p])*F[p]  # deformation gradient update
        h = CUDA.max(0.1, CUDA.min(5, CUDA.exp(10*(1.0 - Jp[p]))))  # Hardening coefficient: snow gets harder when compressed
        if material[p] == 1  # jelly, make it softer
            h = 0.3
        end
        mu = mu_0*h
        la = lambda_0*h
        if material[p] == 0  # liquid
            mu = 0.0
        end
        U, sig, V = svd(F[p])
        J = 1.0
        #for d in 1:2
        #    new_sig = sig[d, d]
        #    if material[p] == 2  # Snow
        #        new_sig = min(max(sig[d, d], 1 - 2.5e-2), 1 + 4.5e-3)  # Plasticity
        #    end
        #    Jp[p] *= sig[d, d] / new_sig
        #    sig[d, d] = new_sig
        #    J *= new_sig
        #end
        if material[p] == 0  # Reset deformation gradient to avoid numerical instability
            F[p] = SMatrix{2, 2}(1.0, 0.0, 0.0, 1.0)*sqrt(J)
        elseif material[p] == 2
            F[p] = U*sig*transpose(V)  # Reconstruct elastic deformation gradient after plasticity
        end
        # stress = 2*mu*(F[p] - transpose(inv(F[p]))) + SMatrix{2,2}(1.0, 0.0, 0.0, 1.0)*la*J*(J - 1)
        stress = (-dt*p_vol*4*inv_dx*inv_dx)*(2*mu*(F[p] - U*transpose(V))*transpose(F[p]) + SMatrix{2,2}(1.0, 0.0, 0.0, 1.0)*la*J*(J - 1))
        affine = stress + p_mass*C[p]
        for i in 1:3
            for j in 1:3
                index = base - SVector{2, Int}(i - 1, j - 1)
                dpos = SVector{2}(i - 1, j - 1) - x[p]*inv_dx + base
                weight = wx[i]*wy[j]
                dv = weight*(p_mass*v[p] + affine*dpos)
                @inbounds @atomic grid_vx[index[1], index[2]] += dv[1]
                @inbounds @atomic grid_vy[index[1], index[2]] += dv[2]
                @inbounds @atomic grid_m[index[1], index[2]] += weight*p_mass
            end
        end
    end
    return nothing
end

@inline function svd(A)
    x = A[1, 1] + A[2, 2]
    y = A[2, 1] - A[1, 2]
    scale = 1/sqrt(x*x + y*y)
    c = x*scale
    s = y*scale
    R = SMatrix{2,2}(c, -s, s, c)
    S = transpose(R)*A
    s1 = 0.0
    s2 = 0.0
    if (abs(S[1, 2]) < 1e-5)
        c = 1.0;
        s = 0.0;
    else
        tao = 0.5*(S[1, 1] - S[2, 2])
        w = CUDA.sqrt(tao*tao + S[1, 2]*S[1, 2])
        t = 0.0
        if (tao > 0.0)
            t = S[1, 2]/(tao + w);
        else
            t = S[1, 2]/(tao - w);
        end
        c = 1/CUDA.sqrt(t*t + 1);
        s = -t*c;
        s1 = c*c*S[1, 1] - 2*c*s*S[1, 2] + s*s*S[2, 2];
        s2 = s*s*S[1, 1] + 2*c*s*S[1, 2] + c*c*S[2, 2];
    end
       
    if (s1 < s2)
        tmp = s1
        s1 = s2;
        s2 = tmp;
        V = SMatrix{2,2}(-s, c, -c, -s);
    else
        V = SMatrix{2,2}(c, s, -s, c);
    end
    U = R*V;
    sig = SMatrix{2,2}(s1, 0.0, 0.0, s2);
    return U, sig, V
end

function grid_update(grid_vx, grid_vy, grid_m, n_grid, gravity_x, gravity_y, attractor_strength, attractor_pos_x, attractor_pos_y, dt, dx)
    ii = (blockIdx().x - 1)*blockDim().x + threadIdx().x

    i = (ii - 1)%n_grid + 1
    j = (ii - 1)÷n_grid + 1
    if i <= n_grid && j <= n_grid && grid_m[i, j] > 0  # No need for epsilon here
        grid_vx[i, j] = 1/grid_m[i, j]*grid_vx[i, j]  # Momentum to velocity
        grid_vy[i, j] = 1/grid_m[i, j]*grid_vy[i, j]  # Momentum to velocity
        grid_vx[i, j] += dt*gravity_x*30  # gravity
        grid_vy[i, j] += dt*gravity_y*30  # gravity
        dist = SVector{2}(attractor_pos_x, attractor_pos_y) - dx*SVector{2}(i, j)
        grid_vx[i, j] += dist[1]/(0.01 + sqrt(dist[1]*dist[1] + dist[2]*dist[2]))*attractor_strength*dt*100
        grid_vy[i, j] += dist[2]/(0.01 + sqrt(dist[1]*dist[1] + dist[2]*dist[2]))*attractor_strength*dt*100
        if i <= 3 && grid_vx[i, j] < 0
            grid_vx[i, j] = 0.0
        end
        if i >= n_grid - 3 && grid_vx[i, j] > 0
            grid_vx[i, j] = 0.0
        end
        if j <= 3 && grid_vy[i, j] < 0
            grid_vy[i, j] = 0.0
        end
        if j >= n_grid - 3 && grid_vy[i, j] > 0
            grid_vy[i, j] = 0.0
        end
    end
    return nothing
end

function G2P(x, v, C, n_particles, grid_vx, grid_vy, dt, inv_dx)
    p = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    if (p <= n_particles)
        base = SVector{2, Int}((x[p][1]*inv_dx - 0.5)÷1, (x[p][2]*inv_dx - 0.5)÷1)
        fx = x[p]*inv_dx - base
        ### Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        wx = SVector{3}(0.5*(1.5 - fx[1])^2, 0.75 - (fx[1] - 1)^2, 0.5*(fx[1] - 0.5)^2)
        wy = SVector{3}(0.5*(1.5 - fx[2])^2, 0.75 - (fx[2] - 1)^2, 0.5*(fx[2] - 0.5)^2)
        new_v = SVector{2}(0.0, 0.0)
        new_C = SMatrix{2,2}(0.0, 0.0, 0.0, 0.0)
        for i in 1:3
            for j in 1:3
                index = base - SVector{2, Int}(i - 1, j - 1)
        #        #dpos = SVector{2}(i, j) - fx
                dpos = SVector{2}(i - 1, j - 1) - x[p]*inv_dx + base
                #        
                if index[1] < 0 || index[1] > 128
                    @cuprintf("p=%ld, index[1]=%ld, x=[%f, %f]\n", p, index[1], x[p][1], x[p][2])
                end
                if index[2] < 0 || index[2] > 128
                    @cuprintf("p=%ld, index[2]=%ld, x=[%f, %f]\n", p, index[2], x[p][1], x[p][2])
                end
                g_v = SVector{2}(grid_vx[index[1], index[2]],grid_vy[index[1], index[2]])
                weight = wx[i] * wy[j]
                new_v += weight * g_v
                new_C += 4*inv_dx*weight*SMatrix{2,2}(g_v[1]*dpos[1], g_v[1]*dpos[2], g_v[2]*dpos[1], g_v[2]*dpos[2])#1.0, 1.0, 1.0, 1.0)
            end
        end
        ##@cuprintf("C[%ld]=[%f, %f, %f, %f]", p, new_C[1], new_C[2], new_C[3], new_C[4])
        v[p] = new_v
        C[p] = new_C
        x[p] += dt*new_v  # advection
    end
    return nothing
end

@krun n_particles reset(x, v, F, Jp, C, t, u, n_particles)
xHost = Array(x)
pp = 165
@printf("x[%ld]=[%f %f]\n", pp, xHost[pp][1], xHost[pp][2])

function substep()
    @krun n_grid^2 reset_grid(grid_vx, grid_vy, grid_m, n_grid)
    
    @krun n_particles P2G(x, v, F, Jp, C, material, n_particles, grid_vx, grid_vy, grid_m, n_grid, dt, inv_dx, mu_0, lambda_0, p_mass, p_vol)

    @krun n_grid^2 grid_update(grid_vx, grid_vy, grid_m, n_grid, gravity_x, gravity_y, attractor_strength, attractor_pos_x, attractor_pos_y, dt, dx)

    @krun n_particles G2P(x, v, C, n_particles, grid_vx, grid_vy, dt, inv_dx)
    xHost = Array(x)
    @printf("x[%ld]=[%f %f]\n", pp, xHost[pp][1], xHost[pp][2])
end

gravity = CUDA.zeros(Float64, 2)

for frame in 1:20#0000
    println("frame=$frame")
    substep()
end

print("C=")
println(Array(C))
#FHost = Array(F)
#println(FHost)
println("Done")
