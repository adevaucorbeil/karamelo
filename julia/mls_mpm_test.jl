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
        $kernel($(args...); threads=$threads, blocks=$blocks)
    end

    return esc(code)
end

quality     = 1  # Use a larger value for higher-res simulations
n_particles = 9000 * quality^2
n_grid      = 128 * quality
dx          = 1 / n_grid
inv_dx      = float(n_grid)
dt          = 1e-4 / quality
p_vol       = (dx * 0.5)^2
p_rho       = 1
p_mass      = p_vol * p_rho
E           = 5e3
nu          = 0.2  # Young's modulus and Poisson's ratio
mu_0        = E / (2 * (1 + nu))
lambda_0    = E * nu / ((1 + nu) * (1 - 2 * nu))  # Lame parameters

x = CuArray{SVector{2, Float64}}(undef, n_particles) # position
v = CuArray{SVector{2, Float64}}(undef, n_particles)  # velocity
F = CuArray{SMatrix{2,2, Float64, 4}}(undef, n_particles)
C = CuArray{SMatrix{2,2, Float64, 4}}(undef, n_particles) # affine velocity field
Jp = CUDA.zeros(Float64, n_particles)  # plastic deformation
material = CUDA.ones(Int64, n_particles)  # material

grid_v = CuArray{SVector{2, Float64}}(undef, n_grid, n_grid) # position
grid_m = CUDA.zeros(Float64, n_grid, n_grid)  # grid node mass


attractor_strength = 0.0
attractor_pos = CUDA.zeros(Float64, 2)


#I = CUDA.ones(Float64, 2, 2)


t = 2*3.141592653589793*CUDA.rand(n_particles)
u = CUDA.rand(n_particles) + CUDA.rand(n_particles)


function reset(x, v, F, Jp, C, t, u, n_particles)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    
    if (i <= n_particles)
        v0 = 2
        group_size = n_particles ÷ 2

        solid = i ÷ group_size
        r = 0.0
        if u[i] > 1
            r = 2-u[i]
        else
            r = u[i]
        end
        r = r * 0.2
        x[i] = SVector{2}(r*CUDA.cos(t[i]) + 0.25 + 0.5 * solid, r*CUDA.sin(t[i]) + 0.25 + 0.5 * solid)
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

function reset_grid(grid_v, grid_m, n_grid)
    ii = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    i = (ii - 1) % n_grid + 1
    j = (ii - 1) ÷ n_grid + 1

    if (i <= n_grid && j <= n_grid)
        grid_v[i,j] = SVector{2}(0, 0)
        grid_m[i,j] = 0
    end
    return nothing
end


function P2G(x, v, F, Jp, C, material, n_particles, grid_v, grid_m, n_grid, dt, inv_dx, mu_0, lambda_0)
    p = (blockIdx().x - 1) * blockDim().x + threadIdx().x

    if (p <= n_particles)
        base = SVector{2}((x[p][1] * inv_dx - 0.5)÷1, (x[p][2] * inv_dx - 0.5)÷1)
        fx = x[p] * inv_dx - base
        ## Quadratic kernels  [http://mpm.graphics   Eqn. 123, with x=fx, fx-1,fx-2]
        wx = SVector{3}(0.5 * (1.5 - fx[1])^2, 0.75 - (fx[1] - 1)^2, 0.5 * (fx[1] - 0.5)^2)
        wy = SVector{3}(0.5 * (1.5 - fx[2])^2, 0.75 - (fx[2] - 1)^2, 0.5 * (fx[2] - 0.5)^2)

        F[p] = (SMatrix{2,2}(1.0, 0.0, 0.0, 1.0) + dt * C[p]) * F[p]  # deformation gradient update
        h = CUDA.max(0.1, CUDA.min(5, CUDA.exp(10 * (1.0 - Jp[p]))))  # Hardening coefficient: snow gets harder when compressed
        if material[p] == 1  # jelly, make it softer
            h = 0.3
        end
        mu = mu_0 * h
        la = lambda_0 * h
        if material[p] == 0  # liquid
            mu = 0.0
        end
        U, sig, V = svd(F[p])
        J = 1.0
        #for d in ti.static(range(2))
        #    new_sig = sig[d, d]
        #    if material[p] == 2  # Snow
        #        new_sig = min(max(sig[d, d], 1 - 2.5e-2),
        #                      1 + 4.5e-3)  # Plasticity
        #    end
        #    Jp[p] *= sig[d, d] / new_sig
        #    sig[d, d] = new_sig
        #    J *= new_sig
        #end
        #if material[p] == 0  # Reset deformation gradient to avoid numerical instability
        #    F[p] = ti.Matrix.identity(float, 2) * ti.sqrt(J)
        #elseif material[p] == 2
        #    F[p] = U @ sig @ V.transpose()  # Reconstruct elastic deformation gradient after plasticity
        #end
        stress = 2 * mu * (F[p] - transpose(inv(F[p]))) + SMatrix{2,2}(1.0, 0.0, 0.0, 1.0) * la * J * (J - 1)
        #stress = 2 * mu * (F[p] - U @ V.transpose()) @ F[p].transpose() + ti.Matrix.identity(float, 2) * la * J * (J - 1)
        #stress = (-dt * p_vol * 4 * inv_dx * inv_dx) * stress
        #affine = stress + p_mass * C[p]
        #for i, j in ti.static(ti.ndrange(3, 3))  # Loop over 3x3 grid node neighborhood
        #    offset = ti.Vector([i, j])
        #    dpos = (offset.cast(float) - fx) * dx
        #    weight = wx[i] * wy[j]
        #    grid_v[base + offset] += weight * (p_mass * v[p] + affine @ dpos)
        #    grid_m[base + offset] += weight * p_mass
        #end
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

@krun n_particles reset(x, v, F, Jp, C, t, u, n_particles)

function substep()
    @krun n_grid^2 reset_grid(grid_v, grid_m, n_grid)
    
    @krun n_particles P2G(x, v, F, Jp, C, material, n_particles, grid_v, grid_m, n_grid, dt, inv_dx, mu_0, lambda_0)

    #@krun n_grid^2 grid_update()

    #@krun n_particles P2G()
end

gravity = CUDA.zeros(Float64, 2)

#for frame in 1:20000
#    substep()
#end
substep()


#FHost = Array(F)
#println(FHost)
println("Done")
