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

function reset(x, v, F, Jp, C, n_particles, quality)
    i = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    if (i <= n_particles)
        group_size = n_particles÷3
        i_in_group = (i - 1)%group_size
        dim_size = 4*quality
        row = Float64(i_in_group%dim_size)
        col = Float64(i_in_group/dim_size)
        x[i] = SVector{2}(row/dim_size*0.2 + 0.3 + 0.1*((i - 1)÷group_size),
                          col/dim_size*0.2 + 0.05 + 0.29*((i - 1)÷group_size))
        v[i] = SVector{2}(0.0, 0.0)
        F[i] = SMatrix{2,2}(1.0, 0.0, 0.0, 1.0)
        Jp[i] = 1
        C[i] = SMatrix{2,2}(0.0, 0.0, 0.0, 0.0)
    end
    return nothing
end

function reset_grid(grid_vx, grid_vy, grid_m, n_grid)
    ii = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    if (ii <= n_grid^2)
        i = (ii - 1)%n_grid + 1
        j = (ii - 1)÷n_grid + 1
        grid_vx[i, j] = 0.0
        grid_vy[i, j] = 0.0
        grid_m[i,j] = 0.0
    end
    return nothing
end

function P2G(x, v, F, Jp, C, material, n_particles, grid_vx, grid_vy, grid_m, dt, dx, inv_dx, mu_0, lambda_0, p_vol, p_mass)
    p = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    if (p <= n_particles)
    #for p in 1:n_particles
        base = SVector{2, Int}(Int((x[p][1]*inv_dx - 0.5)÷1), Int((x[p][2]*inv_dx - 0.5)÷1))
        fx = x[p]*inv_dx - base
        wx = SVector{3}(0.5*(1.5 - fx[1])^2, 0.75 - (fx[1] - 1)^2, 0.5*(fx[1] - 0.5)^2)
        wy = SVector{3}(0.5*(1.5 - fx[2])^2, 0.75 - (fx[2] - 1)^2, 0.5*(fx[2] - 0.5)^2)

        F[p] = (SMatrix{2,2}(1.0, 0.0, 0.0, 1.0) + dt*C[p])*F[p]  # deformation gradient update
        h = CUDA.max(0.1, CUDA.min(5, CUDA.exp(10*(1.0 - Jp[p]))))  # Hardening coefficient: snow gets harder when compressed
        if material[p] == 1  # jelly, make it softer
            h = 0.3
        end
        mu = mu_0*h
        la = lambda_0*h
        #if material[p] == 0  # liquid
        #    mu = 0.0
        #end
        U, sig, V = svd(F[p])
        J = 1.0
        for d in 1:2
            new_sig = sig[d, d]
            if material[p] == 2  # Snow
                new_sig = min(max(sig[d, d], 1 - 2.5e-2), 1 + 4.5e-3)  # Plasticity
            end
            Jp[p] *= sig[d, d]/new_sig
            if d == 1
                sig = SMatrix{2, 2}(new_sig, sig[2, 1], sig[1, 2], sig[2, 2])
            else
                sig = SMatrix{2, 2}(sig[1, 1], sig[2, 1], sig[1, 2], new_sig)
            end
            J *= new_sig
        end
        if material[p] == 0  # Reset deformation gradient to avoid numerical instability
            F[p] = SMatrix{2, 2}(1.0, 0.0, 0.0, 1.0)*sqrt(J)
        elseif material[p] == 2
            F[p] = U*sig*transpose(V)  # Reconstruct elastic deformation gradient after plasticity
        end
        stress = 2*mu*(F[p] - U*transpose(V))*transpose(F[p]) + SMatrix{2,2}(1.0, 0.0, 0.0, 1.0)*la*J*(J - 1)
        stress *= -dt*p_vol*4*inv_dx*inv_dx
        affine = stress + p_mass*C[p]
        for i in 0:2
            for j in 0:2
                index = base + SVector{2, Int}(i + 1, j + 1)
                dpos = (SVector{2}(i, j) - fx)*dx
                weight = wx[i + 1]*wy[j + 1]
                dv = weight*(p_mass*v[p] + affine*dpos)
                @inbounds @atomic grid_vx[index[1], index[2]] += dv[1]
                @inbounds @atomic grid_vy[index[1], index[2]] += dv[2]
                @inbounds @atomic grid_m[index[1], index[2]] += weight*p_mass
            end
        end
    end
    return nothing
end

function grid_update(grid_vx, grid_vy, grid_m, n_grid, gravity_x, gravity_y, dt)
    ii = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    if (ii <= n_grid^2)
        i = (ii - 1)%n_grid + 1
        j = (ii - 1)÷n_grid + 1
        if grid_m[i, j] > 0  # No need for epsilon here
            grid_vx[i, j] /= grid_m[i, j]  # Momentum to velocity
            grid_vy[i, j] /= grid_m[i, j]  # Momentum to velocity
            grid_vx[i, j] += dt*gravity_x*30  # gravity
            grid_vy[i, j] += dt*gravity_y*30  # gravity
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
    end
    return nothing
end

function G2P(x, v, C, n_particles, grid_vx, grid_vy, dt, inv_dx)
    p = (blockIdx().x - 1)*blockDim().x + threadIdx().x
    if p <= n_particles
        base = SVector{2, Int}(Int((x[p][1]*inv_dx - 0.5)÷1), Int((x[p][2]*inv_dx - 0.5)÷1))
        fx = x[p]*inv_dx - base
        wx = SVector{3}(0.5*(1.5 - fx[1])^2, 0.75 - (fx[1] - 1)^2, 0.5*(fx[1] - 0.5)^2)
        wy = SVector{3}(0.5*(1.5 - fx[2])^2, 0.75 - (fx[2] - 1)^2, 0.5*(fx[2] - 0.5)^2)
        new_v = SVector{2}(0.0, 0.0)
        new_C = SMatrix{2,2}(0.0, 0.0, 0.0, 0.0)
        for i in 0:2
            for j in 0:2
                index = base + SVector{2, Int}(i + 1, j + 1)
                dpos = SVector{2}(i, j) - fx
                g_v = SVector{2}(grid_vx[index[1], index[2]], grid_vy[index[1], index[2]])
                weight = wx[i + 1] * wy[j + 1]
                new_v += weight * g_v
                new_C += 4*inv_dx*weight*SMatrix{2,2}(g_v[1]*dpos[1], g_v[2]*dpos[1], g_v[1]*dpos[2], g_v[2]*dpos[2])
            end
        end
        v[p] = new_v
        C[p] = new_C
        x[p] += dt*new_v  # advection
    end
    return nothing
end

@inline function svd(A)
    x = A[1, 1] + A[2, 2]
    y = A[2, 1] - A[1, 2]
    scale = 1/sqrt(x*x + y*y)
    c = x*scale
    s = y*scale
    R = SMatrix{2,2}(c, s, -s, c)
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
        V = SMatrix{2,2}(-s, -c, c, -s);
    else
        V = SMatrix{2,2}(c, -s, s, c);
    end
    U = R*V;
    sig = SMatrix{2,2}(s1, 0.0, 0.0, s2);
    return U, sig, V
end

function test(quality)
    n_particles = 3*(4*quality)^2
    n_grid      = 9*quality

    dx          = 1.0/n_grid
    inv_dx      = float(n_grid)
    dt          = 1e-4/quality
    p_vol       = (dx*0.5)^2
    p_rho       = 1
    p_mass      = p_vol*p_rho
    E           = 5e3
    nu          = 0.2  # Young's modulus and Poisson's ratio
    mu_0        = E/(2*(1 + nu))
    lambda_0    = E*nu/((1 + nu)*(1 - 2*nu))  # Lame parameters

    x = CuArray{SVector{2, Float64}}(undef, n_particles) # position
    v = CuArray{SVector{2, Float64}}(undef, n_particles)  # velocity
    F = CuArray{SMatrix{2,2, Float64, 4}}(undef, n_particles)
    C = CuArray{SMatrix{2,2, Float64, 4}}(undef, n_particles) # affine velocity field
    Jp = CUDA.zeros(Float64, n_particles)  # plastic deformation
    material = CUDA.ones(Int64, n_particles)  # material

    grid_vx = CUDA.zeros(Float64, n_grid, n_grid) # position along x
    grid_vy = CUDA.zeros(Float64, n_grid, n_grid) # position along y
    grid_m = CUDA.zeros(Float64, n_grid, n_grid) # grid node mass

    gravity_x = 0.0
    gravity_y = -1.0

    @krun n_particles reset(x, v, F, Jp, C, n_particles, quality)

    for i in 1:5000
        @krun n_grid^2 reset_grid(grid_vx, grid_vy, grid_m, n_grid)
    
        @krun n_particles P2G(x, v, F, Jp, C, material, n_particles, grid_vx, grid_vy, grid_m, dt, dx, inv_dx, mu_0, lambda_0, p_vol, p_mass)

        @krun n_grid^2 grid_update(grid_vx, grid_vy, grid_m, n_grid, gravity_x, gravity_y, dt)

        @krun n_particles G2P(x, v, C, n_particles, grid_vx, grid_vy, dt, inv_dx)
    end
end

function time_test()
    for quality in 1:40
        @show 3*(4*quality)^2, 9*quality
        test(quality)
        @time test(quality)
    end
end

time_test()

println("Done")
