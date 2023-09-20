module DivClean

#Module that contains functions that perform divergence cleanup on magnetic field components.
export quad_r_int, quad_u_int, Au, Aphi, grid_Au, grid_Aphi, b1_an, b2_an, b3_an

using QuadGK


function quad_r_int(B_eint,xmin,xmax,yo)
    res, err = quadgk(x -> B_eint(x,yo),xmin,xmax)
    return res
end

function quad_r_int(B_eint,xmin,xmax,yo,zo)
    res, err = quadgk(x -> B_eint(x,yo,zo),xmin,xmax)
    return res
end

function quad_u_int(B_eint,ymin,ymax,xo,zo)
    res, err = quadgk(y -> B_eint(xo,y,zo),ymin,ymax)
    return res
end

function quad_u_int(B_eint,ymin,ymax,xo)
    res, err = quadgk(y -> B_eint(xo,y),ymin,ymax)
    return res
end

function Au(B3_eint,ro,uo,rMA,uMA)
    A = quad_r_int(B3_eint,rMA,ro,uo)
    return A
end

function Au(B3_eint,ro,uo,rMA,uMA,zo)
    A = quad_r_int(B3_eint,rMA,ro,uo,zo)
    return A
end

function Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA)
    A = quad_u_int(B1_eint,uMA,uo,ro) - quad_r_int(B2_eint,rMA,ro,uMA)
    return A
end

function Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,zo)
    A = quad_u_int(B1_eint,uMA,uo,ro,zo) - quad_r_int(B2_eint,rMA,ro,uMA,zo)
    return A
end

function grid_Au(B3_eint,rMA,uMA,rn,un,phin)
    Au_arr = []
    for ro in rn
        for uo in un
            for fo in phin
                append!(Au_arr,Au(B3_eint,ro,uo,rMA,uMA,fo))
            end
        end
    end
    Au_arr = permutedims(reshape(Au_arr, size(phin)[1], size(un)[1],size(rn)[1]),(3,2,1))
    return Au_arr
end

function grid_Au(B3_eint,rMA,uMA,rn,un)    
    Au_arr = []
    for ro in rn
        for uo in un
            append!(Au_arr,Au(B3_eint,ro,uo,rMA,uMA))
        end
    end
    Au_arr = permutedims(reshape(Au_arr, size(un)[1],size(rn)[1]),(2,1))
    return Au_arr
end

function grid_Aphi(B1_eint,B2_eint,rMA,uMA,rn,un,phin)
    Aphi_arr = []
    for ro in rn
        for uo in un
            for fo in phin
                append!(Aphi_arr,Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA,fo))
            end
        end
    end
    Aphi_arr = permutedims(reshape(Aphi_arr, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    return Aphi_arr
end

function grid_Aphi(B1_eint,B2_eint,rMA,uMA,rn,un)
    Aphi_arr = []
    for ro in rn
        for uo in un
            append!(Aphi_arr,Aphi(B1_eint,B2_eint,ro,uo,rMA,uMA))
        end
    end
    Aphi_arr = permutedims(reshape(Aphi_arr, size(un)[1],size(rn)[1]),(2,1))
    return Aphi_arr
end

function b1_an(GAp_eint,rn,un)
    b1_an = []
    for r in rn
        for u in un
            append!(b1_an, Interpolations.gradient(GAp_eint,r,u)[2])
        end
    end
    b1_an = permutedims(reshape(b1_an, size(un)[1],size(rn)[1]),(2,1))
    return b1_an
end

function b1_an(GAp_eint,rn,un,phin)
    b1_an = []
    for r in rn
        for u in un
            for f in phin
                append!(b1_an, Interpolations.gradient(GAp_eint,r,u,f)[2])
            end
        end
    end
    b1_an = permutedims(reshape(b1_an, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    return b1_an
end
    
function b2_an(GAp_eint,rn,un)
    b2_an = []
    for r in rn
        for u in un
            append!(b2_an, -Interpolations.gradient(GAp_eint,r,u)[1])
        end
    end
    b2_an = permutedims(reshape(b2_an, size(un)[1],size(rn)[1]),(2,1))
    return b2_an
end

function b2_an(GAp_eint,rn,un,phin)
    b2_an = []
    for r in rn
        for u in un
            for f in phin
                append!(b2_an, -Interpolations.gradient(GAp_eint,r,u,f)[1])
            end
        end
    end
    b2_an = permutedims(reshape(b2_an, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    return b2_an
end 

function b3_an(GAu_eint,rn,un)
    b3_an = []
    for r in rn
        for u in un
            append!(b3_an, Interpolations.gradient(GAu_eint,r,u)[1])
        end
    end
    b3_an = permutedims(reshape(b3_an, size(un)[1],size(rn)[1]),(2,1))
    return b3_an
end

function b3_an(GAu_eint,rn,un,phin)
    b3_an = []
    for r in rn
        for u in un
            for f in phin
                append!(b3_an, Interpolations.gradient(GAu_eint,r,u,f)[1])
            end
        end
    end
    b3_an = permutedims(reshape(b3_an, size(phin)[1],size(un)[1],size(rn)[1]),(3,2,1))
    return b3_an
end































end
