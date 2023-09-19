module fftUtils

"""module that contains utility functions for Fourier analysis."""

using Plots
using LaTeXStrings
using Printf

export n_mode_num, m_mode_num, n_mode_ind, m_mode_ind, q_prime

function n_mode_num(n,Nn)
    """Takes index of the mode array (n) and total number of n-modes (Nn) and returns the n-mode number."""
    if n == 1
        mode = 0
    end
    if n <= ((Nn-1)/2)+1 && n > 1
        mode =  n-1
    end
    if n > ((Nn-1)/2)+1
        mode = n-(Nn+1)
    end
    return mode
end

function m_mode_num(m,Nm)
    """Takes index of the mode array (m) and total number of m-modes (Nm)  and returns the m-mode number."""
    if m == 1
        mode = 0
    end
    if m <= ((Nm-1)/2)+1 && m > 1
        mode = -(m-1) # m's have flipped frequencies
    end
    if m > ((Nm-1)/2)+1
        mode = Nm-m+1
    end
    return mode
end

function m_mode_ind(m,Nm)
    """Takes m-mode number and returns it's index location."""
    if m == 0
        ind = 1
    end
    if m < 0 
        ind = abs(m)+1
    end
    if m > 0
        ind = Nm - m +1
    end
    return ind
end

function n_mode_ind(n,Nn)
    """Takes n-mode number and returns it's index location."""
    if n == 0
        ind = 1
    end
    if n > 0
        ind = n+1
    end
    if n < 0
        ind = Nn - (abs(n)-1)
    end
    return ind
end

function q_prime(q,dpsi)
    """Calculates dq/dpsi from the data."""
    L = length(q);
    qp = circshift(q,1);
    qm = circshift(q,-1);
    qprime = (qp-qm)/(2*dpsi);
    qprime[1] = (q[2]-q[1])/dpsi;
    qprime[end] = (q[end]-q[end-1])/dpsi
    return qprime
end

function n_m_max(F,q,t)
    qPrime = q_prime(q,1/100);
    maxbnm = [];
    isl_width = [];
    ms = [];
    misl = [];
    nisl = [];
    ns = [];
    islands = []
    Nm = size(F)[2];
    Nn = size(F)[3];
    for m in 1:size(F)[2]
        for n in 1:size(F)[3]
            max, maxInd = findmax(abs.(F[:,m,n]))
            append!(maxbnm,max)
            q_factor = (-1 ./q[maxInd,t].^3) .*qPrime[maxInd]
            append!(ns, n_mode_num(n,Nn))
            append!(ms, m_mode_num(m,Nm))
            if m !=1
                append!(nisl,n_mode_num(n,Nn))
                append!(misl,m_mode_num(m,Nm))
                wmn = sqrt((4*max)/abs(m_mode_num(m,Nm)*q_factor))
                append!(isl_width,wmn)
                push!(islands,(m_mode_num(m,Nm),n_mode_num(n,Nn),wmn))
            else
                append!(nisl,n_mode_num(n,Nn))
                append!(misl,m_mode_num(m,Nm))
                wmn = 0
                append!(isl_width,wmn)
                push!(islands,(m_mode_num(m,Nm),n_mode_num(n,Nn),wmn))
                
            end
        end
    end
    return ns,ms,maxbnm,nisl,misl,isl_width,islands
end
    
end
