# I am hard-coding some polynomials of Sp_{2n}-character varieties found in
#
# \bib{Cambo17}{thesis}{
#    author={Camb\`{o}, V.},
#    title={On the $E$-polynomial of parabolic $\Sp_{2n}$-character varieties},
#    year={2017},
#    type={Ph.D.\ thesis},
#    organization={Scuola Internazionale Superiore di Studi Avanzati (SISSA)},
#    note={\href{https://hdl.handle.net/20.500.11767/57152}{SISSA Digital Library}}
#    }
#

q = Pol(:q)
Phi1 = q-1
Phi2 = q+1
Phi4 = q^2+1




# E-polynomial when G=Sp4 and n=1
function ESP4(g::Int)
    return Phi1^(4g-4)*(
        (8-2^(2g+2))*(q^4)^(2g-1)
        +(2^(2g+1)-8)*(q^4*Phi2)^(2g-1)
        +(2^(2g+1)-8)*(q^3*Phi2)^(2g-1)
        +(2^(2g)+2)*(q^3*Phi2^2)^(2g-1)
        +(2^(2g))*(q^2*Phi2^2)^(2g-1)
        +(2^(2g))*(q^4*Phi2^2)^(2g-1)
        +(2^(2g))*(q^3*Phi4)^(2g-1)
        +(1)*(Phi4*Phi2^2)^(2g-1)
        +(1)*(q^4*Phi4*Phi2^2)^(2g-1)
        )
end
